using System;
using System.Collections.Generic;
using System.Text;

using static System.Math;
using static GeographicLib.MathEx;
using static GeographicLib.Macros;
using System.Diagnostics;

namespace GeographicLib
{
    /// <summary>
    /// Implements geodesic calculations.
    /// </summary>
    /// <remarks>
    /// The shortest path between two points on a ellipsoid at (<i>lat1</i>, <i>lon1</i>)
    /// and (<i>lat2</i>, <i>lon2</i>) is called the geodesic.Its length is <i>s12</i> and
    /// the geodesic from point 1 to point 2 has azimuths <i>azi1</i> and <i>azi2</i> at
    /// the two end points.  (The azimuth is the heading measured clockwise from
    /// north.  <i>azi2</i> is the "forward" azimuth, i.e., the heading that takes you
    /// beyond point 2 not back to point 1.)
    /// <para>
    /// Given <i>lat1</i>, <i>lon1</i>, <i>azi1</i>, and <i>s12</i>, we can determine <i>lat2</i>,
    /// <i>lon2</i>, and <i>azi2</i>.This is the <b>direct</b> geodesic problem and its
    /// solution is given by the function <see cref="IGeodesic.Direct(double, double, double, double, out double, out double)"/>.  (If <i>s12</i> is
    /// sufficiently large that the geodesic wraps more than halfway around the
    /// earth, there will be another geodesic between the points with a smaller <i>s12</i>.)
    /// </para>
    /// <para>
    /// Given <i>lat1</i>, <i>lon1</i>, <i>lat2</i>, and <i>lon2</i>, we can determine <i>azi1</i>,
    /// <i>azi2</i>, and <i>s12</i>.This is the <b>inverse</b> geodesic problem, whose solution
    /// is given by <see cref="IGeodesic.Inverse(double, double, double, double, out double)"/>. Usually, the solution to the inverse
    /// problem is unique.In cases where there are multiple solutions (all with
    /// the same <i>s12</i>, of course), all the solutions can be easily generated
    /// once a particular solution is provided.
    /// </para>
    /// <para>
    /// The standard way of specifying the direct problem is the specify the
    /// distance <i>s12</i> to the second point.However it is sometimes useful
    /// instead to specify the arc length <i>a12</i> (in degrees) on the auxiliary
    /// sphere.This is a mathematical construct used in solving the geodesic
    /// problems.  The solution of the direct problem in this form is provided by
    /// <see cref="IGeodesic.ArcDirect(double, double, double, double, out double, out double, out double)"/>. An arc length in excess of 180° indicates that
    /// the geodesic is not a shortest path.In addition, the arc length between
    /// an equatorial crossing and the next extremum of latitude for a geodesic is
    /// 90°.
    /// </para>
    /// <para>This class can also calculate several other quantities related to geodesics. These are:
    /// <list type="bullet">
    /// <item>
    /// <i>reduced length</i>. If we fix the first point and increase <i>azi1</i> by <i>dazi1</i> (radians), 
    /// the second point is displaced <i>m12</i> <i>dazi1</i> in the direction <i>azi2</i> + 90°. 
    /// The quantity <i>m12</i> is called the "reduced length" and is symmetric under interchange of the two points.
    /// On a curved surface the reduced length obeys a symmetry relation, <i>m12</i> + <i>m21</i> = 0. 
    /// On a flat surface, we have <i>m12</i> = <i>s12</i>. The ratio <i>s12</i>/<i>m12</i> gives the azimuthal scale 
    /// for an azimuthal equidistant projection.
    /// </item>
    /// <item>
    /// <i>geodesic scale</i>. Consider a reference geodesic and a second geodesic parallel to this one at point 1 
    /// and separated by a small distance <i>dt</i>. The separation of the two geodesics at point 2 is <i>M12</i> <i>dt</i>
    /// where <i>M12</i> is called the "geodesic scale". <i>M21</i> is defined similarly (with the geodesics being parallel at point 2).
    /// On a flat surface, we have <i>M12</i> = <i>M21</i> = 1. The quantity 1/<i>M12</i> gives the scale of the Cassini-Soldner projection.
    /// </item>
    /// <item>
    /// <i>area</i>. The area between the geodesic from point 1 to point 2 and the equation is represented by <i>S12</i>;
    /// it is the area, measured counter-clockwise, of the geodesic quadrilateral with corners (<i>lat1</i>, <i>lon1</i>), 
    /// (0,<i>lon1</i>), (0,<i>lon2</i>), and (<i>lat2</i>, <i>lon2</i>). It can be used to compute the area of any geodesic polygon.
    /// </item>
    /// </list>
    /// </para>
    /// <para>
    /// Overloaded versions of <see cref="IGeodesic.Direct(double, double, double, double, out double, out double)"/>,
    /// <see cref="IGeodesic.ArcDirect(double, double, double, double, out double, out double, out double)"/>,
    /// and <see cref="IGeodesic.Inverse(double, double, double, double, out double)"/> allow these quantities to be returned.
    /// In addition there are general functions <see cref="GenDirect"/>, 
    /// and <see cref="GenInverse(double, double, double, double, GeodesicFlags, out double, out double, out double, out double, out double, out double, out double)"/> which allow an arbitrary set of results
    /// to be computed. The quantities <i>m12</i>, <i>M12</i>, <i>M21</i> which all specify the behavior of nearby geodesics obey addition rules.
    /// If points 1, 2, and 3 all lie on a single geodesic, then the following rules hold:
    /// <list type="bullet">
    /// <item><i>s13</i> = <i>s12</i> + <i>s23</i></item>
    /// <item><i>a13</i> = <i>a12</i> + <i>a23</i></item>
    /// <item><i>S13</i> = <i>S12</i> + <i>S23</i></item>
    /// <item><i>m13</i> = <i>m12</i> <i>M23</i> + <i>m23</i> <i>M21</i></item>
    /// <item><i>M13</i> = <i>M12</i> <i>M23</i> − (1 − <i>M12</i> <i>M2</i>1) <i>m23</i> / <i>m12</i></item>
    /// <item><i>M31</i> = <i>M32</i> <i>M21</i> − (1 − <i>M23</i> <i>M32</i>) <i>m12</i> / <i>m23</i></item>
    /// </list>
    /// </para>
    /// <para>
    /// Additional functionality is provided by the <see cref="GeodesicLine"/> class, which allows a sequence of points along a geodesic to be computed.
    /// </para>
    /// <para>
    /// The shortest distance returned by the solution of the inverse problem is (obviously) uniquely defined. However, in a few special cases there are
    /// multiple azimuths which yield the same shortest distance. Here is a catalog of those cases:
    /// <list type="bullet">
    /// <item>
    /// <i>lat1</i> = <i>−lat2</i> (with neither point at a pole). If <i>azi1</i> = <i>azi2</i>, the geodesic is unique. 
    /// Otherwise there are two geodesics and the second one is obtained by setting 
    /// [<i>azi1</i>, <i>azi2</i>] → [<i>azi2</i>, <i>azi1</i>], [<i>M12</i>, <i>M21</i>] → [<i>M21</i>, <i>M12</i>], <i>S12</i> → −<i>S12</i>. 
    /// (This occurs when the longitude difference is near ±180° for oblate ellipsoids.)
    /// </item>
    /// <item>
    /// <i>lon2</i> = <i>lon1</i> ± 180° (with neither point at a pole). If <i>azi1</i> = 0° or ±180°, the geodesic is unique.
    /// Otherwise there are two geodesics and the second one is obtained by setting 
    /// [<i>azi1</i>, <i>azi2</i>] → [−<i>azi1</i>, −<i>azi2</i>], <i>S12</i> → −<i>S12</i>. 
    /// (This occurs when <i>lat2</i> is near −<i>lat1</i>  for prolate ellipsoids.)
    /// </item>
    /// <item>
    /// Points 1 and 2 at opposite poles. 
    /// There are infinitely many geodesics which can be generated by setting 
    /// [<i>azi1</i>, <i>azi2</i>] → [<i>azi1</i>, <i>azi2</i>] + [<i>d</i>, -<i>d</i>], for arbitrary <i>d</i>. 
    /// (For spheres, this prescription applies when points 1 and 2 are antipodal.)</item>
    /// <item>
    /// <i>s12</i> = 0 (coincident points). There are infinitely many geodesics which can be generated by setting
    /// [<i>azi1</i>, <i>azi2</i>] → [<i>azi1</i>, <i>azi2</i>] + [<i>d</i>, <i>d</i>], for arbitrary <i>d</i>.</item>
    /// </list>
    /// </para>
    /// <para>
    /// The calculations are accurate to better than 15 nm (15 nanometers) for the WGS84 ellipsoid.
    /// See Sec. 9 of <a href="https://arxiv.org/abs/1102.1215v1">arXiv:1102.1215v1</a> for details.
    /// The algorithms used by this class are based on series expansions using the flattening f as a small parameter.
    /// These are only accurate for |<i>f</i>| &lt; 0.02; however reasonably accurate results will be obtained for |<i>f</i>| &lt; 0.2.
    /// Here is a table of the approximate maximum error (expressed as a distance) for an ellipsoid with the same equatorial 
    /// radius as the WGS84 ellipsoid and different values of the flattening.
    /// <list type="table">
    /// <listheader>
    /// |<i>f</i>| error
    /// </listheader>
    /// <item>0.01 25 nm</item>
    /// <item>0.02 30 nm</item>
    /// <item>0.05 10 um</item>
    /// <item>0.1 1.5 mm</item>
    /// <item>0.2 300 mm</item>
    /// </list>
    /// </para>
    /// <para>For very eccentric ellipsoids, use <see cref="GeodesicExact"/> instead.</para>
    /// <para>The algorithms are described in
    /// <list type="bullet">
    /// <item>
    /// C. F. F. Karney, <a href="https://doi.org/10.1007/s00190-012-0578-z">Algorithms for geodesics</a>,
    /// J. Geodesy 87, 43–55 (2013); DOI: <a href="https://doi.org/10.1007/s00190-012-0578-z">10.1007/s00190-012-0578-z</a>; 
    /// addenda: <a href="https://geographiclib.sourceforge.io/geod-addenda.html">geod-addenda.html</a>.
    /// </item>
    /// </list>
    /// </para>
    /// <para>
    /// For more information on geodesics see <a href="https://geographiclib.sourceforge.io/html/geodesic.html">Geodesics on an ellipsoid of revolution</a>.
    /// </para>
    /// </remarks>
    public class Geodesic : GeodesicBase
    {
        internal const int nA1_ = GEOGRAPHICLIB_GEODESIC_ORDER,
                          nC1_ = GEOGRAPHICLIB_GEODESIC_ORDER,
                          nC1p_ = GEOGRAPHICLIB_GEODESIC_ORDER,
                          nA2_ = GEOGRAPHICLIB_GEODESIC_ORDER,
                          nC2_ = GEOGRAPHICLIB_GEODESIC_ORDER,
                          nA3_ = GEOGRAPHICLIB_GEODESIC_ORDER,
                          nA3x_ = nA3_,
                          nC3_ = GEOGRAPHICLIB_GEODESIC_ORDER,
                          nC3x_ = (nC3_ * (nC3_ - 1)) / 2,
                          nC4_ = GEOGRAPHICLIB_GEODESIC_ORDER,
                          nC4x_ = (nC4_ * (nC4_ + 1)) / 2,
                          nC_ = GEOGRAPHICLIB_GEODESIC_ORDER + 1,
                          maxit1_ = 20;

        private readonly int maxit2_;
        internal readonly double tiny_, tol0_, tol1_, tol2_, tolb_, xthresh_;

        internal readonly double _a, _f, _f1, _e2, _ep2, _n, _b, _c2, _etol2;

        private readonly Memory<double>
            _A3x = new double[nA3x_],
            _C3x = new double[nC3x_],
            _C4x = new double[nC4x_];

        /// <summary>
        /// A global instantiation of <see cref="Geodesic"/> with the parameters for the WGS84 ellipsoid.
        /// </summary>
        public static Geodesic WGS84 { get; } = new Geodesic(Ellipsoid.WGS84);

        /// <inheritdoc/>
        public override double EquatorialRadius => _a;

        /// <inheritdoc/>
        public override double Flattening => _f;

        /// <inheritdoc/>
        public override double EllipsoidArea => 4 * PI * _c2;

        /// <summary>
        /// Constructor for a ellipsoid with equatorial radius and its flattening.
        /// </summary>
        /// <param name="a">equatorial radius (meters).</param>
        /// <param name="f">flattening of ellipsoid.  Setting <i>f</i> = 0 gives a sphere.</param>
        public Geodesic(double a, double f)
        {
            maxit2_ = maxit1_ + DBL_MANT_DIG + 10;

            // Underflow guard.  We require
            //   tiny_ * epsilon() > 0
            //   tiny_ + epsilon() == epsilon()
            tiny_ = Sqrt(DBL_MIN);
            tol0_ = DBL_EPSILON;

            // Increase multiplier in defn of tol1_ from 100 to 200 to fix inverse
            // case 52.784459512564 0 -52.784459512563990912 179.634407464943777557
            // which otherwise failed for Visual Studio 10 (Release and Debug)
            tol1_ = 200 * tol0_;
            tol2_ = Sqrt(tol0_);
            tolb_ = tol0_ * tol2_; // Check on bisection interval
            xthresh_ = 1000 * tol2_;
            _a = a;
            _f = f;
            _f1 = 1 - _f;
            _e2 = _f * (2 - _f);
            _ep2 = _e2 / Sq(_f1); // e2 / (1 - e2)
            _n = _f / (2 - _f);
            _b = _a * _f1;
            _c2 = (Sq(_a) + Sq(_b) *
                   (_e2 == 0 ? 1 :
                    EAtanhE(1, (_f < 0 ? -1 : 1) * Sqrt(Abs(_e2))) / _e2))
                  / 2; // authalic radius squared

            // The sig12 threshold for "really short".  Using the auxiliary sphere
            // solution with dnm computed at (bet1 + bet2) / 2, the relative error in
            // the azimuth consistency check is sig12^2 * abs(f) * min(1, 1-f/2) / 2.
            // (Error measured for 1/100 < b/a < 100 and abs(f) >= 1/1000.  For a
            // given f and sig12, the max error occurs for lines near the pole.  If
            // the old rule for computing dnm = (dn1 + dn2)/2 is used, then the error
            // increases by a factor of 2.)  Setting this equal to epsilon gives
            // sig12 = etol2.  Here 0.1 is a safety factor (error decreased by 100)
            // and max(0.001, abs(f)) stops etol2 getting too large in the nearly
            // spherical case.
            _etol2 = 0.1 * tol2_ /
                Sqrt(Max(0.001, Abs(_f)) * Min(1, 1 - _f / 2) / 2);

            if (!(IsFinite(_a) && _a > 0))
                throw new GeographicException("Equatorial radius is not positive");
            if (!(IsFinite(_b) && _b > 0))
                throw new GeographicException("Polar semi-axis is not positive");

            A3coeff();
            C3coeff();
            C4coeff();
        }

        /// <summary>
        /// Constructor for a ellipsoid with equatorial radius and its flattening copied from another <see cref="IEllipsoid"/> object.
        /// </summary>
        /// <param name="ellipsoid">Source <see cref="IEllipsoid"/> object.</param>
        public Geodesic(IEllipsoid ellipsoid)
            : this(ellipsoid.EquatorialRadius, ellipsoid.Flattening)
        {

        }

        #region Internal methods

        internal static double SinCosSeries(bool sinp, double sinx, double cosx, ReadOnlyMemory<double> c, int n)
            => SinCosSeries(sinp, sinx, cosx, c.Span, n);

        internal static double SinCosSeries(bool sinp, double sinx, double cosx, ReadOnlySpan<double> c, int n)
        {
            // Evaluate
            // y = sinp ? sum(c[i] * sin( 2*i    * x), i, 1, n) :
            //            sum(c[i] * cos((2*i+1) * x), i, 0, n-1)
            // using Clenshaw summation.  N.B. c[0] is unused for sin series
            // Approx operation count = (n + 5) mult and (2 * n + 2) add
            var idx = n + (sinp ? 1 : 0);            // Point to one beyond last element
            double
              ar = 2 * (cosx - sinx) * (cosx + sinx), // 2 * cos(2 * x)
              y0 = (n & 1) == 1 ? c[--idx] : 0, y1 = 0;          // accumulators for sum
                                                                 // Now n is even
            n /= 2;
            while (n-- > 0)
            {
                // Unroll loop x 2, so accumulators return to their original role
                y1 = ar * y0 - y1 + c[--idx];
                y0 = ar * y1 - y0 + c[--idx];
            }
            return sinp
              ? 2 * sinx * cosx * y0    // sin(2 * x) * y0
              : cosx * (y0 - y1);       // cos(x) * (y0 - y1)
        }

        // The static const coefficient arrays in the following functions are
        // generated by Maxima and give the coefficients of the Taylor expansions for
        // the geodesics.  The convention on the order of these coefficients is as
        // follows:
        //
        //   ascending order in the trigonometric expansion,
        //   then powers of eps in descending order,
        //   finally powers of n in descending order.
        //
        // (For some expansions, only a subset of levels occur.)  For each polynomial
        // of order n at the lowest level, the (n+1) coefficients of the polynomial
        // are followed by a divisor which is applied to the whole polynomial.  In
        // this way, the coefficients are expressible with no round off error.  The
        // sizes of the coefficient arrays are:
        //
        //   A1m1f, A2m1f            = floor(N/2) + 2
        //   C1f, C1pf, C2f, A3coeff = (N^2 + 7*N - 2*floor(N/2)) / 4
        //   C3coeff       = (N - 1) * (N^2 + 7*N - 2*floor(N/2)) / 8
        //   C4coeff       = N * (N + 1) * (N + 5) / 6
        //
        // where N = GEOGRAPHICLIB_GEODESIC_ORDER
        //         = nA1 = nA2 = nC1 = nC1p = nA3 = nC4

        /// <summary>
        /// The scale factor A1-1 = mean value of (d/dsigma)I1 - 1
        /// </summary>
        /// <param name="eps"></param>
        /// <returns></returns>
        internal static double A1m1f(double eps)
        {
            if (GEOGRAPHICLIB_GEODESIC_ORDER / 2 != 3)
            {
                throw new GeographicException("Bad value for GEOGRAPHICLIB_GEODESIC_ORDER");
            }

            Span<double> coeff = stackalloc double[] {
                // (1-eps)*A1-1, polynomial in eps2 of order 3
                1,
                4,
                64,
                0,
                256,
            };

            Debug.Assert(
                coeff.Length == nA1_ / 2 + 2,
                "Coefficient array size mismatch in A1m1f"
            );

            var m = nA1_ / 2;
            var t = PolyVal(m, coeff, Sq(eps)) / coeff[m + 1];
            return (t + eps) / (1 - eps);
        }

        /// <summary>
        /// The coefficients C1[l] in the Fourier expansion of B1
        /// </summary>
        /// <param name="eps"></param>
        /// <param name="c"></param>
        internal static void C1f(double eps, Memory<double> c) => C1f(eps, c.Span);

        /// <summary>
        /// The coefficients C1[l] in the Fourier expansion of B1
        /// </summary>
        /// <param name="eps"></param>
        /// <param name="c"></param>
        internal static void C1f(double eps, Span<double> c)
        {
            if (GEOGRAPHICLIB_GEODESIC_ORDER != 6)
            {
                throw new GeographicException("Bad value for GEOGRAPHICLIB_GEODESIC_ORDER");
            }

            Span<double> coeff = stackalloc double[]
            {
                // C1[1]/eps^1, polynomial in eps2 of order 2
                -1,
                6,
                -16,
                32,
                // C1[2]/eps^2, polynomial in eps2 of order 2
                -9,
                64,
                -128,
                2048,
                // C1[3]/eps^3, polynomial in eps2 of order 1
                9,
                -16,
                768,
                // C1[4]/eps^4, polynomial in eps2 of order 1
                3,
                -5,
                512,
                // C1[5]/eps^5, polynomial in eps2 of order 0
                -7,
                1280,
                // C1[6]/eps^6, polynomial in eps2 of order 0
                -7,
                2048,
            };

            Debug.Assert(
                coeff.Length == (nC1_ * nC1_ + 7 * nC1_ - 2 * (nC1_ / 2)) / 4,
                "Coefficient array size mismatch in C1f"
            );

            double
              eps2 = Sq(eps),
              d = eps;

            int o = 0;
            for (int l = 1; l <= nC1_; ++l)
            { // l is index of C1p[l]
                int m = (nC1_ - l) / 2;         // order of polynomial in eps^2
                c[l] = d * PolyVal(m, coeff.Slice(o), eps2) / coeff[o + m + 1];
                o += m + 2;
                d *= eps;
            }
            // Post condition: o == sizeof(coeff) / sizeof(real)
        }

        internal static void C1pf(double eps, Memory<double> c) => C1pf(eps, c.Span);

        internal static void C1pf(double eps, Span<double> c)
        {
            if (GEOGRAPHICLIB_GEODESIC_ORDER != 6)
            {
                throw new GeographicException("Bad value for GEOGRAPHICLIB_GEODESIC_ORDER");
            }

            Span<double> coeff = stackalloc double[]
            {
                // C1p[1]/eps^1, polynomial in eps2 of order 2
                205,
                -432,
                768,
                1536,
                // C1p[2]/eps^2, polynomial in eps2 of order 2
                4005,
                -4736,
                3840,
                12288,
                // C1p[3]/eps^3, polynomial in eps2 of order 1
                -225,
                116,
                384,
                // C1p[4]/eps^4, polynomial in eps2 of order 1
                -7173,
                2695,
                7680,
                // C1p[5]/eps^5, polynomial in eps2 of order 0
                3467,
                7680,
                // C1p[6]/eps^6, polynomial in eps2 of order 0
                38081,
                61440,
            };

            Debug.Assert(
                coeff.Length == (nC1p_ * nC1p_ + 7 * nC1p_ - 2 * (nC1p_ / 2)) / 4,
                "Coefficient array size mismatch in C1pf"
            );

            double
              eps2 = Sq(eps),
              d = eps;
            int o = 0;
            for (int l = 1; l <= nC1p_; ++l)
            { // l is index of C1p[l]
                int m = (nC1p_ - l) / 2;         // order of polynomial in eps^2
                c[l] = d * PolyVal(m, coeff.Slice(o), eps2) / coeff[o + m + 1];
                o += m + 2;
                d *= eps;
            }
            // Post condition: o == sizeof(coeff) / sizeof(real)
        }

        /// <summary>
        /// The scale factor A2-1 = mean value of (d/dsigma)I2 - 1
        /// </summary>
        /// <param name="eps"></param>
        /// <returns></returns>
        internal static double A2m1f(double eps)
        {
            if (GEOGRAPHICLIB_GEODESIC_ORDER / 2 != 3)
            {
                throw new GeographicException("Bad value for GEOGRAPHICLIB_GEODESIC_ORDER");
            }

            Span<double> coeff = stackalloc double[]
            {
                // (eps+1)*A2-1, polynomial in eps2 of order 3
                -11,
                -28,
                -192,
                0,
                256,
            };

            Debug.Assert(
                coeff.Length == nA2_ / 2 + 2,
                "Coefficient array size mismatch in A2m1f"
            );

            int m = nA2_ / 2;
            var t = PolyVal(m, coeff, Sq(eps)) / coeff[m + 1];
            return (t - eps) / (1 + eps);
        }

        /// <summary>
        /// The coefficients C2[l] in the Fourier expansion of B2
        /// </summary>
        /// <param name="eps"></param>
        /// <param name="c"></param>
        internal static void C2f(double eps, Memory<double> c) => C2f(eps, c.Span);

        internal static void C2f(double eps, Span<double> c)
        {
            if (GEOGRAPHICLIB_GEODESIC_ORDER != 6)
            {
                throw new GeographicException("Bad value for GEOGRAPHICLIB_GEODESIC_ORDER");
            }

            Span<double> coeff = stackalloc double[]
            {
                // C2[1]/eps^1, polynomial in eps2 of order 2
                1,
                2,
                16,
                32,
                // C2[2]/eps^2, polynomial in eps2 of order 2
                35,
                64,
                384,
                2048,
                // C2[3]/eps^3, polynomial in eps2 of order 1
                15,
                80,
                768,
                // C2[4]/eps^4, polynomial in eps2 of order 1
                7,
                35,
                512,
                // C2[5]/eps^5, polynomial in eps2 of order 0
                63,
                1280,
                // C2[6]/eps^6, polynomial in eps2 of order 0
                77,
                2048,
            };

            Debug.Assert(
                coeff.Length == (nC2_ * nC2_ + 7 * nC2_ - 2 * (nC2_ / 2)) / 4,
                "Coefficient array size mismatch in C2f"
            );

            double
              eps2 = Sq(eps),
              d = eps;
            int o = 0;
            for (int l = 1; l <= nC2_; ++l)
            { // l is index of C2[l]
                int m = (nC2_ - l) / 2;         // order of polynomial in eps^2
                c[l] = d * PolyVal(m, coeff.Slice(o), eps2) / coeff[o + m + 1];
                o += m + 2;
                d *= eps;
            }
            // Post condition: o == sizeof(coeff) / sizeof(real)
        }

        internal double A3f(double eps) => PolyVal(nA3_ - 1, _A3x, eps);

        internal void C3f(double eps, Memory<double> c) => C3f(eps, c.Span);

        internal void C3f(double eps, Span<double> c)
        {
            // Evaluate C3 coeffs
            // Elements c[1] thru c[nC3_ - 1] are set
            var mult = 1d;
            var o = 0;
            for (int l = 1; l < nC3_; ++l)
            { // l is index of C3[l]
                int m = nC3_ - l - 1;          // order of polynomial in eps
                mult *= eps;
                c[l] = mult * PolyVal(m, _C3x.Slice(o), eps);
                o += m + 1;
            }
            // Post condition: o == nC3x_
        }

        internal void C4f(double eps, Memory<double> c) => C4f(eps, c.Span);

        internal void C4f(double eps, Span<double> c)
        {
            // Evaluate C4 coeffs
            // Elements c[0] thru c[nC4_ - 1] are set
            var mult = 1d;
            int o = 0;
            for (int l = 0; l < nC4_; ++l)
            { // l is index of C4[l]
                int m = nC4_ - l - 1;          // order of polynomial in eps
                c[l] = mult * PolyVal(m, _C4x.Slice(o), eps);
                o += m + 1;
                mult *= eps;
            }
            // Post condition: o == nC4x_
        }

        #endregion

        #region Private methods

        /// <summary>
        /// The scale factor A3 = mean value of (d/dsigma)I3
        /// </summary>
        private void A3coeff()
        {
            if (GEOGRAPHICLIB_GEODESIC_ORDER != 6)
            {
                throw new GeographicException("Bad value for GEOGRAPHICLIB_GEODESIC_ORDER");
            }

            Span<double> coeff = stackalloc double[]
            {
                // A3, coeff of eps^5, polynomial in n of order 0
                -3,
                128,
                // A3, coeff of eps^4, polynomial in n of order 1
                -2,
                -3,
                64,
                // A3, coeff of eps^3, polynomial in n of order 2
                -1,
                -3,
                -1,
                16,
                // A3, coeff of eps^2, polynomial in n of order 2
                3,
                -1,
                -2,
                8,
                // A3, coeff of eps^1, polynomial in n of order 1
                1,
                -1,
                2,
                // A3, coeff of eps^0, polynomial in n of order 0
                1,
                1,
            };

            Debug.Assert(
                coeff.Length == (nA3_ * nA3_ + 7 * nA3_ - 2 * (nA3_ / 2)) / 4,
                "Coefficient array size mismatch in A3coeff"
            );

            int o = 0, k = 0;
            for (int j = nA3_ - 1; j >= 0; --j)
            { // coeff of eps^j
                int m = Min(nA3_ - j - 1, j);       // order of polynomial in n
                _A3x.Span[k++] = PolyVal(m, coeff.Slice(o), _n) / coeff[o + m + 1];
                o += m + 2;
            }
            // Post condition: o == sizeof(coeff) / sizeof(real) && k == nA3x_
        }

        /// <summary>
        /// The coefficients C3[l] in the Fourier expansion of B3
        /// </summary>
        private void C3coeff()
        {
            if (GEOGRAPHICLIB_GEODESIC_ORDER != 6)
            {
                throw new GeographicException("Bad value for GEOGRAPHICLIB_GEODESIC_ORDER");
            }

            Span<double> coeff = stackalloc double[]
            {
                // C3[1], coeff of eps^5, polynomial in n of order 0
                3,
                128,
                // C3[1], coeff of eps^4, polynomial in n of order 1
                2,
                5,
                128,
                // C3[1], coeff of eps^3, polynomial in n of order 2
                -1,
                3,
                3,
                64,
                // C3[1], coeff of eps^2, polynomial in n of order 2
                -1,
                0,
                1,
                8,
                // C3[1], coeff of eps^1, polynomial in n of order 1
                -1,
                1,
                4,
                // C3[2], coeff of eps^5, polynomial in n of order 0
                5,
                256,
                // C3[2], coeff of eps^4, polynomial in n of order 1
                1,
                3,
                128,
                // C3[2], coeff of eps^3, polynomial in n of order 2
                -3,
                -2,
                3,
                64,
                // C3[2], coeff of eps^2, polynomial in n of order 2
                1,
                -3,
                2,
                32,
                // C3[3], coeff of eps^5, polynomial in n of order 0
                7,
                512,
                // C3[3], coeff of eps^4, polynomial in n of order 1
                -10,
                9,
                384,
                // C3[3], coeff of eps^3, polynomial in n of order 2
                5,
                -9,
                5,
                192,
                // C3[4], coeff of eps^5, polynomial in n of order 0
                7,
                512,
                // C3[4], coeff of eps^4, polynomial in n of order 1
                -14,
                7,
                512,
                // C3[5], coeff of eps^5, polynomial in n of order 0
                21,
                2560,
            };

            Debug.Assert(
                coeff.Length == ((nC3_ - 1) * (nC3_ * nC3_ + 7 * nC3_ - 2 * (nC3_ / 2))) / 8,
                "Coefficient array size mismatch in C3coeff"
            );

            int o = 0, k = 0;
            for (int l = 1; l < nC3_; ++l)
            {        // l is index of C3[l]
                for (int j = nC3_ - 1; j >= l; --j)
                { // coeff of eps^j
                    int m = Min(nC3_ - j - 1, j);       // order of polynomial in n
                    _C3x.Span[k++] = PolyVal(m, coeff.Slice(o), _n) / coeff[o + m + 1];
                    o += m + 2;
                }
            }
            // Post condition: o == sizeof(coeff) / sizeof(real) && k == nC3x_
        }

        private void C4coeff()
        {
            if (GEOGRAPHICLIB_GEODESIC_ORDER != 6)
            {
                throw new GeographicException("Bad value for GEOGRAPHICLIB_GEODESIC_ORDER");
            }

            Span<double> coeff = stackalloc double[]
            {
                // C4[0], coeff of eps^5, polynomial in n of order 0
                97,
                15015,
                // C4[0], coeff of eps^4, polynomial in n of order 1
                1088,
                156,
                45045,
                // C4[0], coeff of eps^3, polynomial in n of order 2
                -224,
                -4784,
                1573,
                45045,
                // C4[0], coeff of eps^2, polynomial in n of order 3
                -10656,
                14144,
                -4576,
                -858,
                45045,
                // C4[0], coeff of eps^1, polynomial in n of order 4
                64,
                624,
                -4576,
                6864,
                -3003,
                15015,
                // C4[0], coeff of eps^0, polynomial in n of order 5
                100,
                208,
                572,
                3432,
                -12012,
                30030,
                45045,
                // C4[1], coeff of eps^5, polynomial in n of order 0
                1,
                9009,
                // C4[1], coeff of eps^4, polynomial in n of order 1
                -2944,
                468,
                135135,
                // C4[1], coeff of eps^3, polynomial in n of order 2
                5792,
                1040,
                -1287,
                135135,
                // C4[1], coeff of eps^2, polynomial in n of order 3
                5952,
                -11648,
                9152,
                -2574,
                135135,
                // C4[1], coeff of eps^1, polynomial in n of order 4
                -64,
                -624,
                4576,
                -6864,
                3003,
                135135,
                // C4[2], coeff of eps^5, polynomial in n of order 0
                8,
                10725,
                // C4[2], coeff of eps^4, polynomial in n of order 1
                1856,
                -936,
                225225,
                // C4[2], coeff of eps^3, polynomial in n of order 2
                -8448,
                4992,
                -1144,
                225225,
                // C4[2], coeff of eps^2, polynomial in n of order 3
                -1440,
                4160,
                -4576,
                1716,
                225225,
                // C4[3], coeff of eps^5, polynomial in n of order 0
                -136,
                63063,
                // C4[3], coeff of eps^4, polynomial in n of order 1
                1024,
                -208,
                105105,
                // C4[3], coeff of eps^3, polynomial in n of order 2
                3584,
                -3328,
                1144,
                315315,
                // C4[4], coeff of eps^5, polynomial in n of order 0
                -128,
                135135,
                // C4[4], coeff of eps^4, polynomial in n of order 1
                -2560,
                832,
                405405,
                // C4[5], coeff of eps^5, polynomial in n of order 0
                128,
                99099,
            };

            Debug.Assert(
                coeff.Length == (nC4_ * (nC4_ + 1) * (nC4_ + 5)) / 6,
                "Coefficient array size mismatch in C4coeff"
            );

            int o = 0, k = 0;
            for (int l = 0; l < nC4_; ++l)
            {        // l is index of C4[l]
                for (int j = nC4_ - 1; j >= l; --j)
                { // coeff of eps^j
                    int m = nC4_ - j - 1;               // order of polynomial in n
                    _C4x.Span[k++] = PolyVal(m, coeff.Slice(o), _n) / coeff[o + m + 1];
                    o += m + 2;
                }
            }
            // Post condition: o == sizeof(coeff) / sizeof(real) && k == nC4x_
        }

        private double Astroid(double x, double y)
        {
            // Solve k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0 for positive root k.
            // This solution is adapted from Geocentric::Reverse.
            double k;
            double
              p = Sq(x),
              q = Sq(y),
              r = (p + q - 1) / 6;
            if (!(q == 0 && r <= 0))
            {
                double
                  // Avoid possible division by zero when r = 0 by multiplying equations
                  // for s and t by r^3 and r, resp.
                  S = p * q / 4,            // S = r^3 * s
                  r2 = Sq(r),
                  r3 = r * r2,
                  // The discriminant of the quadratic equation for T3.  This is zero on
                  // the evolute curve p^(1/3)+q^(1/3) = 1
                  disc = S * (S + 2 * r3);
                var u = r;
                if (disc >= 0)
                {
                    var T3 = S + r3;
                    // Pick the sign on the sqrt to maximize abs(T3).  This minimizes loss
                    // of precision due to cancellation.  The result is unchanged because
                    // of the way the T is used in definition of u.
                    T3 += T3 < 0 ? -Sqrt(disc) : Sqrt(disc); // T3 = (r * t)^3
                                                             // N.B. cbrt always returns the real root.  cbrt(-8) = -2.
                    var T = Cbrt(T3); // T = r * t
                                      // T can be zero; but then r2 / T -> 0.
                    u += T + (T != 0 ? r2 / T : 0);
                }
                else
                {
                    // T is complex, but the way u is defined the result is real.
                    var ang = Atan2(Sqrt(-disc), -(S + r3));
                    // There are three possible cube roots.  We choose the root which
                    // avoids cancellation.  Note that disc < 0 implies that r < 0.
                    u += 2 * r * Cos(ang / 3);
                }
                double
                  v = Sqrt(Sq(u) + q),    // guaranteed positive
                                          // Avoid loss of accuracy when u < 0.
                  uv = u < 0 ? q / (v - u) : u + v, // u+v, guaranteed positive
                  w = (uv - q) / (2 * v);           // positive?
                                                    // Rearrange expression for k to avoid loss of accuracy due to
                                                    // subtraction.  Division by 0 not possible because uv > 0, w >= 0.
                k = uv / (Sqrt(uv + Sq(w)) + w);   // guaranteed positive
            }
            else
            {               // q == 0 && r <= 0
                            // y = 0 with |x| <= 1.  Handle this case directly.
                            // for y small, positive root is k = abs(y)/sqrt(1-x^2)
                k = 0;
            }
            return k;
        }

        private void Lengths(double eps, double sig12,
                 double ssig1, double csig1, double dn1,
                 double ssig2, double csig2, double dn2,
                 double cbet1, double cbet2, GeodesicFlags outmask,
                 out double s12b, out double m12b, out double m0,
                 out double M12, out double M21, Span<double> Ca)
        {
            s12b = m12b = m0 = M12 = M21 = double.NaN;

            // Return m12b = (reduced length)/_b; also calculate s12b = distance/_b,
            // and m0 = coefficient of secular term in expression for reduced length.

            outmask = outmask.Flags();
            // outmask & DISTANCE: set s12b
            // outmask & REDUCEDLENGTH: set m12b & m0
            // outmask & GEODESICSCALE: set M12 & M21

            double m0x = 0, J12 = 0, A1 = 0, A2 = 0;
            Span<double> Cb = stackalloc double[nC2_ + 1];

            if (outmask.HasAny(GeodesicFlags.Distance | GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale))
            {
                A1 = A1m1f(eps);
                C1f(eps, Ca);
                if (outmask.HasAny(GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale))
                {
                    A2 = A2m1f(eps);
                    C2f(eps, Cb);
                    m0x = A1 - A2;
                    A2 = 1 + A2;
                }
                A1 = 1 + A1;
            }
            if (outmask.HasAny(GeodesicFlags.Distance))
            {
                var B1 = SinCosSeries(true, ssig2, csig2, Ca, nC1_) -
                  SinCosSeries(true, ssig1, csig1, Ca, nC1_);
                // Missing a factor of _b
                s12b = A1 * (sig12 + B1);
                if (outmask.HasAny(GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale))
                {
                    var B2 = SinCosSeries(true, ssig2, csig2, Cb, nC2_) -
                      SinCosSeries(true, ssig1, csig1, Cb, nC2_);
                    J12 = m0x * sig12 + (A1 * B1 - A2 * B2);
                }
            }
            else if (outmask.HasAny(GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale))
            {
                // Assume here that nC1_ >= nC2_
                for (int l = 1; l <= nC2_; ++l)
                    Cb[l] = A1 * Ca[l] - A2 * Cb[l];
                J12 = m0x * sig12 + (SinCosSeries(true, ssig2, csig2, Cb, nC2_) -
                                     SinCosSeries(true, ssig1, csig1, Cb, nC2_));
            }
            if (outmask.HasAny(GeodesicFlags.ReducedLength))
            {
                m0 = m0x;
                // Missing a factor of _b.
                // Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure
                // accurate cancellation in the case of coincident points.
                m12b = dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) -
                  csig1 * csig2 * J12;
            }
            if (outmask.HasAny(GeodesicFlags.GeodesicScale))
            {
                var csig12 = csig1 * csig2 + ssig1 * ssig2;
                var t = _ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2);
                M12 = csig12 + (t * ssig2 - csig2 * J12) * ssig1 / dn1;
                M21 = csig12 - (t * ssig1 - csig1 * J12) * ssig2 / dn2;
            }
        }

        private double InverseStart(double sbet1, double cbet1, double dn1,
                  double sbet2, double cbet2, double dn2,
                  double lam12, double slam12, double clam12,
                  out double salp1, out double calp1,
                  out double salp2, out double calp2, out double dnm,
                  Span<double> Ca)

        {
            salp2 = calp2 = dnm = double.NaN;

            // Return a starting point for Newton's method in salp1 and calp1 (function
            // value is -1).  If Newton's method doesn't need to be used, return also
            // salp2 and calp2 and function value is sig12.
            double
              sig12 = -1,               // Return value
                                        // bet12 = bet2 - bet1 in [0, pi); bet12a = bet2 + bet1 in (-pi, 0]
              sbet12 = sbet2 * cbet1 - cbet2 * sbet1,
              cbet12 = cbet2 * cbet1 + sbet2 * sbet1;
            var sbet12a = sbet2 * cbet1 + cbet2 * sbet1;
            bool shortline = cbet12 >= 0 && sbet12 < 0.5 &&
              cbet2 * lam12 < 0.5;
            double somg12, comg12;
            if (shortline)
            {
                var sbetm2 = Sq(sbet1 + sbet2);
                // sin((bet1+bet2)/2)^2
                // =  (sbet1 + sbet2)^2 / ((sbet1 + sbet2)^2 + (cbet1 + cbet2)^2)
                sbetm2 /= sbetm2 + Sq(cbet1 + cbet2);
                dnm = Sqrt(1 + _ep2 * sbetm2);
                var omg12 = lam12 / (_f1 * dnm);
                somg12 = Sin(omg12); comg12 = Cos(omg12);
            }
            else
            {
                somg12 = slam12; comg12 = clam12;
            }

            salp1 = cbet2 * somg12;
            calp1 = comg12 >= 0 ?
              sbet12 + cbet2 * sbet1 * Sq(somg12) / (1 + comg12) :
              sbet12a - cbet2 * sbet1 * Sq(somg12) / (1 - comg12);

            double
              ssig12 = Hypot(salp1, calp1),
              csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12;

            if (shortline && ssig12 < _etol2)
            {
                // really short lines
                salp2 = cbet1 * somg12;
                calp2 = sbet12 - cbet1 * sbet2 *
                  (comg12 >= 0 ? Sq(somg12) / (1 + comg12) : 1 - comg12);
                Norm(ref salp2, ref calp2);
                // Set return value
                sig12 = Atan2(ssig12, csig12);
            }
            else if (Abs(_n) > 0.1 || // Skip astroid calc if too eccentric
                     csig12 >= 0 ||
                     ssig12 >= 6 * Abs(_n) * PI * Sq(cbet1))
            {
                // Nothing to do, zeroth order spherical approximation is OK
            }
            else
            {
                // Scale lam12 and bet2 to x, y coordinate system where antipodal point
                // is at origin and singular point is at y = 0, x = -1.
                double x, y, lamscale, betscale;
                double lam12x = Atan2(-slam12, -clam12); // lam12 - pi
                if (_f >= 0)
                {            // In fact f == 0 does not get here
                             // x = dlong, y = dlat
                    {
                        double
                          k2 = Sq(sbet1) * _ep2,
                          eps = k2 / (2 * (1 + Sqrt(1 + k2)) + k2);
                        lamscale = _f * cbet1 * A3f(eps) * PI;
                    }
                    betscale = lamscale * cbet1;

                    x = lam12x / lamscale;
                    y = sbet12a / betscale;
                }
                else
                {                  // _f < 0
                                   // x = dlat, y = dlong
                    double
                      cbet12a = cbet2 * cbet1 - sbet2 * sbet1,
                      bet12a = Atan2(sbet12a, cbet12a);

                    // In the case of lon12 = 180, this repeats a calculation made in
                    // Inverse.
                    Lengths(_n, PI + bet12a,
                            sbet1, -cbet1, dn1, sbet2, cbet2, dn2,
                            cbet1, cbet2,
                            GeodesicFlags.ReducedLength, out _, out var m12b, out var m0, out _, out _, Ca);
                    x = -1 + m12b / (cbet1 * cbet2 * m0 * PI);
                    betscale = x < -0.01 ? sbet12a / x :
                      -_f * Sq(cbet1) * PI;
                    lamscale = betscale / cbet1;
                    y = lam12x / lamscale;
                }

                if (y > -tol1_ && x > -1 - xthresh_)
                {
                    // strip near cut
                    // Need real(x) here to cast away the volatility of x for min/max
                    if (_f >= 0)
                    {
                        salp1 = Min(1, -x); calp1 = -Sqrt(1 - Sq(salp1));
                    }
                    else
                    {
                        calp1 = Max(x > -tol1_ ? 0 : -1, x);
                        salp1 = Sqrt(1 - Sq(calp1));
                    }
                }
                else
                {
                    // Estimate alp1, by solving the astroid problem.
                    //
                    // Could estimate alpha1 = theta + pi/2, directly, i.e.,
                    //   calp1 = y/k; salp1 = -x/(1+k);  for _f >= 0
                    //   calp1 = x/(1+k); salp1 = -y/k;  for _f < 0 (need to check)
                    //
                    // However, it's better to estimate omg12 from astroid and use
                    // spherical formula to compute alp1.  This reduces the mean number of
                    // Newton iterations for astroid cases from 2.24 (min 0, max 6) to 2.12
                    // (min 0 max 5).  The changes in the number of iterations are as
                    // follows:
                    //
                    // change percent
                    //    1       5
                    //    0      78
                    //   -1      16
                    //   -2       0.6
                    //   -3       0.04
                    //   -4       0.002
                    //
                    // The histogram of iterations is (m = number of iterations estimating
                    // alp1 directly, n = number of iterations estimating via omg12, total
                    // number of trials = 148605):
                    //
                    //  iter    m      n
                    //    0   148    186
                    //    1 13046  13845
                    //    2 93315 102225
                    //    3 36189  32341
                    //    4  5396      7
                    //    5   455      1
                    //    6    56      0
                    //
                    // Because omg12 is near pi, estimate work with omg12a = pi - omg12
                    var k = Astroid(x, y);
                    var
                      omg12a = lamscale * (_f >= 0 ? -x * k / (1 + k) : -y * (1 + k) / k);
                    somg12 = Sin(omg12a); comg12 = -Cos(omg12a);
                    // Update spherical estimate of alp1 using omg12 instead of lam12
                    salp1 = cbet2 * somg12;
                    calp1 = sbet12a - cbet2 * sbet1 * Sq(somg12) / (1 - comg12);
                }
            }
            // Sanity check on starting guess.  Backwards check allows NaN through.
            if (!(salp1 <= 0))
                Norm(ref salp1, ref calp1);
            else
            {
                salp1 = 1; calp1 = 0;
            }
            return sig12;
        }

        private double Lambda12(double sbet1, double cbet1, double dn1,
                  double sbet2, double cbet2, double dn2,
                  double salp1, double calp1, double slam120, double clam120,
                  out double salp2, out double calp2, out double sig12,
                  out double ssig1, out double csig1, out double ssig2, out double csig2,
                  out double eps, out double domg12,
                  bool diffp, out double dlam12, Span<double> Ca)
        {
            dlam12 = double.NaN;

            if (sbet1 == 0 && calp1 == 0)
                // Break degeneracy of equatorial line.  This case has already been
                // handled.
                calp1 = -tiny_;

            double
              // sin(alp1) * cos(bet1) = sin(alp0)
              salp0 = salp1 * cbet1,
              calp0 = Hypot(calp1, salp1 * sbet1); // calp0 > 0

            double somg1, comg1, somg2, comg2, somg12, comg12, lam12;
            // tan(bet1) = tan(sig1) * cos(alp1)
            // tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1)
            ssig1 = sbet1; somg1 = salp0 * sbet1;
            csig1 = comg1 = calp1 * cbet1;
            Norm(ref ssig1, ref csig1);
            // Math::norm(somg1, comg1); -- don't need to normalize!

            // Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
            // about this case, since this can yield singularities in the Newton
            // iteration.
            // sin(alp2) * cos(bet2) = sin(alp0)
            salp2 = cbet2 != cbet1 ? salp0 / cbet2 : salp1;
            // calp2 = sqrt(1 - sq(salp2))
            //       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
            // and subst for calp0 and rearrange to give (choose positive sqrt
            // to give alp2 in [0, pi/2]).
            calp2 = cbet2 != cbet1 || Abs(sbet2) != -sbet1 ?
              Sqrt(Sq(calp1 * cbet1) +
                   (cbet1 < -sbet1 ?
                    (cbet2 - cbet1) * (cbet1 + cbet2) :
                    (sbet1 - sbet2) * (sbet1 + sbet2))) / cbet2 :
              Abs(calp1);
            // tan(bet2) = tan(sig2) * cos(alp2)
            // tan(omg2) = sin(alp0) * tan(sig2).
            ssig2 = sbet2; somg2 = salp0 * sbet2;
            csig2 = comg2 = calp2 * cbet2;
            Norm(ref ssig2, ref csig2);
            // Math::norm(somg2, comg2); -- don't need to normalize!

            // sig12 = sig2 - sig1, limit to [0, pi]
            sig12 = Atan2(Max(0, csig1 * ssig2 - ssig1 * csig2) + 0d,
                                       csig1 * csig2 + ssig1 * ssig2);

            // omg12 = omg2 - omg1, limit to [0, pi]
            somg12 = Max(0, comg1 * somg2 - somg1 * comg2) + 0d;
            comg12 = comg1 * comg2 + somg1 * somg2;
            // eta = omg12 - lam120
            var eta = Atan2(somg12 * clam120 - comg12 * slam120,
                             comg12 * clam120 + somg12 * slam120);
            double B312;
            var k2 = Sq(calp0) * _ep2;
            eps = k2 / (2 * (1 + Sqrt(1 + k2)) + k2);
            C3f(eps, Ca);
            B312 = (SinCosSeries(true, ssig2, csig2, Ca, nC3_ - 1) -
                    SinCosSeries(true, ssig1, csig1, Ca, nC3_ - 1));
            domg12 = -_f * A3f(eps) * salp0 * (sig12 + B312);
            lam12 = eta + domg12;

            if (diffp)
            {
                if (calp2 == 0)
                    dlam12 = -2 * _f1 * dn1 / sbet1;
                else
                {
                    Lengths(eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                            cbet1, cbet2, GeodesicFlags.ReducedLength,
                            out _, out dlam12, out _, out _, out _, Ca);
                    dlam12 *= _f1 / (calp2 * cbet2);
                }
            }

            return lam12;
        }

        private double GenInverse(double lat1, double lon1, double lat2, double lon2,
                    GeodesicFlags outmask, out double s12,
                    out double salp1, out double calp1, out double salp2, out double calp2,
                    out double m12, out double M12, out double M21, out double S12)
        {
            salp1 = salp2 = calp1 = calp2 = double.NaN;
            s12 = S12 = m12 = M12 = M21 = double.NaN;

            // Compute longitude difference (AngDiff does this carefully).  Result is
            // in [-180, 180] but -180 is only for west-going geodesics.  180 is for
            // east-going and meridional geodesics.
            var lon12 = AngDiff(lon1, lon2, out var lon12s);

            // Make longitude difference positive.
            int lonsign = SignBit(lon12) ? -1 : 1;
            // If very close to being on the same half-meridian, then make it so.
            lon12 = lonsign * AngRound(lon12);
            lon12s = AngRound((180 - lon12) - lonsign * lon12s);
            double
              lam12 = lon12 * Degree,
              slam12, clam12;

            if (lon12 > 90)
            {
                SinCosd(lon12s, out slam12, out clam12);
                clam12 = -clam12;
            }
            else
                SinCosd(lon12, out slam12, out clam12);

            // If really close to the equator, treat as on equator.
            lat1 = AngRound(LatFix(lat1));
            lat2 = AngRound(LatFix(lat2));

            // Swap points so that point with higher (abs) latitude is point 1.
            // If one latitude is a nan, then it becomes lat1.
            int swapp = Abs(lat1) < Abs(lat2) || double.IsNaN(lat2) ? -1 : 1;
            if (swapp < 0)
            {
                lonsign *= -1;
                Swap(ref lat1, ref lat2);
            }

            // Make lat1 <= 0
            int latsign = SignBit(lat1) ? 1 : -1;
            lat1 *= latsign;
            lat2 *= latsign;

            // Now we have
            //
            //     0 <= lon12 <= 180
            //     -90 <= lat1 <= 0
            //     lat1 <= lat2 <= -lat1
            //
            // longsign, swapp, latsign register the transformation to bring the
            // coordinates to this canonical form.  In all cases, 1 means no change was
            // made.  We make these transformations so that there are few cases to
            // check, e.g., on verifying quadrants in atan2.  In addition, this
            // enforces some symmetries in the results returned.

            double sbet1, cbet1, sbet2, cbet2, s12x = double.NaN, m12x = double.NaN;

            SinCosd(lat1, out sbet1, out cbet1); sbet1 *= _f1;
            // Ensure cbet1 = +epsilon at poles; doing the fix on beta means that sig12
            // will be <= 2*tiny for two points at the same pole.
            Norm(ref sbet1, ref cbet1); cbet1 = Max(tiny_, cbet1);

            SinCosd(lat2, out sbet2, out cbet2); sbet2 *= _f1;
            // Ensure cbet2 = +epsilon at poles
            Norm(ref sbet2, ref cbet2); cbet2 = Max(tiny_, cbet2);

            // If cbet1 < -sbet1, then cbet2 - cbet1 is a sensitive measure of the
            // |bet1| - |bet2|.  Alternatively (cbet1 >= -sbet1), abs(sbet2) + sbet1 is
            // a better measure.  This logic is used in assigning calp2 in Lambda12.
            // Sometimes these quantities vanish and in that case we force bet2 = +/-
            // bet1 exactly.  An example where is is necessary is the inverse problem
            // 48.522876735459 0 -48.52287673545898293 179.599720456223079643
            // which failed with Visual Studio 10 (Release and Debug)

            if (cbet1 < -sbet1)
            {
                if (cbet2 == cbet1)
                    sbet2 = SignBit(sbet2) ? sbet1 : -sbet1;
            }
            else
            {
                if (Abs(sbet2) == -sbet1)
                    cbet2 = cbet1;
            }

            double
              dn1 = Sqrt(1 + _ep2 * Sq(sbet1)),
              dn2 = Sqrt(1 + _ep2 * Sq(sbet2));

            double a12 = double.NaN, sig12;
            // index zero element of this array is unused
            Span<double> Ca = stackalloc double[nC_];

            bool meridian = lat1 == -90 || slam12 == 0;

            if (meridian)
            {

                // Endpoints are on a single full meridian, so the geodesic might lie on
                // a meridian.

                calp1 = clam12; salp1 = slam12; // Head to the target longitude
                calp2 = 1; salp2 = 0;           // At the target we're heading north

                double
                  // tan(bet) = tan(sig) * cos(alp)
                  ssig1 = sbet1, csig1 = calp1 * cbet1,
                  ssig2 = sbet2, csig2 = calp2 * cbet2;

                // sig12 = sig2 - sig1
                sig12 = Atan2(Max(0, csig1 * ssig2 - ssig1 * csig2) + 0d,
                                           csig1 * csig2 + ssig1 * ssig2);
                {
                    Lengths(_n, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2,
                            outmask | GeodesicFlags.Distance | GeodesicFlags.ReducedLength,
                            out s12x, out m12x, out _, out M12, out M21, Ca);
                }

                // Add the check for sig12 since zero length geodesics might yield m12 <
                // 0.  Test case was
                //
                //    echo 20.001 0 20.001 0 | GeodSolve -i
                //
                // In fact, we will have sig12 > pi/2 for meridional geodesic which is
                // not a shortest path.
                // TODO: investigate m12 < 0 result for aarch/ppc (with -f -p 20)
                // 20.001000000000001 0.000000000000000 180.000000000000000
                // 20.001000000000001 0.000000000000000 180.000000000000000
                // 0.0000000002 0.000000000000001 -0.0000000001
                // 0.99999999999999989 0.99999999999999989 0.000
                if (sig12 < 1 || m12x >= 0)
                {
                    // Need at least 2, to handle 90 0 90 180
                    if (sig12 < 3 * tiny_ ||
                            // Prevent negative s12 or m12 for short lines
                            (sig12 < tol0_ && (s12x < 0 || m12x < 0)))
                        sig12 = m12x = s12x = 0;
                    m12x *= _b;
                    s12x *= _b;
                    a12 = sig12 / Degree;
                }
                else
                    // m12 < 0, i.e., prolate and too close to anti-podal
                    meridian = false;
            }

            // somg12 > 1 marks that it needs to be calculated
            double omg12 = 0, somg12 = 2, comg12 = 0;
            if (!meridian &&
                sbet1 == 0 &&   // and sbet2 == 0
                (_f <= 0 || lon12s >= _f * 180))
            {

                // Geodesic runs along equator
                calp1 = calp2 = 0; salp1 = salp2 = 1;
                s12x = _a * lam12;
                sig12 = omg12 = lam12 / _f1;
                m12x = _b * Sin(sig12);
                if (outmask.Flags().HasAny(GeodesicFlags.GeodesicScale))
                    M12 = M21 = Cos(sig12);
                a12 = lon12 / _f1;

            }
            else if (!meridian)
            {

                // Now point1 and point2 belong within a hemisphere bounded by a
                // meridian and geodesic is neither meridional or equatorial.

                // Figure a starting point for Newton's method
                sig12 = InverseStart(sbet1, cbet1, dn1, sbet2, cbet2, dn2,
                                     lam12, slam12, clam12,
                                     out salp1, out calp1, out salp2, out calp2, out var dnm,
                                     Ca);

                if (sig12 >= 0)
                {
                    // Short lines (InverseStart sets salp2, calp2, dnm)
                    s12x = sig12 * _b * dnm;
                    m12x = Sq(dnm) * _b * Sin(sig12 / dnm);
                    if (outmask.Flags().HasAny(GeodesicFlags.GeodesicScale))
                        M12 = M21 = Cos(sig12 / dnm);
                    a12 = sig12 / Degree;
                    omg12 = lam12 / (_f1 * dnm);
                }
                else
                {

                    // Newton's method.  This is a straightforward solution of f(alp1) =
                    // lambda12(alp1) - lam12 = 0 with one wrinkle.  f(alp) has exactly one
                    // root in the interval (0, pi) and its derivative is positive at the
                    // root.  Thus f(alp) is positive for alp > alp1 and negative for alp <
                    // alp1.  During the course of the iteration, a range (alp1a, alp1b) is
                    // maintained which brackets the root and with each evaluation of
                    // f(alp) the range is shrunk, if possible.  Newton's method is
                    // restarted whenever the derivative of f is negative (because the new
                    // value of alp1 is then further from the solution) or if the new
                    // estimate of alp1 lies outside (0,pi); in this case, the new starting
                    // guess is taken to be (alp1a + alp1b) / 2.
                    //
                    // initial values to suppress warnings (if loop is executed 0 times)
                    double ssig1 = 0, csig1 = 0, ssig2 = 0, csig2 = 0, eps = 0, domg12 = 0;
                    var numit = 0u;
                    // Bracketing range
                    double salp1a = tiny_, calp1a = 1, salp1b = tiny_, calp1b = -1;

                    for (bool tripn = false, tripb = false;
                         numit < maxit2_ || GEOGRAPHICLIB_PANIC;
                         ++numit)
                    {
                        // the WGS84 test set: mean = 1.47, sd = 1.25, max = 16
                        // WGS84 and random input: mean = 2.85, sd = 0.60
                        double dv;
                        var v = Lambda12(sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1,
                                          slam12, clam12,
                                          out salp2, out calp2, out sig12, out ssig1, out csig1, out ssig2, out csig2,
                                          out eps, out domg12, numit < maxit1_, out dv, Ca);
                        // Reversed test to allow escape with NaNs
                        if (tripb || !(Abs(v) >= (tripn ? 8 : 1) * tol0_)) break;
                        // Update bracketing values
                        if (v > 0 && (numit > maxit1_ || calp1 / salp1 > calp1b / salp1b))
                        { salp1b = salp1; calp1b = calp1; }
                        else if (v < 0 && (numit > maxit1_ || calp1 / salp1 < calp1a / salp1a))
                        { salp1a = salp1; calp1a = calp1; }
                        if (numit < maxit1_ && dv > 0)
                        {
                            var dalp1 = -v / dv;
                            double
                              sdalp1 = Sin(dalp1), cdalp1 = Cos(dalp1),
                              nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;

                            if (nsalp1 > 0 && Abs(dalp1) < PI)
                            {
                                calp1 = calp1 * cdalp1 - salp1 * sdalp1;
                                salp1 = nsalp1;
                                Norm(ref salp1, ref calp1);
                                // In some regimes we don't get quadratic convergence because
                                // slope -> 0.  So use convergence conditions based on epsilon
                                // instead of sqrt(epsilon).
                                tripn = Abs(v) <= 16 * tol0_;
                                continue;
                            }
                        }
                        // Either dv was not positive or updated value was outside legal
                        // range.  Use the midpoint of the bracket as the next estimate.
                        // This mechanism is not needed for the WGS84 ellipsoid, but it does
                        // catch problems with more eccentric ellipsoids.  Its efficacy is
                        // such for the WGS84 test set with the starting guess set to alp1 =
                        // 90deg:
                        // the WGS84 test set: mean = 5.21, sd = 3.93, max = 24
                        // WGS84 and random input: mean = 4.74, sd = 0.99
                        salp1 = (salp1a + salp1b) / 2;
                        calp1 = (calp1a + calp1b) / 2;
                        Norm(ref salp1, ref calp1);
                        tripn = false;
                        tripb = (Abs(salp1a - salp1) + (calp1a - calp1) < tolb_ ||
                                 Abs(salp1 - salp1b) + (calp1 - calp1b) < tolb_);
                    }
                    {
                        // Ensure that the reduced length and geodesic scale are computed in
                        // a "canonical" way, with the I2 integral.
                        var lengthmask = outmask |
                          (outmask.Flags().HasAny(GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale) ? GeodesicFlags.Distance : GeodesicFlags.None);

                        Lengths(eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                                cbet1, cbet2, lengthmask, out s12x, out m12x, out _, out M12, out M21, Ca);
                    }
                    m12x *= _b;
                    s12x *= _b;
                    a12 = sig12 / Degree;
                    if (outmask.Flags().HasAny(GeodesicFlags.Area))
                    {
                        // omg12 = lam12 - domg12
                        double sdomg12 = Sin(domg12), cdomg12 = Cos(domg12);
                        somg12 = slam12 * cdomg12 - clam12 * sdomg12;
                        comg12 = clam12 * cdomg12 + slam12 * sdomg12;
                    }
                }
            }

            if (outmask.Flags().HasAny(GeodesicFlags.Distance))
                s12 = 0d + s12x;           // Convert -0 to 0

            if (outmask.Flags().HasAny(GeodesicFlags.ReducedLength))
                m12 = 0d + m12x;           // Convert -0 to 0

            if (outmask.Flags().HasAny(GeodesicFlags.Area))
            {
                double
                  // From Lambda12: sin(alp1) * cos(bet1) = sin(alp0)
                  salp0 = salp1 * cbet1,
                  calp0 = Hypot(calp1, salp1 * sbet1); // calp0 > 0
                double alp12;
                if (calp0 != 0 && salp0 != 0)
                {
                    double
                      // From Lambda12: tan(bet) = tan(sig) * cos(alp)
                      ssig1 = sbet1, csig1 = calp1 * cbet1,
                      ssig2 = sbet2, csig2 = calp2 * cbet2,
                      k2 = Sq(calp0) * _ep2,
                      eps = k2 / (2 * (1 + Sqrt(1 + k2)) + k2),
                      // Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0).
                      A4 = Sq(_a) * calp0 * salp0 * _e2;
                    Norm(ref ssig1, ref csig1);
                    Norm(ref ssig2, ref csig2);
                    C4f(eps, Ca);
                    double
                      B41 = SinCosSeries(false, ssig1, csig1, Ca, nC4_),
                      B42 = SinCosSeries(false, ssig2, csig2, Ca, nC4_);
                    S12 = A4 * (B42 - B41);
                }
                else
                    // Avoid problems with indeterminate sig1, sig2 on equator
                    S12 = 0;

                if (!meridian && somg12 > 1)
                {
                    somg12 = Sin(omg12); comg12 = Cos(omg12);
                }

                if (!meridian &&
                    // omg12 < 3/4 * pi
                    comg12 > -0.7071 &&     // Long difference not too big
                    sbet2 - sbet1 < 1.75)
                { // Lat difference not too big
                  // Use tan(Gamma/2) = tan(omg12/2)
                  // * (tan(bet1/2)+tan(bet2/2))/(1+tan(bet1/2)*tan(bet2/2))
                  // with tan(x/2) = sin(x)/(1+cos(x))
                    double domg12 = 1 + comg12, dbet1 = 1 + cbet1, dbet2 = 1 + cbet2;
                    alp12 = 2 * Atan2(somg12 * (sbet1 * dbet2 + sbet2 * dbet1),
                                       domg12 * (sbet1 * sbet2 + dbet1 * dbet2));
                }
                else
                {
                    // alp12 = alp2 - alp1, used in atan2 so no need to normalize
                    double
                      salp12 = salp2 * calp1 - calp2 * salp1,
                      calp12 = calp2 * calp1 + salp2 * salp1;
                    // The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
                    // salp12 = -0 and alp12 = -180.  However this depends on the sign
                    // being attached to 0 correctly.  The following ensures the correct
                    // behavior.
                    if (salp12 == 0 && calp12 < 0)
                    {
                        salp12 = tiny_ * calp1;
                        calp12 = -1;
                    }
                    alp12 = Atan2(salp12, calp12);
                }
                S12 += _c2 * alp12;
                S12 *= swapp * lonsign * latsign;
                // Convert -0 to 0
                S12 += 0;
            }

            // Convert calp, salp to azimuth accounting for lonsign, swapp, latsign.
            if (swapp < 0)
            {
                Swap(ref salp1, ref salp2);
                Swap(ref calp1, ref calp2);
                if (outmask.Flags().HasAny(GeodesicFlags.GeodesicScale))
                    Swap(ref M12, ref M21);
            }

            salp1 *= swapp * lonsign; calp1 *= swapp * latsign;
            salp2 *= swapp * lonsign; calp2 *= swapp * latsign;

            // Returned value in [0, 180]
            return a12;
        }

        #endregion

        /// <inheritdoc/>
        public override double GenDirect(double lat1, double lon1, double azi1,
                                     bool arcmode, double s12_a12, GeodesicFlags outmask,
                                     out double lat2, out double lon2, out double azi2,
                                     out double s12, out double m12, out double M12, out double M21,
                                     out double S12)
        {
            // Automatically supply DISTANCE_IN if necessary
            if (!arcmode) outmask |= GeodesicFlags.DistanceIn;

            return new GeodesicLine(this, lat1, lon1, azi1, outmask)
              .                         // Note the dot!
              GenPosition(arcmode, s12_a12, outmask,
                          out lat2, out lon2, out azi2, out s12, out m12, out M12, out M21, out S12);
        }

        /// <inheritdoc/>
        public override double GenInverse(double lat1, double lon1, double lat2, double lon2,
                          GeodesicFlags outmask,
                          out double s12, out double azi1, out double azi2,
                          out double m12, out double M12, out double M21, out double S12)
        {
            outmask = outmask.Flags();

            double salp1, calp1, salp2, calp2,
              a12 = GenInverse(lat1, lon1, lat2, lon2,
                                outmask, out s12, out salp1, out calp1, out salp2, out calp2,
                                 out m12, out M12, out M21, out S12);
            if (outmask.HasAny(GeodesicFlags.Azimuth))
            {
                azi1 = Atan2d(salp1, calp1);
                azi2 = Atan2d(salp2, calp2);
            }
            else
            {
                azi1 = azi2 = double.NaN;
            }

            return a12;
        }

        /// <inheritdoc/>
        public override IGeodesicLine Line(double lat1, double lon1, double azi1, GeodesicFlags caps = GeodesicFlags.All)
            => new GeodesicLine(this, lat1, lon1, azi1, caps);

        /// <inheritdoc/>
        public override IGeodesicLine InverseLine(double lat1, double lon1, double lat2, double lon2,
            GeodesicFlags caps = GeodesicFlags.All)
        {
            double salp1, calp1,
                   a12 = GenInverse(lat1, lon1, lat2, lon2,
                                    // No need to specify AZIMUTH here
                                    0u, out _, out salp1, out calp1, out _, out _,
                                    out _, out _, out _, out _),
                   azi1 = Atan2d(salp1, calp1);

            // Ensure that a12 can be converted to a distance
            if (caps.Flags().HasAny(GeodesicFlags.DistanceIn)) caps |= GeodesicFlags.Distance;
            return new GeodesicLine(this, lat1, lon1, azi1, salp1, calp1, caps, true, a12);
        }

        /// <inheritdoc/>
        public override IGeodesicLine GenDirectLine(double lat1, double lon1, double azi1,
                           bool arcmode, double s12_a12,
                           GeodesicFlags caps = GeodesicFlags.All)
        {
            azi1 = AngNormalize(azi1);

            // Guard against underflow in salp0.  Also -0 is converted to +0.
            SinCosd(AngRound(azi1), out var salp1, out var calp1);
            // Automatically supply DISTANCE_IN if necessary
            if (!arcmode) caps |= GeodesicFlags.DistanceIn;
            return new GeodesicLine(this, lat1, lon1, azi1, salp1, calp1,
                                caps, arcmode, s12_a12);
        }
    }
}
