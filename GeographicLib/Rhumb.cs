using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using static System.Math;
using static GeographicLib.MathEx;
using static GeographicLib.Macros;
using System.Diagnostics;

namespace GeographicLib
{
    /// <summary>
    /// Solve of the direct and inverse rhumb problems.
    /// </summary>
    /// <remarks>
    /// The path of constant azimuth between two points on an ellipsoid at (<i>lat1</i>, <i>lon1</i>) and (<i>lat2</i>, <i>lon2</i>)
    /// is called the rhumb line (also called the loxodrome). Its length is <i>s12</i> and its azimuth is <i>azi12</i>.
    /// (The azimuth is the heading measured clockwise from north.)
    /// <para>
    /// Given <i>lat1</i>, <i>lon1</i>, <i>azi12</i>, and <i>s12</i>, we can determine <i>lat2</i>, and <i>lon2</i>.
    /// This is the direct rhumb problem and its solution is given by the function <see cref="Direct(double, double, double, double, out double, out double)"/>.
    /// </para>
    /// <para>
    /// Given <i>lat1</i>, <i>lon1</i>, <i>lat2</i>, and <i>lon2</i>, we can determine <i>azi12</i> and <i>s12</i>.
    /// This is the inverse rhumb problem, whose solution is given by <see cref="Inverse(double, double, double, double, out double, out double)"/>.
    /// This finds the shortest such rhumb line, i.e., the one that wraps no more than half way around the earth.
    /// If the end points are on opposite meridians, there are two shortest rhumb lines and the east-going one is chosen.
    /// </para>
    /// <para>
    /// These routines also optionally calculate the area under the rhumb line, <i>S12</i>.
    /// This is the area, measured counter-clockwise,
    /// of the rhumb line quadrilateral with corners (l<i>lat1</i>, <i>lon1</i>), (0,<i>lon1</i>), (0,<i>lon2</i>), and (<i>lat2</i>, <i>lon2</i>).
    /// </para>
    /// <para>
    /// Note that rhumb lines may be appreciably longer (up to 50%) than the corresponding <see cref="Geodesic"/>.
    /// For example the distance between London Heathrow and Tokyo Narita via the rhumb line is 11400 km which is 18% longer than the
    /// geodesic distance 9600 km.
    /// </para>
    /// <para>
    /// For more information on rhumb lines see <a href="https://geographiclib.sourceforge.io/html/rhumb.html">Rhumb lines</a>.
    /// </para>
    /// </remarks>
    public class Rhumb : IGeodesicLike
    {
        private const int tm_maxord = GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER;
        private const int maxpow_ = GEOGRAPHICLIB_RHUMBAREA_ORDER;

        internal readonly Ellipsoid _ell;
        private readonly bool _exact;
        internal readonly double _c2;

        // _R[0] unused
        private readonly Memory<double> _R = new double[maxpow_ + 1];

        /// <inheritdoc/>
        public double EquatorialRadius => _ell.EquatorialRadius;

        /// <inheritdoc/>
        public double Flattening => _ell.Flattening;

        /// <inheritdoc/>
        public double EllipsoidArea => _ell.Area;

        /// <summary>
        /// Gets a value representing whether current <see cref="Rhumb"/> instance performs exact calculation for arbitrary flattening.
        /// </summary>
        public bool IsExact => _exact;

        /// <summary>
        /// A global instantiation of <see cref="Rhumb"/> with the parameters for the WGS84 ellipsoid.
        /// </summary>
        public static Rhumb WGS84 { get; } = new Rhumb(Ellipsoid.WGS84, false);

        /// <summary>
        /// Initialize a new <see cref="Rhumb"/> instance with specified equatorial radius and flattening of the ellipsoid.
        /// </summary>
        /// <param name="a">equatorial radius (meters)</param>
        /// <param name="f">flattening of ellipsoid. Setting <i>f</i> = 0 gives a sphere. Negative <i>f</i> gives a prolate ellipsoid.</param>
        /// <param name="exact">
        /// if <see langword="true"/> (the default) use an addition theorem for elliptic integrals to compute divided differences;
        /// otherwise use series expansion (accurate for |<i>f</i>| &lt; 0.01).
        /// </param>
        /// <remarks>
        /// See <a href="https://geographiclib.sourceforge.io/html/rhumb.html">Rhumb lines</a>,
        /// for a detailed description of the <paramref name="exact"/> parameter.
        /// </remarks>
        public Rhumb(double a, double f, bool exact = true) : this(new Ellipsoid(a, f), exact) { }

        /// <summary>
        /// Initialize a new <see cref="Rhumb"/> instance with specified ellipsoid.
        /// </summary>
        /// <param name="ellipsoid">the ellipoid to be used for calculation.</param>
        /// <param name="exact">
        /// if <see langword="true"/> (the default) use an addition theorem for elliptic integrals to compute divided differences;
        /// otherwise use series expansion (accurate for |<i>f</i>| &lt; 0.01).
        /// </param>
        /// <remarks>
        /// See <a href="https://geographiclib.sourceforge.io/html/rhumb.html">Rhumb lines</a>,
        /// for a detailed description of the <paramref name="exact"/> parameter.
        /// </remarks>
        public Rhumb(IEllipsoid ellipsoid, bool exact = true) : this(ellipsoid.EquatorialRadius, ellipsoid.Flattening, exact) { }

        /// <summary>
        /// Initialize a new <see cref="Rhumb"/> instance with specified ellipsoid.
        /// </summary>
        /// <param name="ellipsoid">the ellipoid to be used for calculation.</param>
        /// <param name="exact">
        /// if <see langword="true"/> (the default) use an addition theorem for elliptic integrals to compute divided differences;
        /// otherwise use series expansion (accurate for |<i>f</i>| &lt; 0.01).
        /// </param>
        /// <remarks>
        /// See <a href="https://geographiclib.sourceforge.io/html/rhumb.html">Rhumb lines</a>,
        /// for a detailed description of the <paramref name="exact"/> parameter.
        /// </remarks>
        public Rhumb(Ellipsoid ellipsoid, bool exact = true)
        {
            Span<double> coeff = stackalloc[] {
                // R[0]/n^0, polynomial in n of order 6
                128346268d,
                -107884140,
                31126095,
                354053700,
                -908107200,
                851350500,
                0,
                2554051500,
                // R[1]/n^1, polynomial in n of order 5
                -114456994,
                56868630,
                79819740,
                -240540300,
                312161850,
                -212837625,
                638512875,
                // R[2]/n^2, polynomial in n of order 4
                51304574,
                24731070,
                -78693615,
                71621550,
                -28378350,
                212837625,
                // R[3]/n^3, polynomial in n of order 3
                1554472,
                -6282003,
                4684680,
                -1396395,
                14189175,
                // R[4]/n^4, polynomial in n of order 2
                -4913956,
                3205800,
                -791505,
                8108100,
                // R[5]/n^5, polynomial in n of order 1
                1092376,
                -234468,
                2027025,
                // R[6]/n^6, polynomial in n of order 0
                -313076,
                2027025,
            };

            Debug.Assert(GEOGRAPHICLIB_RHUMBAREA_ORDER == 6, "GEOGRAPHICLIB_RHUMBAREA_ORDER must be 6.");
            Debug.Assert(coeff.Length == ((maxpow_ + 1) * (maxpow_ + 4)) / 2, "Coefficient array size mismatch for Rhumb");

            _ell = ellipsoid;
            _exact = exact;
            _c2 = _ell.Area / (2 * TD);

            var d = 1.0;
            int o = 0;
            for (int l = 0; l <= maxpow_; ++l)
            {
                int m = maxpow_ - l;
                // R[0] is just an integration constant so it cancels when evaluating a
                // definite integral.  So don't bother computing it.  It won't be used
                // when invoking SinCosSeries.
                if (l != 0)
                    _R.Span[l] = d * PolyVal(m, coeff.Slice(o), _ell._n) / coeff[o + m + 1];
                o += m + 2;
                d *= _ell._n;
            }
            // Post condition: o == sizeof(alpcoeff) / sizeof(real)
        }

        #region Private Methods

        private static double Gd(double x) => Atan(Sinh(x));

        // Use divided differences to determine (mu2 - mu1) / (psi2 - psi1)
        // accurately
        //
        // Definition: Df(x,y,d) = (f(x) - f(y)) / (x - y)
        // See:
        //   W. M. Kahan and R. J. Fateman,
        //   Symbolic computation of divided differences,
        //   SIGSAM Bull. 33(2), 7-28 (1999)
        //   https://doi.org/10.1145/334714.334716
        //   http://www.cs.berkeley.edu/~fateman/papers/divdiff.pdf
        private static double Dlog(double x, double y)
        {
            var t = x - y;
            // Change
            //
            //   atanh(t / (x + y))
            //
            // to
            //
            //   asinh(t / (2 * sqrt(x*y)))
            //
            // to avoid taking atanh(1) when x is large and y is 1.  N.B., this
            // routine is invoked with positive x and y, so no need to guard against
            // taking the sqrt of a negative quantity.  This fixes bogus results for
            // the area being returning when an endpoint is at a pole.
            return t != 0 ? 2 * Asinh(t / (2 * Sqrt(x * y))) / t : 1 / x;
        }

        // N.B., x and y are in degrees
        private static double Dtan(double x, double y)
        {
            double d = x - y, tx = Tand(x), ty = Tand(y), txy = tx * ty;
            return d != 0 ?
              (2 * txy > -1 ? (1 + txy) * Tand(d) : tx - ty) /
              (d * Degree) :
              1 + txy;
        }

        private static double Datan(double x, double y)
        {
            double d = x - y, xy = x * y;
            return d != 0 ?
              (2 * xy > -1 ? Atan(d / (1 + xy)) : Atan(x) - Atan(y)) / d :
              1 / (1 + xy);
        }

        private static double Dsin(double x, double y)
        {
            double d = (x - y) / 2;
            return Cos((x + y) / 2) * (d != 0 ? Sin(d) / d : 1);
        }

        private static double Dsinh(double x, double y)
        {
            double d = (x - y) / 2;
            return Cosh((x + y) / 2) * (d != 0 ? Sinh(d) / d : 1);
        }

        private static double Dcosh(double x, double y)
        {
            double d = (x - y) / 2;
            return Sinh((x + y) / 2) * (d != 0 ? Sinh(d) / d : 1);
        }

        private static double Dasinh(double x, double y)
        {
            double d = x - y,
              hx = Hypot(1, x), hy = Hypot(1, y);
            return d != 0 ?
              Asinh(x * y > 0 ? d * (x + y) / (x * hy + y * hx) : x * hy - y * hx) / d :
              1 / hx;
        }

        private static double Dgd(double x, double y) => Datan(Sinh(x), Sinh(y)) * Dsinh(x, y);

        // N.B., x and y are the tangents of the angles
        private static double Dgdinv(double x, double y) => Dasinh(x, y) / Datan(x, y);

        // Copied from LambertConformalConic...
        // Deatanhe(x,y) = eatanhe((x-y)/(1-e^2*x*y))/(x-y)
        private double Deatanhe(double x, double y)
        {
            double t = x - y, d = 1 - _ell._e2 * x * y;
            return t != 0 ? EAtanhE(t / d, _ell._es) / t : _ell._e2 / d;
        }

        // (E(x) - E(y)) / (x - y) -- E = incomplete elliptic integral of 2nd kind
        private double DE(double x, double y)
        {
            var ei = _ell._ell;
            var d = x - y;
            if (x * y <= 0)
                return d != 0 ? (ei.E(x) - ei.E(y)) / d : 1;
            // See DLMF: Eqs (19.11.2) and (19.11.4) letting
            // theta -> x, phi -> -y, psi -> z
            //
            // (E(x) - E(y)) / d = E(z)/d - k2 * sin(x) * sin(y) * sin(z)/d
            //
            // tan(z/2) = (sin(x)*Delta(y) - sin(y)*Delta(x)) / (cos(x) + cos(y))
            //          = d * Dsin(x,y) * (sin(x) + sin(y))/(cos(x) + cos(y)) /
            //             (sin(x)*Delta(y) + sin(y)*Delta(x))
            //          = t = d * Dt
            // sin(z) = 2*t/(1+t^2); cos(z) = (1-t^2)/(1+t^2)
            // Alt (this only works for |z| <= pi/2 -- however, this conditions holds
            // if x*y > 0):
            // sin(z) = d * Dsin(x,y) * (sin(x) + sin(y))/
            //          (sin(x)*cos(y)*Delta(y) + sin(y)*cos(x)*Delta(x))
            // cos(z) = sqrt((1-sin(z))*(1+sin(z)))
            double sx = Sin(x), sy = Sin(y), cx = Cos(x), cy = Cos(y);
            double Dt = Dsin(x, y) * (sx + sy) /
              ((cx + cy) * (sx * ei.Delta(sy, cy) + sy * ei.Delta(sx, cx))),
              t = d * Dt, Dsz = 2 * Dt / (1 + t * t),
              sz = d * Dsz, cz = (1 - t) * (1 + t) / (1 + t * t);
            return ((sz != 0 ? ei.E(sz, cz, ei.Delta(sz, cz)) / sz : 1)
                    - ei.K2 * sx * sy) * Dsz;
        }

        // (mux - muy) / (phix - phiy) using elliptic integrals
        private double DRectifying(double latx, double laty)
        {
            double
              tbetx = _ell._f1 * Tand(latx),
              tbety = _ell._f1 * Tand(laty);

            return (PI / 2) * _ell._b * _ell._f1 * DE(Atan(tbetx), Atan(tbety))
              * Dtan(latx, laty) * Datan(tbetx, tbety) / _ell.QuarterMeridian;
        }

        // (psix - psiy) / (phix - phiy)
        private double DIsometric(double latx, double laty)
        {
            double
              phix = latx * Degree, tx = Tand(latx),
              phiy = laty * Degree, ty = Tand(laty);
            return Dasinh(tx, ty) * Dtan(latx, laty)
              - Deatanhe(Sin(phix), Sin(phiy)) * Dsin(phix, phiy);
        }

        // (sum(c[j]*sin(2*j*x),j=1..n) - sum(c[j]*sin(2*j*x),j=1..n)) / (x - y)
        private static double SinCosSeries(bool sinp,
                                 double x, double y, ReadOnlySpan<double> c, int n)
        {
            // N.B. n >= 0 and c[] has n+1 elements 0..n, of which c[0] is ignored.
            //
            // Use Clenshaw summation to evaluate
            //   m = (g(x) + g(y)) / 2         -- mean value
            //   s = (g(x) - g(y)) / (x - y)   -- average slope
            // where
            //   g(x) = sum(c[j]*SC(2*j*x), j = 1..n)
            //   SC = sinp ? sin : cos
            //   CS = sinp ? cos : sin
            //
            // This function returns only s; m is discarded.
            //
            // Write
            //   t = [m; s]
            //   t = sum(c[j] * f[j](x,y), j = 1..n)
            // where
            //   f[j](x,y) = [ (SC(2*j*x)+SC(2*j*y))/2 ]
            //               [ (SC(2*j*x)-SC(2*j*y))/d ]
            //
            //             = [       cos(j*d)*SC(j*p)    ]
            //               [ +/-(2/d)*sin(j*d)*CS(j*p) ]
            // (+/- = sinp ? + : -) and
            //    p = x+y, d = x-y
            //
            //   f[j+1](x,y) = A * f[j](x,y) - f[j-1](x,y)
            //
            //   A = [  2*cos(p)*cos(d)      -sin(p)*sin(d)*d]
            //       [ -4*sin(p)*sin(d)/d   2*cos(p)*cos(d)  ]
            //
            // Let b[n+1] = b[n+2] = [0 0; 0 0]
            //     b[j] = A * b[j+1] - b[j+2] + c[j] * I for j = n..1
            //    t =  (c[0] * I  - b[2]) * f[0](x,y) + b[1] * f[1](x,y)
            // c[0] is not accessed for s = t[2]
            double p = x + y, d = x - y,
              cp = Cos(p), cd = Cos(d),
              sp = Sin(p), sd = d != 0 ? Sin(d) / d : 1,
              m = 2 * cp * cd, s = sp * sd;
            // 2x2 matrices stored in row-major order
            Span<double>
                a = stackalloc[] { m, -s * d * d, -4 * s, m },
                ba = stackalloc[] { 0d, 0, 0, 0 },
                bb = stackalloc[] { 0d, 0, 0, 0 };

            var b1 = ba;
            var b2 = bb;
            if (n > 0) b1[0] = b1[3] = c[n];
            for (int j = n - 1; j > 0; --j)
            { // j = n-1 .. 1
                var _t = b1;
                b1 = b2;
                b2 = _t;

                // b1 = A * b2 - b1 + c[j] * I
                b1[0] = a[0] * b2[0] + a[1] * b2[2] - b1[0] + c[j];
                b1[1] = a[0] * b2[1] + a[1] * b2[3] - b1[1];
                b1[2] = a[2] * b2[0] + a[3] * b2[2] - b1[2];
                b1[3] = a[2] * b2[1] + a[3] * b2[3] - b1[3] + c[j];
            }
            // Here are the full expressions for m and s
            // m =   (c[0] - b2[0]) * f01 - b2[1] * f02 + b1[0] * f11 + b1[1] * f12;
            // s = - b2[2] * f01 + (c[0] - b2[3]) * f02 + b1[2] * f11 + b1[3] * f12;
            if (sinp)
            {
                // real f01 = 0, f02 = 0;
                double f11 = cd * sp, f12 = 2 * sd * cp;
                // m = b1[0] * f11 + b1[1] * f12;
                s = b1[2] * f11 + b1[3] * f12;
            }
            else
            {
                // real f01 = 1, f02 = 0;
                double f11 = cd * cp, f12 = -2 * sd * sp;
                // m = c[0] - b2[0] + b1[0] * f11 + b1[1] * f12;
                s = -b2[2] + b1[2] * f11 + b1[3] * f12;
            }
            return s;
        }

        // (mux - muy) / (chix - chiy) using Krueger's series
        private double DConformalToRectifying(double chix, double chiy)
            => 1 + SinCosSeries(true, chix, chiy, _ell.ConformalToRectifyingCoeffs.Span, tm_maxord);

        // (chix - chiy) / (mux - muy) using Krueger's series
        private double DRectifyingToConformal(double mux, double muy)
            => 1 - SinCosSeries(true, mux, muy, _ell.RectifyingToConformalCoeffs.Span, tm_maxord);

        // (mux - muy) / (psix - psiy)
        // N.B., psix and psiy are in degrees
        private double DIsometricToRectifying(double psix, double psiy)
        {
            if (_exact)
            {
                double
                  latx = _ell.InverseIsometricLatitude(psix),
                  laty = _ell.InverseIsometricLatitude(psiy);
                return DRectifying(latx, laty) / DIsometric(latx, laty);
            }
            else
            {
                psix *= Degree;
                psiy *= Degree;
                return DConformalToRectifying(Gd(psix), Gd(psiy)) * Dgd(psix, psiy);
            }
        }

        // (psix - psiy) / (mux - muy)
        internal double DRectifyingToIsometric(double mux, double muy)
        {
            double
              latx = _ell.InverseRectifyingLatitude(mux / Degree),
              laty = _ell.InverseRectifyingLatitude(muy / Degree);
            return _exact ?
              DIsometric(latx, laty) / DRectifying(latx, laty) :
              Dgdinv(Taupf(Tand(latx), _ell._es),
                     Taupf(Tand(laty), _ell._es)) *
              DRectifyingToConformal(mux, muy);
        }

        internal double MeanSinXi(double psix, double psiy)
            => Dlog(Cosh(psix), Cosh(psiy)) * Dcosh(psix, psiy)
                   + SinCosSeries(false, Gd(psix), Gd(psiy), _R.Span, maxpow_) * Dgd(psix, psiy);

        #endregion

        #region IGeodesicLike Members

        /// <summary>
        /// This function is not supposed to be called from user code directly.
        /// Use <see cref="GenDirect(double, double, double, double, GeodesicFlags, out double, out double, out double)"/> or
        /// <see cref="Direct(double, double, double, double, out double, out double)"/> instead.
        /// </summary>
        /// <param name="lat1"></param>
        /// <param name="lon1"></param>
        /// <param name="azi12"></param>
        /// <param name="_1"></param>
        /// <param name="s12"></param>
        /// <param name="outmask"></param>
        /// <param name="lat2"></param>
        /// <param name="lon2"></param>
        /// <param name="_2"></param>
        /// <param name="_3"></param>
        /// <param name="_4"></param>
        /// <param name="_5"></param>
        /// <param name="_6"></param>
        /// <param name="S12"></param>
        /// <returns></returns>
        public double GenDirect(double lat1, double lon1, double azi12, bool _1, double s12, GeodesicFlags outmask,
            out double lat2, out double lon2, out double _2, out double _3, out double _4, out double _5, out double _6, out double S12)
        {
            _2 = _3 = _4 = _5 = _6 = 0;
            GenDirect(lat1, lon1, azi12, s12, outmask, out lat2, out lon2, out S12);
            return 0;
        }

        /// <summary>
        /// This function is not supposed to be called from user code directly.
        /// Use <see cref="GenInverse(double, double, double, double, GeodesicFlags, out double, out double, out double)"/> or
        /// <see cref="Inverse(double, double, double, double, out double, out double)"/> instead.
        /// </summary>
        /// <param name="lat1"></param>
        /// <param name="lon1"></param>
        /// <param name="lat2"></param>
        /// <param name="lon2"></param>
        /// <param name="outmask"></param>
        /// <param name="s12"></param>
        /// <param name="azi12"></param>
        /// <param name="_1"></param>
        /// <param name="_2"></param>
        /// <param name="_3"></param>
        /// <param name="_4"></param>
        /// <param name="S12"></param>
        /// <returns></returns>
        public double GenInverse(double lat1, double lon1, double lat2, double lon2, GeodesicFlags outmask,
            out double s12, out double azi12, out double _1, out double _2, out double _3, out double _4, out double S12)
        {
            _1 = _2 = _3 = _4 = 0;
            GenInverse(lat1, lon1, lat2, lon2, outmask, out s12, out azi12, out S12);
            return 0;
        }

        #endregion

        /// <summary>
        /// The general direct rhumb problem.
        /// <see cref="Direct(double, double, double, double, out double, out double)"/> is defined in terms of this function.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi12">azimuth of the rhumb line (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="S12">area under the rhumb line (meters^2).</param>
        /// <param name="outmask">
        /// a bitor'ed combination of <see cref="GeodesicFlags"/> values specifying which of the following parameters should be set.
        /// </param>
        /// <remarks>
        /// The <see cref="GeodesicFlags"/> values possible for <paramref name="outmask"/> are
        /// <list type="bullet">
        /// <item>outmask |= <see cref="GeodesicFlags.Latitude"/> for the latitude <paramref name="lat2"/>;</item>
        /// <item>outmask |= <see cref="GeodesicFlags.Longitude"/> for the longitude <paramref name="lon2"/>;</item>
        /// <item>outmask |= <see cref="GeodesicFlags.Area"/> for the area <paramref name="S12"/>;</item>
        /// <item>outmask |= <see cref="GeodesicFlags.All"/> for all of the above;</item>
        /// <item>outmask |= <see cref="GeodesicFlags.LongUnroll"/> to unroll <paramref name="lon2"/> instead of wrapping it into the range [−180°, 180°].</item>
        /// </list>
        /// With the <see cref="GeodesicFlags.LongUnroll"/> bit set, the quantity <paramref name="lon2"/> − <paramref name="lon1"/> indicates
        /// how many times and in what sense the rhumb line encircles the ellipsoid.
        /// </remarks>
        public void GenDirect(double lat1, double lon1, double azi12, double s12,
            GeodesicFlags outmask, out double lat2, out double lon2, out double S12)
            => Line(lat1, lon1, azi12).GenPosition(s12, outmask, out lat2, out lon2, out S12);

        /// <summary>
        /// Solve the direct rhumb problem.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi12">azimuth of the rhumb line (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
        /// <param name="outmask">
        /// a bitor'ed combination of <see cref="GeodesicFlags"/> values specifying
        /// which of the properties in returned <see cref="DirectRhumbResult"/> instance should be set.
        /// </param>
        /// <returns>
        /// A <see cref="DirectRhumbResult"/> instance containing the result of the calcutation.
        /// </returns>
        /// <remarks>
        /// The <see cref="GeodesicFlags"/> values possible for <paramref name="outmask"/> are
        /// <list type="bullet">
        /// <item><i>outmask</i> |= <see cref="GeodesicFlags.Latitude"/> for the latitude returned in <see cref="DirectRhumbResult.Latitude"/>;</item>
        /// <item><i>outmask</i> |= <see cref="GeodesicFlags.Longitude"/> for the longitude returned in <see cref="DirectRhumbResult.Longitude"/>;</item>
        /// <item><i>outmask</i> |= <see cref="GeodesicFlags.Area"/> for the area returned in <see cref="RhumbResult.Area"/>;</item>
        /// <item><i>outmask</i> |= <see cref="GeodesicFlags.All"/> for all of the above;</item>
        /// <item><i>outmask</i> |= <see cref="GeodesicFlags.LongUnroll"/> to unroll <see cref="DirectRhumbResult.Longitude"/> instead of wrapping it into the range [−180°, 180°].</item>
        /// </list>
        /// With the <see cref="GeodesicFlags.LongUnroll"/> bit set, the quantity <see cref="DirectRhumbResult.Longitude"/> − <paramref name="lon1"/> indicates
        /// how many times and in what sense the rhumb line encircles the ellipsoid.
        /// </remarks>
        public DirectRhumbResult Direct(double lat1, double lon1, double azi12, double s12, GeodesicFlags outmask = GeodesicFlags.All)
            => Line(lat1, lon1, azi12).Position(s12, outmask);

        /// <summary>
        /// Solve the direct rhumb problem returning also the area.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi12">azimuth of the rhumb line (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="S12">area under the rhumb line (meters^2).</param>
        /// <remarks>
        /// <paramref name="lat1"/> should be in the range [−90°, 90°]. The value of <paramref name="lon2"/> returned is in the range [−180°, 180°].
        /// <para>
        /// If point 1 is a pole, the cosine of its latitude is taken to be 1/ε2 (where ε is 2^-52).
        /// This position, which is extremely close to the actual pole, allows the calculation to be carried out in finite terms.
        /// If <paramref name="s12"/> is large enough that the rhumb line crosses a pole, the longitude of point 2 is indeterminate
        /// (a <see cref="double.NaN"/> is returned for <paramref name="lon2"/> and <paramref name="S12"/>).
        /// </para>
        /// </remarks>
        public void Direct(double lat1, double lon1, double azi12, double s12,
            out double lat2, out double lon2, out double S12)
            => GenDirect(lat1, lon1, azi12, s12,
                      GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Area, out lat2, out lon2, out S12);

        /// <summary>
        /// Solve the direct rhumb problem without the area.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi12">azimuth of the rhumb line (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        public void Direct(double lat1, double lon1, double azi12, double s12,
                out double lat2, out double lon2)
            => GenDirect(lat1, lon1, azi12, s12, GeodesicFlags.Latitude | GeodesicFlags.Longitude, out lat2, out lon2, out _);

        /// <summary>
        /// The general inverse rhumb problem.
        /// <see cref="Inverse(double, double, double, double, out double, out double)"/> is defined in terms of this function.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi12">azimuth of the rhumb line (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters).</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="S12">area under the rhumb line (meters^2).</param>
        /// <param name="outmask">
        /// a bitor'ed combination of <see cref="GeodesicFlags"/> values specifying which of the following parameters should be set.
        /// </param>
        /// <remarks>
        /// The <see cref="GeodesicFlags"/> values possible for <paramref name="outmask"/> are
        /// <list type="bullet">
        /// <item>outmask |= <see cref="GeodesicFlags.Distance"/> for the distance <paramref name="s12"/>;</item>
        /// <item>outmask |= <see cref="GeodesicFlags.Azimuth"/> for the rhumb line azimuth <paramref name="azi12"/>;</item>
        /// <item>outmask |= <see cref="GeodesicFlags.Area"/> for the area <paramref name="S12"/>;</item>
        /// <item>outmask |= <see cref="GeodesicFlags.All"/> for all of the above;</item>
        /// </list>
        /// </remarks>
        public void GenInverse(double lat1, double lon1, double lat2, double lon2,
                          GeodesicFlags outmask,
                          out double s12, out double azi12, out double S12)
        {
            double
              lon12 = AngDiff(lon1, lon2),
              psi1 = _ell.IsometricLatitude(lat1),
              psi2 = _ell.IsometricLatitude(lat2),
              psi12 = psi2 - psi1,
              h = Hypot(lon12, psi12);

            s12 = azi12 = S12 = double.NaN;

            if (outmask.HasFlag(GeodesicFlags.Azimuth))
                azi12 = Atan2d(lon12, psi12);
            if (outmask.HasFlag(GeodesicFlags.Distance))
            {
                var dmudpsi = DIsometricToRectifying(psi2, psi1);
                s12 = h * dmudpsi * _ell.QuarterMeridian / QD;
            }
            if (outmask.HasFlag(GeodesicFlags.Area))
                S12 = _c2 * lon12 *
                  MeanSinXi(psi2 * Degree, psi1 * Degree);
        }

        /// <summary>
        /// Solve the inverse rhumb problem.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="outmask">
        /// a bitor'ed combination of <see cref="GeodesicFlags"/> values specifying
        /// which of the properties in returned <see cref="InverseRhumbResult"/> instance should be set.
        /// </param>
        /// <returns>A <see cref="InverseRhumbResult"/> instance containing the result of the calcutation.</returns>
        /// <remarks>
        /// The <see cref="GeodesicFlags"/> values possible for <paramref name="outmask"/> are
        /// <list type="bullet">
        /// <item>outmask |= <see cref="GeodesicFlags.Distance"/> for the distance, <see cref="InverseRhumbResult.Distance"/>;</item>
        /// <item>outmask |= <see cref="GeodesicFlags.Azimuth"/> for the rhumb line azimuth, <see cref="InverseRhumbResult.Azimuth"/>;</item>
        /// <item>outmask |= <see cref="GeodesicFlags.Area"/> for the area, <see cref="RhumbResult.Area"/>;</item>
        /// <item>outmask |= <see cref="GeodesicFlags.All"/> for all of the above;</item>
        /// </list>
        /// </remarks>
        public InverseRhumbResult Inverse(double lat1, double lon1, double lat2, double lon2, GeodesicFlags outmask = GeodesicFlags.All)
        {
            GenInverse(lat1, lon1, lat2, lon2, outmask, out var s12, out var azi12, out var S12);
            return new InverseRhumbResult
            {
                Area = S12,
                Azimuth = azi12,
                Distance = s12
            };
        }

        /// <summary>
        /// Solve the inverse rhumb problem returning also the area.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi12">azimuth of the rhumb line (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters).</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="S12">area under the rhumb line (meters^2).</param>
        /// <remarks>
        /// The shortest rhumb line is found. If the end points are on opposite meridians, there are two shortest rhumb lines and the east-going one
        /// is chosen. <paramref name="lat1"/> and <paramref name="lat2"/> should be in the range [−90°, 90°].
        /// The value of <paramref name="azi12"/> returned is in the range [−180°, 180°].
        /// <para>
        /// If either point is a pole, the cosine of its latitude is taken to be 1/ε2 (where ε is 2^-52).
        /// This position, which is extremely close to the actual pole, allows the calculation to be carried out in finite terms.
        /// </para>
        /// </remarks>
        public void Inverse(double lat1, double lon1, double lat2, double lon2,
             out double s12, out double azi12, out double S12)
            => GenInverse(lat1, lon1, lat2, lon2,
                       GeodesicFlags.Distance | GeodesicFlags.Azimuth | GeodesicFlags.Area, out s12, out azi12, out S12);

        /// <summary>
        /// Solve the inverse rhumb problem without the area.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi12">azimuth of the rhumb line (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters).</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        public void Inverse(double lat1, double lon1, double lat2, double lon2,
                 out double s12, out double azi12)
            => GenInverse(lat1, lon1, lat2, lon2, GeodesicFlags.Distance | GeodesicFlags.Azimuth, out s12, out azi12, out _);

        /// <summary>
        /// Set up to compute several points on a single rhumb line.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi12">azimuth of the rhumb line (degrees).</param>
        /// <returns>a <see cref="RhumbLine"/> object.</returns>
        /// <remarks>
        /// <paramref name="lat1"/> should be in the range [−90°, 90°].
        /// If point 1 is a pole, the cosine of its latitude is taken to be 1/ε2 (where ε is 2^-52).
        /// This position, which is extremely close to the actual pole, allows the calculation to be carried out in finite terms.
        /// </remarks>
        public RhumbLine Line(double lat1, double lon1, double azi12)
            => new RhumbLine(this, lat1, lon1, azi12);

    }
}
