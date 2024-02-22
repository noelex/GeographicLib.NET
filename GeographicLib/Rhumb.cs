using System;
using static GeographicLib.Macros;
using static GeographicLib.MathEx;
using static System.Math;

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
    /// This implementation is described in
    /// <para>
    /// - C. F. F. Karney,
    ///   <a href="https://arxiv.org/abs/2303.03219" > The area of rhumb
    /// polygons</a>,
    /// Technical Report, SRI International, March 2023.
    ///   <a href="https://arxiv.org/abs/2303.03219" > arxiv:2303.03219</a>
    /// .
    /// </para>
    /// <para>
    /// For more information on rhumb lines see <a href="https://geographiclib.sourceforge.io/html/rhumb.html">Rhumb lines</a>.
    /// </para>
    /// </remarks>
    public class Rhumb : IGeodesicLike
    {
        internal readonly DAuxLatitude _aux;
        private readonly bool _exact;
        internal readonly double _a, _f, _n, _rm, _c2;

        private readonly int _lL; // N.B. names of the form _[A-Z].* are reserved in C++
        private readonly double[] _pP; // The Fourier coefficients P_l

        private const int Lmax_ = GEOGRAPHICLIB_RHUMBAREA_ORDER;

        private static readonly double[] s_coeffs = new[] {
            // Coefficients in matrix Q
            138734126/(double) (638512875), -102614/(double) (467775), 596/(double) (2025),
            -398/(double) (945), 22/(double) (45), -1/(double) (3),
            17749373/(double) (425675250), -24562/(double) (155925), 1543/(double) (4725),
            -118/(double) (315), 1/(double) (5),
            1882432/(double) (8513505), -38068/(double) (155925), 152/(double) (945),
            -17/(double) (315),
            268864/(double) (2027025), -752/(double) (10395), 5/(double) (252),
            62464/(double) (2027025), -101/(double) (17325),
            11537/(double) (4054050),
        };

        /// <inheritdoc/>
        public double EquatorialRadius => _a;

        /// <inheritdoc/>
        public double Flattening => _f;

        /// <inheritdoc/>
        public double EllipsoidArea
            // _c2 contains a Math::degrees() factor, so 4*pi -> 2*Math::td.
            => 2 * TD * _c2;

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
        /// if <see langword="true"/> use the exact expressions for the auxiliary latitudes; otherwise use series expansion(accurate for |<i>f</i>| &lt; 0.01)
        /// </param>
        public Rhumb(double a, double f, bool exact = false)
        {
            _aux = new DAuxLatitude(a, f);
            _exact = exact;
            _a = a;
            _f = f;
            _n = _f / (2 - _f);
            _rm = _aux.RectifyingRadius(_exact);
            _c2 = _aux.AuthalicRadiusSquared(_exact) * Degree;
            _lL = _exact ? 8 : Lmax_;   // 8 is starting size for DFT fit
            _pP = new double[_lL];

            AreaCoeffs(ref _lL, ref _pP);
        }

        /// <summary>
        /// Initialize a new <see cref="Rhumb"/> instance with specified ellipsoid.
        /// </summary>
        /// <param name="ellipsoid">the ellipoid to be used for calculation.</param>
        /// <param name="exact">
        /// if <see langword="true"/> use the exact expressions for the auxiliary latitudes; otherwise use series expansion(accurate for |<i>f</i>| &lt; 0.01)
        /// </param>
        public Rhumb(IEllipsoid ellipsoid, bool exact = false) : this(ellipsoid.EquatorialRadius, ellipsoid.Flattening, exact) { }

        /// <summary>
        /// Initialize a new <see cref="Rhumb"/> instance with specified ellipsoid.
        /// </summary>
        /// <param name="ellipsoid">the ellipoid to be used for calculation.</param>
        /// <param name="exact">
        /// if <see langword="true"/> use the exact expressions for the auxiliary latitudes; otherwise use series expansion(accurate for |<i>f</i>| &lt; 0.01)
        /// </param>
        public Rhumb(Ellipsoid ellipsoid, bool exact = false)
            : this(ellipsoid.EquatorialRadius, ellipsoid.Flattening, exact)
        {
        }

        #region Private Methods

        private double qIntegrand(double beta)
        {
            // pbeta(beta) = integrate(q(beta), beta)
            //   q(beta) = (1-f) * (sin(xi) - sin(chi)) / cos(phi)
            //           = (1-f) * (cos(chi) - cos(xi)) / cos(phi) *
            //   (cos(xi) + cos(chi)) / (sin(xi) + sin(chi))
            // Fit q(beta)/cos(beta) with Fourier transform
            //   q(beta)/cos(beta) = sum(c[k] * sin((2*k+1)*beta), k, 0, K-1)
            // then the integral is
            //   pbeta = sum(d[k] * cos((2*k+2)*beta), k, 0, K-1)
            // where
            //   d[k] = -1/(4*(k+1)) * (c[k] + c[k+1]) for k in 0..K-2
            //   d[K-1] = -1/(4*K) * c[K-1]
            AuxAngle betaa = AuxAngle.FromRadians(beta),
                     phia = _aux.Convert(AuxLatitudeType.Beta, AuxLatitudeType.Phi,
                                        betaa, true).Normalized(),
                     chia = _aux.Convert(AuxLatitudeType.Phi, AuxLatitudeType.Chi,
                                        phia, true).Normalized(),
                     xia = _aux.Convert(AuxLatitudeType.Phi, AuxLatitudeType.Xi,
                                        phia, true).Normalized();
            double schi = chia.Y, cchi = chia.X, sxi = xia.Y, cxi = xia.X,
              cphi = phia.X, cbeta = betaa.X;
            return (1 - _aux.Flattening) *
              (Abs(schi) < Abs(cchi) ? sxi - schi :
                (cchi - cxi) * (cxi + cchi) / (sxi + schi)) / (cphi * cbeta);
        }

        private void AreaCoeffs(ref int _lL, ref double[] _pP)
        {
            // Set up coefficients for area calculation
            if (_exact)
            {
                // Compute coefficients by Fourier transform of integrand
                const double eps = DBL_EPSILON / 2;
                Func<double, double> f = qIntegrand;
                int L = 4;
                var c = new double[L];
                var fft = new DST(L); fft.Transform(f, c.AsSpan()); L *= 2;
                // For |n| <= 0.99, actual max for doubles is 2163.  This scales as
                // Math::digits() and for long doubles (GEOGRAPHICLIB_PRECISION = 3,
                // digits = 64), this becomes 2163 * 64 / 53 = 2612.  Round this up to
                // 2^12 = 4096 and scale this by Math::digits()/64 if digits() > 64.
                //
                // 64 = digits for long double, 6 = 12 - log2(64)
                int Lmax = 1 << ((int)(Ceiling(Log2(Max(DBL_MANT_DIG, 64)))) + 6);
                for (_lL = 0; L <= Lmax && _lL == 0; L *= 2)
                {

                    fft = new DST(L / 2); Array.Resize(ref c, L); fft.Refine(f, c.AsSpan());
                    Array.Resize(ref _pP, L);
                    for (int l = 0, k = -1; l < L; ++l)
                    {
                        // Compute Fourier coefficients of integral
                        _pP[l] = (c[l] + (l + 1 < L ? c[l + 1] : 0)) / (-4 * (l + 1));
                        if (Abs(_pP[l]) <= eps)
                        {
                            if (k < 0) k = l;   // mark as first small value
                        }
                        else
                            k = -1;             // run interrupted
                        if (k >= 0 && l - k + 1 >= (l + 1 + 7) / 8)
                        {
                            // run of small values of at least l/8?
                            _lL = l + 1; Array.Resize(ref _pP, _lL); break;
                        }
                    }
                    // loop exits if _lL > 0
                }
                if (_lL == 0)          // Hasn't converged -- just use the values we have
                    _lL = _pP.Length;
            }
            else
            {
                // Use series expansions in n for Fourier coeffients of the integral
                // See "Series expansions for computing rhumb areas"
                // https://doi.org/10.5281/zenodo.7685484
                double d = 1;
                int o = 0;
                for (int l = 0; l < Lmax_; ++l)
                {
                    int m = Lmax_ - l - 1;
                    d *= _n;
                    _pP[l] = d * PolyVal(m, s_coeffs.AsSpan().Slice(o), _n);
                    o += m + 1;
                }
            }
            // Post condition: o == sizeof(coeffs) / sizeof(real)
        }

        internal double MeanSinXi(AuxAngle chix, AuxAngle chiy)
        {
            AuxAngle
              phix = _aux.Convert(AuxLatitudeType.Chi, AuxLatitudeType.Phi, chix, _exact),
              phiy = _aux.Convert(AuxLatitudeType.Chi, AuxLatitudeType.Phi, chiy, _exact),
              betax = _aux.Convert(AuxLatitudeType.Phi, AuxLatitudeType.Beta, phix, _exact).Normalized(),
              betay = _aux.Convert(AuxLatitudeType.Phi, AuxLatitudeType.Beta, phiy, _exact).Normalized();
            double DpbetaDbeta =
              DAuxLatitude.DClenshaw(false,
                                betay.Radians - betax.Radians,
                                betax.Y, betax.X, betay.Y, betay.X,
                                _pP.AsSpan(), _lL),
              tx = chix.Tan, ty = chiy.Tan,
              DbetaDpsi = _exact ?
              _aux.DParametric(phix, phiy) / _aux.DIsometric(phix, phiy) :
              _aux.DConvert(AuxLatitudeType.Chi, AuxLatitudeType.Beta, chix, chiy) /
              DAuxLatitude.Dlam(tx, ty);
            return DAuxLatitude.Dp0Dpsi(tx, ty) + DpbetaDbeta * DbetaDpsi;
        }

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
        double IGeodesicLike.GenDirect(double lat1, double lon1, double azi12, bool _1, double s12, GeodesicFlags outmask,
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
        double IGeodesicLike.GenInverse(double lat1, double lon1, double lat2, double lon2, GeodesicFlags outmask,
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
            s12 = azi12 = S12 = double.NaN;

            AuxAngle phi1 = AuxAngle.FromDegrees(lat1), phi2 = AuxAngle.FromDegrees(lat2),
                     chi1 = _aux.Convert(AuxLatitudeType.Phi, AuxLatitudeType.Chi, phi1, _exact),
                     chi2 = _aux.Convert(AuxLatitudeType.Phi, AuxLatitudeType.Chi, phi2, _exact);
            double
              lon12 = AngDiff(lon1, lon2),
              lam12 = lon12 * Degree,
              psi1 = chi1.Lam,
              psi2 = chi2.Lam,
              psi12 = psi2 - psi1;
            if (outmask.HasFlag(GeodesicFlags.Azimuth))
                azi12 = Atan2d(lam12, psi12);
            if (outmask.HasFlag(GeodesicFlags.Distance))
            {
                if (double.IsInfinity(psi1) || double.IsInfinity(psi2))
                {
                    s12 = Abs(_aux.Convert(AuxLatitudeType.Phi, AuxLatitudeType.Mu,
                                            phi2, _exact).Radians -
                               _aux.Convert(AuxLatitudeType.Phi, AuxLatitudeType.Mu,
                                            phi1, _exact).Radians) * _rm;
                }
                else
                {
                    double h = Hypot(lam12, psi12);
                    // dmu/dpsi = dmu/dchi / dpsi/dchi
                    double dmudpsi = _exact ?
                      _aux.DRectifying(phi1, phi2) / _aux.DIsometric(phi1, phi2) :
                      _aux.DConvert(AuxLatitudeType.Chi, AuxLatitudeType.Mu, chi1, chi2)
                      / DAuxLatitude.Dlam(chi1.Tan, chi2.Tan);
                    s12 = h * dmudpsi * _rm;
                }
            }
            if (outmask.HasFlag(GeodesicFlags.Area))
                S12 = _c2 * lon12 * MeanSinXi(chi1, chi2);
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
