using System;
using System.Collections.Generic;
using System.Text;

using static System.Math;
using static GeographicLib.Macros;
using static GeographicLib.MathEx;

namespace GeographicLib
{
    /// <summary>
    /// Provide exact geodesic calculations.
    /// </summary>
    public partial class GeodesicExact : GeodesicBase
    {
        internal const int
            nC4_ = GEOGRAPHICLIB_GEODESICEXACT_ORDER,
            nC4x_ = (nC4_ * (nC4_ + 1)) / 2;
        private const uint
            maxit1_ = 20,
            maxit2_ = maxit1_ + DBL_MANT_DIG + 10;

        private const double tol0_ = DBL_EPSILON;
        internal static readonly double
            tiny_ = Sqrt(DBL_MIN),
            tol1_ = 200 * tol0_,
            tol2_ = Sqrt(tol0_),
            tolb_ = tol0_ * tol2_, // Check on bisection interval
            xthresh_ = 1000 * tol2_;

        internal readonly double _a, _f, _f1, _e2, _ep2, _n, _b, _c2, _etol2;
        private readonly Memory<double> _C4x = new double[nC4x_];

        /// <summary>
        /// Initialize a new <see cref="GeodesicExact"/> instance with specified ellipoid.
        /// </summary>
        /// <param name="ellipsoid">Source <see cref="IEllipsoid"/> object.</param>
        public GeodesicExact(IEllipsoid ellipsoid) : this(ellipsoid.EquatorialRadius, ellipsoid.Flattening) { }

        /// <summary>
        /// Initialize a new <see cref="GeodesicExact"/> instance with specified equatorial radius and flattening of the ellipsoid.
        /// </summary>
        /// <param name="a">equatorial radius (meters).</param>
        /// <param name="f">flattening of ellipsoid.  Setting <i>f</i> = 0 gives a sphere.</param>
        public GeodesicExact(double a, double f)
        {
            _a = a;
            _f = f;
            _f1 = 1 - f;
            _e2 = _f * (2 - _f);
            _ep2 = _e2 / Sq(_f1); // e2 / (1 - e2)
            _n = _f / (2 - _f);
            _b = _a * _f1;

            // The Geodesic class substitutes atanh(sqrt(e2)) for asinh(sqrt(ep2)) in
            // the definition of _c2.  The latter is more accurate for very oblate
            // ellipsoids (which the Geodesic class does not attempt to handle).  Of
            // course, the area calculation in GeodesicExact is still based on a
            // series and so only holds for moderately oblate (or prolate)
            // ellipsoids.
            _c2 = (Sq(_a) + Sq(_b) *
                  (_f == 0 ? 1 :
                   (_f > 0 ? Asinh(Sqrt(_ep2)) : Atan(Sqrt(-_e2))) /
                   Sqrt(Abs(_e2)))) / 2; // authalic radius squared

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
                    Sqrt(Max(0.001, Abs(_f)) * Min(1d, 1 - _f / 2) / 2);

            if (!(IsFinite(_a) && _a > 0))
                throw new GeographicException("Equatorial radius is not positive");

            if (!(IsFinite(_b) && _b > 0))
                throw new GeographicException("Polar semi-axis is not positive");

            C4coeff();
        }

        /// <summary>
        /// A global instantiation of <see cref="GeodesicExact"/> with the parameters for the WGS84 ellipsoid.
        /// </summary>
        public static GeodesicExact WGS84 { get; } = new GeodesicExact(Ellipsoid.WGS84);

        /// <inheritdoc/>
        public override double EquatorialRadius => _a;

        /// <inheritdoc/>
        public override double Flattening => _f;

        /// <inheritdoc/>
        public override double EllipsoidArea => 4 * PI * _c2;

        #region Private methods

        internal static double CosSeries(double sinx, double cosx, ReadOnlySpan<double> c, int n)
        {
            // Evaluate
            // y = sum(c[i] * cos((2*i+1) * x), i, 0, n-1)
            // using Clenshaw summation.
            // Approx operation count = (n + 5) mult and (2 * n + 2) add
            var i = n;                    // Point to one beyond last element
            double
              ar = 2 * (cosx - sinx) * (cosx + sinx), // 2 * cos(2 * x)
              y0 = (n & 1) != 0 ? c[--i] : 0, y1 = 0;          // accumulators for sum
                                                               // Now n is even
            n /= 2;
            while (n-- > 0)
            {
                // Unroll loop x 2, so accumulators return to their original role
                y1 = ar * y0 - y1 + c[--i];
                y0 = ar * y1 - y0 + c[--i];
            }
            return cosx * (y0 - y1);    // cos(x) * (y0 - y1)
        }

        private static double Astroid(double x, double y)
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

        private void Lengths(EllipticFunction E,
                 double sig12,
                 double ssig1, double csig1, double dn1,
                 double ssig2, double csig2, double dn2,
                 double cbet1, double cbet2, GeodesicFlags outmask,
                 out double s12b, out double m12b, out double m0,
                 out double M12, out double M21)
        {
            // Return m12b = (reduced length)/_b; also calculate s12b = distance/_b,
            // and m0 = coefficient of secular term in expression for reduced length.
            s12b = m12b = m0 = M12 = M21 = double.NaN;

            outmask = outmask.Flags();
            // outmask & DISTANCE: set s12b
            // outmask & REDUCEDLENGTH: set m12b & m0
            // outmask & GEODESICSCALE: set M12 & M21

            // It's OK to have repeated dummy arguments,
            // e.g., s12b = m0 = M12 = M21 = dummy

            if (outmask.HasAny(GeodesicFlags.Distance))
                // Missing a factor of _b
                s12b = E.E() / (PI / 2) *
                  (sig12 + (E.DeltaE(ssig2, csig2, dn2) - E.DeltaE(ssig1, csig1, dn1)));
            if (outmask.HasAny(GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale))
            {
                double
                  m0x = -E.K2 * E.D() / (PI / 2),
                  J12 = m0x *
                  (sig12 + (E.DeltaD(ssig2, csig2, dn2) - E.DeltaD(ssig1, csig1, dn1)));
                if (outmask.HasAny(GeodesicFlags.ReducedLength))
                {
                    m0 = m0x;
                    // Missing a factor of _b.  Add parens around (csig1 * ssig2) and
                    // (ssig1 * csig2) to ensure accurate cancellation in the case of
                    // coincident points.
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
        }

        private double InverseStart(EllipticFunction E,
                  double sbet1, double cbet1, double dn1,
                  double sbet2, double cbet2, double dn2,
                  double lam12, double slam12, double clam12,
                  out double salp1, out double calp1,
                  out double salp2, out double calp2, out double dnm)
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

                var lam12x = Atan2(-slam12, -clam12); // lam12 - pi
                if (_f >= 0)
                {            // In fact f == 0 does not get here
                             // x = dlong, y = dlat
                    {
                        var k2 = Sq(sbet1) * _ep2;
                        E.Reset(-k2, -_ep2, 1 + k2, 1 + _ep2);
                        lamscale = _e2 / _f1 * cbet1 * 2 * E.H();
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
                    Lengths(E, PI + bet12a,
                            sbet1, -cbet1, dn1, sbet2, cbet2, dn2,
                            cbet1, cbet2, GeodesicFlags.ReducedLength, out _, out var m12b, out var m0, out _, out _);
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
                        calp1 = Max(x > -tol1_ ? 0d : -1d, x);
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
                      EllipticFunction E,
                      out double domg12, bool diffp, out double dlam12)
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

            double somg1, comg1, somg2, comg2, somg12, comg12, cchi1, cchi2, lam12;
            // tan(bet1) = tan(sig1) * cos(alp1)
            // tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1)
            ssig1 = sbet1; somg1 = salp0 * sbet1;
            csig1 = comg1 = calp1 * cbet1;
            // Without normalization we have schi1 = somg1.
            cchi1 = _f1 * dn1 * comg1;
            Norm(ref ssig1, ref csig1);
            // Math::norm(somg1, comg1); -- don't need to normalize!
            // Math::norm(schi1, cchi1); -- don't need to normalize!

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
            // Without normalization we have schi2 = somg2.
            cchi2 = _f1 * dn2 * comg2;
            Norm(ref ssig2, ref csig2);
            // Math::norm(somg2, comg2); -- don't need to normalize!
            // Math::norm(schi2, cchi2); -- don't need to normalize!

            // sig12 = sig2 - sig1, limit to [0, pi]
            sig12 = Atan2(Max(0, csig1 * ssig2 - ssig1 * csig2),
                                       csig1 * csig2 + ssig1 * ssig2);

            // omg12 = omg2 - omg1, limit to [0, pi]
            somg12 = Max(0, comg1 * somg2 - somg1 * comg2);
            comg12 = comg1 * comg2 + somg1 * somg2;
            var k2 = Sq(calp0) * _ep2;
            E.Reset(-k2, -_ep2, 1 + k2, 1 + _ep2);
            // chi12 = chi2 - chi1, limit to [0, pi]
            double
              schi12 = Max(0, cchi1 * somg2 - somg1 * cchi2),
              cchi12 = cchi1 * cchi2 + somg1 * somg2;
            // eta = chi12 - lam120
            var eta = Atan2(schi12 * clam120 - cchi12 * slam120,
                             cchi12 * clam120 + schi12 * slam120);
            var deta12 = -_e2 / _f1 * salp0 * E.H() / (PI / 2) *
              (sig12 + (E.DeltaH(ssig2, csig2, dn2) - E.DeltaH(ssig1, csig1, dn1)));
            lam12 = eta + deta12;
            // domg12 = deta12 + chi12 - omg12
            domg12 = deta12 + Atan2(schi12 * comg12 - cchi12 * somg12,
                                    cchi12 * comg12 + schi12 * somg12);
            if (diffp)
            {
                if (calp2 == 0)
                    dlam12 = -2 * _f1 * dn1 / sbet1;
                else
                {
                    Lengths(E, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                            cbet1, cbet2, GeodesicFlags.ReducedLength,
                            out _, out dlam12, out _, out _, out _);
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
            m12 = M12 = M21 = s12 = S12 = double.NaN;
            salp1 = calp1 = salp2 = calp2 = double.NaN;

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
            // Swap points so that point with higher (abs) latitude is point 1
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

            double s12x = 0, m12x = 0;
            // Initialize for the meridian.  No longitude calculation is done in this
            // case to let the parameter default to 0.
            var E = new EllipticFunction(-_ep2);

            SinCosd(lat1, out var sbet1, out var cbet1); sbet1 *= _f1;
            // Ensure cbet1 = +epsilon at poles; doing the fix on beta means that sig12
            // will be <= 2*tiny for two points at the same pole.
            Norm(ref sbet1, ref cbet1); cbet1 = Max(tiny_, cbet1);

            SinCosd(lat2, out var sbet2, out var cbet2); sbet2 *= _f1;
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
              dn1 = (_f >= 0 ? Sqrt(1 + _ep2 * Sq(sbet1)) :
                     Sqrt(1 - _e2 * Sq(cbet1)) / _f1),
              dn2 = (_f >= 0 ? Sqrt(1 + _ep2 * Sq(sbet2)) :
                     Sqrt(1 - _e2 * Sq(cbet2)) / _f1);

            double a12=double.NaN, sig12;

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
                sig12 = Atan2(Max(0, csig1 * ssig2 - ssig1 * csig2),
                                           csig1 * csig2 + ssig1 * ssig2);
                {
                    Lengths(E, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                            cbet1, cbet2, outmask | GeodesicFlags.ReducedLength,
                            out s12x, out m12x, out _, out M12, out M21);
                }
                // Add the check for sig12 since zero length geodesics might yield m12 <
                // 0.  Test case was
                //
                //    echo 20.001 0 20.001 0 | GeodSolve -i
                //
                // In fact, we will have sig12 > pi/2 for meridional geodesic which is
                // not a shortest path.
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
                if (outmask.HasAny(GeodesicFlags.GeodesicScale))
                    M12 = M21 = Cos(sig12);
                a12 = lon12 / _f1;

            }
            else if (!meridian)
            {

                // Now point1 and point2 belong within a hemisphere bounded by a
                // meridian and geodesic is neither meridional or equatorial.

                // Figure a starting point for Newton's method
                sig12 = InverseStart(E, sbet1, cbet1, dn1, sbet2, cbet2, dn2,
                                     lam12, slam12, clam12,
                                     out salp1, out calp1, out salp2, out calp2, out var dnm);

                if (sig12 >= 0)
                {
                    // Short lines (InverseStart sets salp2, calp2, dnm)
                    s12x = sig12 * _b * dnm;
                    m12x = Sq(dnm) * _b * Sin(sig12 / dnm);
                    if (outmask.HasAny(GeodesicFlags.GeodesicScale))
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
                    double ssig1 = 0, csig1 = 0, ssig2 = 0, csig2 = 0, domg12 = 0;
                    uint numit = 0;
                    // Bracketing range
                    double salp1a = tiny_, calp1a = 1, salp1b = tiny_, calp1b = -1;
                    for (bool tripn = false, tripb = false;
                         numit < maxit2_ || GEOGRAPHICLIB_PANIC;
                         ++numit)
                    {
                        // 1/4 meridian = 10e6 m and random input.  max err is estimated max
                        // error in nm (checking solution of inverse problem by direct
                        // solution).  iter is mean and sd of number of iterations
                        //
                        //           max   iter
                        // log2(b/a) err mean  sd
                        //    -7     387 5.33 3.68
                        //    -6     345 5.19 3.43
                        //    -5     269 5.00 3.05
                        //    -4     210 4.76 2.44
                        //    -3     115 4.55 1.87
                        //    -2      69 4.35 1.38
                        //    -1      36 4.05 1.03
                        //     0      15 0.01 0.13
                        //     1      25 5.10 1.53
                        //     2      96 5.61 2.09
                        //     3     318 6.02 2.74
                        //     4     985 6.24 3.22
                        //     5    2352 6.32 3.44
                        //     6    6008 6.30 3.45
                        //     7   19024 6.19 3.30
                        double dv;
                        var v = Lambda12(sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1,
                                          slam12, clam12,
                                          out salp2, out calp2, out sig12, out ssig1, out csig1, out ssig2, out csig2,
                                          E, out domg12, numit < maxit1_, out dv);
                        // Reversed test to allow escape with NaNs
                        if (tripb || !(Abs(v) >= (tripn ? 8 : 1) * tol0_)) break;
                        // Update bracketing values
                        if (v > 0 && (numit > maxit1_ || calp1 / salp1 > calp1b / salp1b))
                        { salp1b = salp1; calp1b = calp1; }
                        else if (v < 0 && (numit > maxit1_ || calp1 / salp1 < calp1a / salp1a))
                        { salp1a = salp1; calp1a = calp1; }
                        if (numit < maxit1_ && dv > 0)
                        {
                            var
                              dalp1 = -v / dv;
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
                        Lengths(E, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                                cbet1, cbet2, outmask, out s12x, out m12x, out _, out M12, out M21);
                    }
                    m12x *= _b;
                    s12x *= _b;
                    a12 = sig12 / Degree;
                    if (outmask.HasAny(GeodesicFlags.GeodesicScale))
                    {
                        // omg12 = lam12 - domg12
                        double sdomg12 = Sin(domg12), cdomg12 = Cos(domg12);
                        somg12 = slam12 * cdomg12 - clam12 * sdomg12;
                        comg12 = clam12 * cdomg12 + slam12 * sdomg12;
                    }
                }
            }

            if (outmask.HasAny(GeodesicFlags.Distance))
                s12 = 0d + s12x;           // Convert -0 to 0

            if (outmask.HasAny(GeodesicFlags.ReducedLength))
                m12 = 0d + m12x;           // Convert -0 to 0

            if (outmask.HasAny(GeodesicFlags.Area))
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
                    Span<double> C4a = stackalloc double[nC4_];
                    C4f(eps, C4a);
                    double
                      B41 = CosSeries(ssig1, csig1, C4a, nC4_),
                      B42 = CosSeries(ssig2, csig2, C4a, nC4_);
                    S12 = A4 * (B42 - B41);
                }
                else
                    // Avoid problems with indeterminate sig1, sig2 on equator
                    S12 = 0;

                if (!meridian)
                {
                    if (somg12 > 1)
                    {
                        somg12 = Sin(omg12); comg12 = Cos(omg12);
                    }
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
                if (outmask.HasAny(GeodesicFlags.GeodesicScale))
                    Swap(ref M12, ref M21);
            }

            salp1 *= swapp * lonsign; calp1 *= swapp * latsign;
            salp2 *= swapp * lonsign; calp2 *= swapp * latsign;

            // Returned value in [0, 180]
            return a12;
        }

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
            if (!(o == nC4x_))
                throw new GeographicException("C4 misalignment");
        }

        #endregion

        /// <inheritdoc/>
        public override double GenInverse(double lat1, double lon1, double lat2, double lon2,
                       GeodesicFlags outmask,
                       out double s12, out double azi1, out double azi2,
                       out double m12, out double M12, out double M21, out double S12)
        {
            outmask = outmask.Flags();
            var a12 = GenInverse(lat1, lon1, lat2, lon2,
                               outmask, out s12, out var salp1, out var calp1, out var salp2, out var calp2,
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
        public override double GenDirect(double lat1, double lon1, double azi1,
                         bool arcmode, double s12_a12, GeodesicFlags outmask,
                         out double lat2, out double lon2, out double azi2,
                         out double s12, out double m12, out double M12, out double M21,
                         out double S12)
        {
            // Automatically supply DISTANCE_IN if necessary
            if (!arcmode) outmask |= GeodesicFlags.DistanceIn;

            return new GeodesicLineExact(this, lat1, lon1, azi1, outmask)
              .                         // Note the dot!
              GenPosition(arcmode, s12_a12, outmask,
                          out lat2, out lon2, out azi2, out s12, out m12, out M12, out M21, out S12);
        }

        /// <inheritdoc/>
        public override IGeodesicLine GenDirectLine(double lat1, double lon1, double azi1, bool arcmode, double s12_a12, GeodesicFlags caps = GeodesicFlags.All)
        {
            azi1 = AngNormalize(azi1);

            // Guard against underflow in salp0.  Also -0 is converted to +0.
            SinCosd(AngRound(azi1), out var salp1, out var calp1);

            // Automatically supply DISTANCE_IN if necessary
            if (!arcmode) caps |= GeodesicFlags.DistanceIn;

            return new GeodesicLineExact(this, lat1, lon1, azi1, salp1, calp1,
                                     caps, arcmode, s12_a12);
        }

        /// <inheritdoc/>
        public override IGeodesicLine Line(double lat1, double lon1, double azi1, GeodesicFlags caps = GeodesicFlags.All)
            => new GeodesicLineExact(this, lat1, lon1, azi1, caps);

        /// <inheritdoc/>
        public override IGeodesicLine InverseLine(double lat1, double lon1, double lat2, double lon2, GeodesicFlags caps = GeodesicFlags.All)
        {
            double  
                a12 = GenInverse(lat1, lon1, lat2, lon2,
                   // No need to specify AZIMUTH here
                   0u,out _, out var salp1, out var calp1, out var salp2, out var calp2,
                   out _, out _, out _, out _),
                azi1 = Atan2d(salp1, calp1);
            // Ensure that a12 can be converted to a distance
            if (caps.Flags().HasAny(GeodesicFlags.DistanceIn)) caps |= GeodesicFlags.Distance;
            return new GeodesicLineExact(this, lat1, lon1, azi1, salp1, calp1, caps, true, a12);
        }
    }
}
