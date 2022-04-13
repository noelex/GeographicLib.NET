using System;
using System.Collections.Generic;
using System.Text;

using static System.Math;
using static GeographicLib.MathEx;
using static GeographicLib.Macros;

namespace GeographicLib.Projections
{
    /// <summary>
    /// An exact implementation of the transverse Mercator projection.
    /// </summary>
    /// <remarks>
    /// Implementation of the Transverse Mercator Projection given in
    /// <list type="bullet">
    /// <item>
    /// L. P. Lee,
    /// <a href="https://doi.org/10.3138/X687-1574-4325-WM62">Conformal Projections Based On Jacobian Elliptic Functions</a>,
    /// Part V of Conformal Projections Based on Elliptic Functions,
    /// (B. V. Gutsell, Toronto, 1976), 128pp.,
    /// ISBN: 0919870163 (also appeared as: Monograph 16,
    /// Suppl. No. 1 to Canadian Cartographer, Vol 13).
    /// </item>
    /// <item>
    /// C. F. F. Karney,
    /// <a href="https://doi.org/10.1007/s00190-011-0445-3">Transverse Mercator with an accuracy of a few nanometers</a>,
    /// J. Geodesy 85(8), 475–485 (Aug. 2011);
    /// preprint <a href="https://arxiv.org/abs/1002.1417">arXiv:1002.1417</a>.
    /// </item>
    /// </list>
    /// Lee gives the correct results for forward and reverse transformations subject to the branch cut rules
    /// (see the description of the extendp argument to the constructor). The maximum error is about 8 nm (8 nanometers), ground distance,
    /// for the forward and reverse transformations. The error in the convergence is 2 × 10^-15", the relative error in the scale is 7 × 10^−12%%.
    /// See Sec. 3 of <a href="https://arxiv.org/abs/1002.1417">arXiv:1002.1417</a>for details.
    /// The method is "exact" in the sense that the errors are close to the round-off limit and that no changes are needed in the algorithms for
    /// them to be used with reals of a higher precision. Thus the errors using long double (with a 64-bit fraction) are about 2000 times smaller
    /// than using double (with a 53-bit fraction).
    /// <para>
    /// This algorithm is about 4.5 times slower than the 6th-order Krüger method, <see cref="TransverseMercator"/>, taking about 11 us for a
    /// combined forward and reverse projection on a 2.66 GHz Intel machine (g++, version 4.3.0, -O3).
    /// </para>
    /// <para>
    /// The ellipsoid parameters and the central scale are set in the constructor. The central meridian (which is a trivial shift of the longitude)
    /// is specified as the lon0 argument of the <see cref="Forward(double, double, double)"/> and <see cref="Reverse(double, double, double)"/> 
    /// functions. The latitude of origin is taken to be the equator. See the documentation on TransverseMercator for how to include a false easting,
    /// false northing, or a latitude of origin.
    /// </para>
    /// <para>
    /// See <a href="https://geographiclib.sourceforge.io/tm-grid.kmz">tm-grid.kmz</a>, for an illustration of the transverse Mercator grid in Google Earth.
    /// </para>
    /// <para>
    /// This class also returns the meridian convergence <i>gamma</i> and scale <i>k</i>.
    /// The meridian convergence is the bearing of grid north (the <i>y</i> axis) measured clockwise from true north.
    /// </para>
    /// <para>
    /// See <a href="https://geographiclib.sourceforge.io/html/transversemercator.html">Transverse Mercator projection</a> for a discussion of this projection.
    /// </para>
    /// </remarks>
    public class TransverseMercatorExact : IEllipsoid
    {
        private const int numit_ = 10;

        private static readonly double
            tol_ = DBL_EPSILON,
            tol2_ = 0.1 * tol_,
            taytol_ = Pow(tol_, 0.6);

        private readonly double _a, _f, _k0, _mu, _mv, _e;
        private readonly bool _extendp;
        private readonly EllipticFunction _Eu, _Ev;

        /// <summary>
        /// Initialize a new <see cref="TransverseMercatorExact"/> instance with specified equatorial radius, flattening and central scale factor.
        /// </summary>
        /// <param name="a">equatorial radius (meters).</param>
        /// <param name="f">flattening of ellipsoid.</param>
        /// <param name="k0">central scale factor.</param>
        /// <param name="extendp">use extended domain.</param>
        /// <remarks>
        /// <para>
        /// The transverse Mercator projection has a branch point singularity at <i>lat</i> = 0 and <i>lon</i> − <i>lon0</i> = 90 (1 − <i>e</i>)
        /// or (for <see cref="UTM"/>) <i>x</i> = 18381 km, <i>y</i> = 0m. The <paramref name="extendp"/> argument governs where the branch cut is placed.
        /// With <paramref name="extendp"/> = <see langword="false"/>, the "standard" convention is followed, namely the cut is placed along
        /// <i>x</i> > 18381 km, <i>y</i> = 0m.
        /// <see cref="Forward(double, double, double)"/> can be called with any <i>lat</i> and <i>lon</i> then produces the transformation shown
        /// in Lee, Fig 46.
        /// <see cref="Reverse(double, double, double)"/> analytically continues this in the ± <i>x</i> direction.
        /// As a consequence, Reverse may map multiple points to the same geographic location; for example, for <see cref="UTM"/>,
        /// <i>x</i> = 22051449.037349 m, <i>y</i> = −7131237.022729 m and <i>x</i> = 29735142.378357 m, <i>y</i> = 4235043.607933 m both map to
        /// <i>lat</i> = −2°, <i>lon</i> = 88°.
        /// </para>
        /// <para>
        /// With <paramref name="extendp"/> = <see langword="true"/>, the branch cut is moved to the lower left quadrant.
        /// The various symmetries of the transverse Mercator projection can be used to explore the projection on any sheet.
        /// In this mode the domains of <i>lat</i>, <i>lon</i>, <i>x</i>, and <i>y</i> are restricted to
        /// <list type="bullet">
        /// <item>
        /// the union of
        /// <list type="bullet">
        /// <item><i>lat</i> in [0, 90] and <i>lon</i> − <i>lon0</i> in [0, 90]</item>
        /// <item><i>lat</i> in (-90, 0] and <i>lon</i> − <i>lon0</i> in [90 (1 − e), 90]</item>
        /// </list>
        /// </item>
        /// <item>
        /// the union of
        /// <list type="bullet">
        /// <item><i>x</i>/(<i>k0</i> <i>a</i>) in [0, ∞) and <i>y</i>/(<i>k0</i> <i>a</i>) in [0, E(e2)]</item>
        /// <item><i>x</i>/(<i>k0</i> <i>a</i>) in [K(1 − <i>e</i>^2) − E(1 − <i>e</i>^2), ∞) and <i>y</i>/(<i>k0</i> <i>a</i>) in (−∞, 0]</item>
        /// </list>
        /// </item>
        /// </list>
        /// See Sec. 5 of <a href="https://arxiv.org/abs/1002.1417">arXiv:1002.1417</a> for a full discussion of the treatment of the branch cut.
        /// </para>
        /// <para>
        /// The method will work for all ellipsoids used in terrestrial geodesy.
        /// The method cannot be applied directly to the case of a sphere (<i>f</i> = 0) because some the constants characterizing this method diverge in
        /// that limit, and in practice, <i>f</i> should be larger than about <see cref="DBL_EPSILON"/>.
        /// However, <see cref="TransverseMercator"/> treats the sphere exactly.
        /// </para>
        /// </remarks>
        public TransverseMercatorExact(double a, double f, double k0, bool extendp = false)
        {
            _a = a;
            _f = f;
            _k0 = k0;
            _mu = _f * (2 - _f);
            _mv = 1 - _mu;
            _e = Sqrt(_mu);
            _extendp = extendp;
            _Eu = new EllipticFunction(_mu);
            _Ev = new EllipticFunction(_mv);

            if (!(IsFinite(_a) && _a > 0))
                throw new GeographicException("Equatorial radius is not positive");
            if (!(_f > 0))
                throw new GeographicException("Flattening is not positive");
            if (!(_f < 1))
                throw new GeographicException("Polar semi-axis is not positive");
            if (!(IsFinite(_k0) && _k0 > 0))
                throw new GeographicException("Scale is not positive");
        }

        /// <summary>
        /// Initialize a new <see cref="TransverseMercatorExact"/> instance with specified <see cref="IEllipsoid"/> object and central scale factor.
        /// </summary>
        /// <param name="ellipsoid">the <see cref="IEllipsoid"/> object.</param>
        /// <param name="k0">central scale factor.</param>
        /// <param name="extendp">use extended domain.</param>
        /// <remarks>
        /// See <see cref="TransverseMercatorExact(double, double, double, bool)"/> for detailed explanation.
        /// </remarks>
        public TransverseMercatorExact(IEllipsoid ellipsoid, double k0, bool extendp = false)
            : this(ellipsoid.EquatorialRadius, ellipsoid.Flattening, k0, extendp) { }

        #region Private Methods

        private void Zeta(double u, double snu, double cnu, double dnu,
                  double v, double snv, double cnv, double dnv,
                  out double taup, out double lam)
        {
            // Lee 54.17 but write
            // atanh(snu * dnv) = asinh(snu * dnv / sqrt(cnu^2 + _mv * snu^2 * snv^2))
            // atanh(_e * snu / dnv) =
            //         asinh(_e * snu / sqrt(_mu * cnu^2 + _mv * cnv^2))
            // Overflow value s.t. atan(overflow) = pi/2
            const double overflow = 1 / (DBL_EPSILON * DBL_EPSILON);
            double
              d1 = Sqrt(Sq(cnu) + _mv * Sq(snu * snv)),
              d2 = Sqrt(_mu * Sq(cnu) + _mv * Sq(cnv)),
              t1 = (d1 != 0 ? snu * dnv / d1 : (snu < 0 ? -overflow : overflow)),
              t2 = (d2 != 0 ? Sinh(_e * Asinh(_e * snu / d2)) :
                    (snu < 0 ? -overflow : overflow));
            // psi = asinh(t1) - asinh(t2)
            // taup = sinh(psi)
            taup = t1 * Hypot(1, t2) - t2 * Hypot(1, t1);
            lam = (d1 != 0 && d2 != 0) ?
              Atan2(dnu * snv, cnu * cnv) - _e * Atan2(_e * cnu * snv, dnu * cnv) :
              0;
        }

        private void DwdZeta(double u, double snu, double cnu, double dnu,
                     double v, double snv, double cnv, double dnv,
                     out double du, out double dv)
        {
            // Lee 54.21 but write (1 - dnu^2 * snv^2) = (cnv^2 + _mu * snu^2 * snv^2)
            // (see A+S 16.21.4)
            var d = _mv * Sq(Sq(cnv) + _mu * Sq(snu * snv));
            du = cnu * dnu * dnv * (Sq(cnv) - _mu * Sq(snu * snv)) / d;
            dv = -snu * snv * cnv * (Sq(dnu * dnv) + _mu * Sq(cnu)) / d;
        }

        private bool ZetaInv0(double psi, double lam, out double u, out double v)
        {
            bool retval = false;
            if (psi < -_e * PI / 4 &&
                lam > (1 - 2 * _e) * PI / 2 &&
                psi < lam - (1 - _e) * PI / 2)
            {
                // N.B. this branch is normally not taken because psi < 0 is converted
                // psi > 0 by Forward.
                //
                // There's a log singularity at w = w0 = Eu.K() + i * Ev.K(),
                // corresponding to the south pole, where we have, approximately
                //
                //   psi = _e + i * pi/2 - _e * atanh(cos(i * (w - w0)/(1 + _mu/2)))
                //
                // Inverting this gives:
                double
                  psix = 1 - psi / _e,
                  lamx = (PI / 2 - lam) / _e;
                u = Asinh(Sin(lamx) / Hypot(Cos(lamx), Sinh(psix))) *
                  (1 + _mu / 2);
                v = Atan2(Cos(lamx), Sinh(psix)) * (1 + _mu / 2);
                u = _Eu.K() - u;
                v = _Ev.K() - v;
            }
            else if (psi < _e * PI / 2 &&
                     lam > (1 - 2 * _e) * PI / 2)
            {
                // At w = w0 = i * Ev.K(), we have
                //
                //     zeta = zeta0 = i * (1 - _e) * pi/2
                //     zeta' = zeta'' = 0
                //
                // including the next term in the Taylor series gives:
                //
                // zeta = zeta0 - (_mv * _e) / 3 * (w - w0)^3
                //
                // When inverting this, we map arg(w - w0) = [-90, 0] to
                // arg(zeta - zeta0) = [-90, 180]
                double
                  dlam = lam - (1 - _e) * PI / 2,
                  rad = Hypot(psi, dlam),
                  // atan2(dlam-psi, psi+dlam) + 45d gives arg(zeta - zeta0) in range
                  // [-135, 225).  Subtracting 180 (since multiplier is negative) makes
                  // range [-315, 45).  Multiplying by 1/3 (for cube root) gives range
                  // [-105, 15).  In particular the range [-90, 180] in zeta space maps
                  // to [-90, 0] in w space as required.
                  ang = Atan2(dlam - psi, psi + dlam) - 0.75 * PI;
                // Error using this guess is about 0.21 * (rad/e)^(5/3)
                retval = rad < _e * taytol_;
                rad = Cbrt(3 / (_mv * _e) * rad);
                ang /= 3;
                u = rad * Cos(ang);
                v = rad * Sin(ang) + _Ev.K();
            }
            else
            {
                // Use spherical TM, Lee 12.6 -- writing atanh(sin(lam) / cosh(psi)) =
                // asinh(sin(lam) / hypot(cos(lam), sinh(psi))).  This takes care of the
                // log singularity at zeta = Eu.K() (corresponding to the north pole)
                v = Asinh(Sin(lam) / Hypot(Cos(lam), Sinh(psi)));
                u = Atan2(Sinh(psi), Cos(lam));
                // But scale to put 90,0 on the right place
                u *= _Eu.K() / (PI / 2);
                v *= _Eu.K() / (PI / 2);
            }
            return retval;
        }

        private void ZetaInv(double taup, double lam, out double u, out double v)
        {
            double
                psi = Asinh(taup),
                scal = 1 / Hypot(1, taup);
            if (ZetaInv0(psi, lam, out u, out v))
                return;
            var stol2 = tol2_ / Sq(Max(psi, 1));
            // min iterations = 2, max iterations = 6; mean = 4.0
            for (int i = 0, trip = 0; i < numit_ || GEOGRAPHICLIB_PANIC; ++i)
            {
                _Eu.Sncndn(u, out var snu, out var cnu, out var dnu);
                _Ev.Sncndn(v, out var snv, out var cnv, out var dnv);

                Zeta(u, snu, cnu, dnu, v, snv, cnv, dnv, out var tau1, out var lam1);
                DwdZeta(u, snu, cnu, dnu, v, snv, cnv, dnv, out var du1, out var dv1);

                tau1 -= taup;
                lam1 -= lam;
                tau1 *= scal;
                double
                  delu = tau1 * du1 - lam1 * dv1,
                  delv = tau1 * dv1 + lam1 * du1;
                u -= delu;
                v -= delv;
                if (trip != 0)
                    break;
                var delw2 = Sq(delu) + Sq(delv);
                if (!(delw2 >= stol2))
                    ++trip;
            }
        }

        private void Sigma(double u, double snu, double cnu, double dnu,
                   double v, double snv, double cnv, double dnv,
                   out double xi, out double eta)
        {
            // Lee 55.4 writing
            // dnu^2 + dnv^2 - 1 = _mu * cnu^2 + _mv * cnv^2
            var d = _mu * Sq(cnu) + _mv * Sq(cnv);
            xi = _Eu.E(snu, cnu, dnu) - _mu * snu * cnu * dnu / d;
            eta = v - _Ev.E(snv, cnv, dnv) + _mv * snv * cnv * dnv / d;
        }

        private void DwdSigma(double u, double snu, double cnu, double dnu,
                      double v, double snv, double cnv, double dnv,
                      out double du, out double dv)
        {
            // Reciprocal of 55.9: dw/ds = dn(w)^2/_mv, expanding complex dn(w) using
            // A+S 16.21.4
            var d = _mv * Sq(Sq(cnv) + _mu * Sq(snu * snv));
            double
              dnr = dnu * cnv * dnv,
              dni = -_mu * snu * cnu * snv;
            du = (Sq(dnr) - Sq(dni)) / d;
            dv = 2 * dnr * dni / d;
        }

        private bool SigmaInv0(double xi, double eta, out double u, out double v)
        {
            bool retval = false;
            if (eta > 1.25 * _Ev.KE() ||
                (xi < -0.25 * _Eu.E() && xi < eta - _Ev.KE()))
            {
                // sigma as a simple pole at w = w0 = Eu.K() + i * Ev.K() and sigma is
                // approximated by
                //
                // sigma = (Eu.E() + i * Ev.KE()) + 1/(w - w0)
                double
                  x = xi - _Eu.E(),
                  y = eta - _Ev.KE(),
                  r2 = Sq(x) + Sq(y);
                u = _Eu.K() + x / r2;
                v = _Ev.K() - y / r2;
            }
            else if ((eta > 0.75 * _Ev.KE() && xi < 0.25 * _Eu.E())
                     || eta > _Ev.KE())
            {
                // At w = w0 = i * Ev.K(), we have
                //
                //     sigma = sigma0 = i * Ev.KE()
                //     sigma' = sigma'' = 0
                //
                // including the next term in the Taylor series gives:
                //
                // sigma = sigma0 - _mv / 3 * (w - w0)^3
                //
                // When inverting this, we map arg(w - w0) = [-pi/2, -pi/6] to
                // arg(sigma - sigma0) = [-pi/2, pi/2]
                // mapping arg = [-pi/2, -pi/6] to [-pi/2, pi/2]
                double
                  deta = eta - _Ev.KE(),
                  rad = Hypot(xi, deta),
                  // Map the range [-90, 180] in sigma space to [-90, 0] in w space.  See
                  // discussion in zetainv0 on the cut for ang.
                  ang = Atan2(deta - xi, xi + deta) - 0.75 * PI;
                // Error using this guess is about 0.068 * rad^(5/3)
                retval = rad < 2 * taytol_;
                rad = Cbrt(3 / _mv * rad);
                ang /= 3;
                u = rad * Cos(ang);
                v = rad * Sin(ang) + _Ev.K();
            }
            else
            {
                // Else use w = sigma * Eu.K/Eu.E (which is correct in the limit _e -> 0)
                u = xi * _Eu.K() / _Eu.E();
                v = eta * _Eu.K() / _Eu.E();
            }
            return retval;
        }

        private void SigmaInv(double xi, double eta, out double u, out double v)
        {
            if (SigmaInv0(xi, eta, out u, out v))
                return;
            // min iterations = 2, max iterations = 7; mean = 3.9
            for (int i = 0, trip = 0; i < numit_ || GEOGRAPHICLIB_PANIC; ++i)
            {
                _Eu.Sncndn(u, out var snu, out var cnu, out var dnu);
                _Ev.Sncndn(v, out var snv, out var cnv, out var dnv);

                Sigma(u, snu, cnu, dnu, v, snv, cnv, dnv, out var xi1, out var eta1);
                DwdSigma(u, snu, cnu, dnu, v, snv, cnv, dnv, out var du1, out var dv1);

                xi1 -= xi;
                eta1 -= eta;
                double
                  delu = xi1 * du1 - eta1 * dv1,
                  delv = xi1 * dv1 + eta1 * du1;
                u -= delu;
                v -= delv;
                if (trip != 0)
                    break;
                var delw2 = Sq(delu) + Sq(delv);
                if (!(delw2 >= tol2_))
                    ++trip;
            }
        }

        private void Scale(double tau, double lam,
                   double snu, double cnu, double dnu,
                   double snv, double cnv, double dnv,
                   out double gamma, out double k)
        {
            var sec2 = 1 + Sq(tau);    // sec(phi)^2
                                       // Lee 55.12 -- negated for our sign convention.  gamma gives the bearing
                                       // (clockwise from true north) of grid north
            gamma = Atan2(_mv * snu * snv * cnv, cnu * dnu * dnv);
            // Lee 55.13 with nu given by Lee 9.1 -- in sqrt change the numerator
            // from
            //
            //    (1 - snu^2 * dnv^2) to (_mv * snv^2 + cnu^2 * dnv^2)
            //
            // to maintain accuracy near phi = 90 and change the denomintor from
            //
            //    (dnu^2 + dnv^2 - 1) to (_mu * cnu^2 + _mv * cnv^2)
            //
            // to maintain accuracy near phi = 0, lam = 90 * (1 - e).  Similarly
            // rewrite sqrt term in 9.1 as
            //
            //    _mv + _mu * c^2 instead of 1 - _mu * sin(phi)^2
            k = Sqrt(_mv + _mu / sec2) * Sqrt(sec2) *
              Sqrt((_mv * Sq(snv) + Sq(cnu * dnv)) /
                    (_mu * Sq(cnu) + _mv * Sq(cnv)));
        }

        #endregion

        /// <inheritdoc/>
        public double EquatorialRadius => _a;

        /// <inheritdoc/>
        public double Flattening => _f;

        /// <summary>
        /// Gets a value representing central scale for the projection.
        /// This is the value of <i>k0</i> used in the constructor and is the scale on the central meridian.
        /// </summary>
        public double CentralScale => _k0;

        /// <summary>
        /// A global instantiation of <see cref="TransverseMercatorExact"/> with the WGS84 ellipsoid and the UTM scale factor.
        /// However, unlike UTM, no false easting or northing is added.
        /// </summary>
        public static TransverseMercatorExact UTM { get; } = new TransverseMercatorExact(Ellipsoid.WGS84, Constants.UTM_k0);

        /// <summary>
        /// Forward projection, from geographic to transverse Mercator.
        /// </summary>
        /// <param name="lon0">central meridian of the projection (degrees).</param>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <param name="gamma">meridian convergence at point (degrees).</param>
        /// <param name="k">scale of projection at point.</param>
        /// <returns><i>x</i>, easting of point and <i>y</i>, northing of point, in meters.</returns>
        /// <remarks>
        /// No false easting or northing is added. <paramref name="lat"/> should be in the range [−90°, 90°].
        /// </remarks>
        public (double x, double y) Forward(double lon0, double lat, double lon, out double gamma, out double k)
        {
            lat = LatFix(lat);
            lon = AngDiff(lon0, lon);
            // Explicitly enforce the parity
            int
              latsign = (!_extendp && SignBit(lat)) ? -1 : 1,
              lonsign = (!_extendp && SignBit(lon)) ? -1 : 1;
            lon *= lonsign;
            lat *= latsign;
            bool backside = !_extendp && lon > 90;
            if (backside)
            {
                if (lat == 0)
                    latsign = -1;
                lon = 180 - lon;
            }
            double
              lam = lon * Degree,
              tau = Tand(lat);

            // u,v = coordinates for the Thompson TM, Lee 54
            double u, v;
            if (lat == 90)
            {
                u = _Eu.K();
                v = 0;
            }
            else if (lat == 0 && lon == 90 * (1 - _e))
            {
                u = 0;
                v = _Ev.K();
            }
            else
                // tau = tan(phi), taup = sinh(psi)
                ZetaInv(Taupf(tau, _e), lam, out u, out v);

            _Eu.Sncndn(u, out var snu, out var cnu, out var dnu);
            _Ev.Sncndn(v, out var snv, out var cnv, out var dnv);

            Sigma(u, snu, cnu, dnu, v, snv, cnv, dnv, out var xi, out var eta);

            if (backside)
                xi = 2 * _Eu.E() - xi;

            var (y, x) = (xi * _a * _k0 * latsign, eta * _a * _k0 * lonsign);

            if (lat == 90)
            {
                gamma = lon;
                k = 1;
            }
            else
            {
                // Recompute (tau, lam) from (u, v) to improve accuracy of Scale
                Zeta(u, snu, cnu, dnu, v, snv, cnv, dnv,out tau, out lam);
                tau = Tauf(tau, _e);
                Scale(tau, lam, snu, cnu, dnu, snv, cnv, dnv, out gamma, out k);
                gamma /= Degree;
            }
            if (backside)
                gamma = 180 - gamma;
            gamma *= latsign * lonsign;
            k *= _k0;

            return (x, y);
        }

        /// <summary>
        /// Reverse projection, from transverse Mercator to geographic.
        /// </summary>
        /// <param name="lon0">central meridian of the projection (degrees).</param>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <param name="gamma">meridian convergence at point (degrees).</param>
        /// <param name="k">scale of projection at point.</param>
        /// <returns><i>lat</i>, latitude of point and <i>lon</i>, longitude of point, in degrees.</returns>
        /// <remarks>
        /// No false easting or northing is added. The value of <i>lon</i> returned is in the range [−180°, 180°].
        /// </remarks>
        public (double lat, double lon) Reverse(double lon0, double x, double y, out double gamma, out double k)
        {
            // This undoes the steps in Forward.
            double
              xi = y / (_a * _k0),
              eta = x / (_a * _k0);
            // Explicitly enforce the parity
            int
              xisign = (!_extendp && SignBit(xi)) ? -1 : 1,
              etasign = (!_extendp && SignBit(eta)) ? -1 : 1;
            xi *= xisign;
            eta *= etasign;
            bool backside = !_extendp && xi > _Eu.E();
            if (backside)
                xi = 2 * _Eu.E() - xi;

            // u,v = coordinates for the Thompson TM, Lee 54
            double u, v;
            if (xi == 0 && eta == _Ev.KE())
            {
                u = 0;
                v = _Ev.K();
            }
            else
                SigmaInv(xi, eta, out u,out v);

            _Eu.Sncndn(u, out var snu, out var cnu, out var dnu);
            _Ev.Sncndn(v, out var snv, out var cnv, out var dnv);

            double phi, lam, tau, lat, lon;
            if (v != 0 || u != _Eu.K())
            {
                Zeta(u, snu, cnu, dnu, v, snv, cnv, dnv, out tau,out lam);
                tau = Tauf(tau, _e);
                phi = Atan(tau);
                lat = phi / Degree;
                lon = lam / Degree;
                Scale(tau, lam, snu, cnu, dnu, snv, cnv, dnv, out gamma,out k);
                gamma /= Degree;
            }
            else
            {
                lat = 90;
                lon = lam = gamma = 0;
                k = 1;
            }

            if (backside)
                lon = 180 - lon;
            lon *= etasign;
            lon = AngNormalize(lon + AngNormalize(lon0));
            lat *= xisign;
            if (backside)
                gamma = 180 - gamma;
            gamma *= xisign * etasign;
            k *= _k0;

            return (lat, lon);
        }

        /// <summary>
        /// <see cref="Forward(double, double, double, out double, out double)"/> without returning the convergence and scale.
        /// </summary>
        /// <param name="lon0">central meridian of the projection (degrees).</param>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <returns><i>x</i>, easting of point and <i>y</i>, northing of point, in meters.</returns>
        public (double x, double y) Forward(double lon0, double lat, double lon) => Forward(lon0, lat, lon, out _, out _);

        /// <summary>
        /// <see cref="Reverse(double, double, double, out double, out double)"/> without returning the convergence and scale.
        /// </summary>
        /// <param name="lon0">central meridian of the projection (degrees).</param>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <returns><i>lat</i>, latitude of point and <i>lon</i>, longitude of point, in degrees.</returns>
        public (double lat, double lon) Reverse(double lon0, double x, double y) => Reverse(lon0, x, y, out _, out _);

    }
}
