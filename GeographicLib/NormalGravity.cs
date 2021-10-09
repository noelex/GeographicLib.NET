using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using static System.Math;
using static GeographicLib.Macros;
using static GeographicLib.MathEx;

namespace GeographicLib
{
    /// <summary>
    /// The normal gravity of the earth.
    /// </summary>
    /// <remarks>
    /// "Normal" gravity refers to an idealization of the earth which is modeled as an rotating ellipsoid.
    /// The eccentricity of the ellipsoid, the rotation speed, and the distribution of mass within the ellipsoid are such that
    /// the ellipsoid is a "level ellipoid", a surface of constant potential (gravitational plus centrifugal). The acceleration
    /// due to gravity is therefore perpendicular to the surface of the ellipsoid.
    /// <para>
    /// Because the distribution of mass within the ellipsoid is unspecified, only the potential exterior to the ellipsoid is well
    /// defined. In this class, the mass is assumed to be to concentrated on a "focal disc" of radius, (<i>a</i>^2 − <i>b</i>^2)^(1/2),
    /// where a is the equatorial radius of the ellipsoid and b is its polar semi-axis. In the case of an oblate ellipsoid,
    /// the mass is concentrated on a "focal rod" of length 2(<i>b</i>^2 − <i>a</i>^2)^(1/2). As a result the potential is well
    /// defined everywhere.
    /// </para>
    /// <para>
    /// There is a closed solution to this problem which is implemented here.
    /// Series "approximations" are only used to evaluate certain combinations of elementary functions where use of the closed
    /// expression results in a loss of accuracy for small arguments due to cancellation of the leading terms. However these series
    /// include sufficient terms to give full machine precision.
    /// </para>
    /// <para>
    /// Although the formulation used in this class applies to ellipsoids with arbitrary flattening, in practice, its use should be
    /// limited to about <i>b</i>/<i>a</i> ∈ [0.01, 100] or <i>f</i> ∈ [−99, 0.99].
    /// </para>
    /// Definitions:
    /// <list type="bullet">
    /// <item><i>V</i>0, the gravitational contribution to the normal potential;</item>
    /// <item>Φ, the rotational contribution to the normal potential;</item>
    /// <item><i>U</i> = <i>V</i>0 + Φ, the total potential;</item>
    /// <item><b>Γ</b> = ∇<i>V</i>0, the acceleration due to mass of the earth;</item>
    /// <item><b>f</b> = ∇Φ, the centrifugal acceleration;</item>
    /// <item><b>γ</b> = ∇<i>U</i> = <b>Γ</b> + <b>f</b>, the normal acceleration;</item>
    /// <item><i>X</i>, <i>Y</i>, <i>Z</i>, geocentric coordinates;</item>
    /// <item><i>x</i>, <i>y</i>, <i>z</i>, local cartesian coordinates used to denote the east, north and up directions.</item>
    /// </list>
    /// References:
    /// <list type="bullet">
    /// <item>
    /// C. Somigliana,
    /// Teoria generale del campo gravitazionale dell'ellissoide di rotazione,
    /// Mem. Soc. Astron. Ital, 4, 541–599 (1929).
    /// </item>
    /// <item>
    /// W. A. Heiskanen and H. Moritz,
    /// Physical Geodesy (Freeman, San Francisco, 1967),
    /// Secs. 1-19, 2-7, 2-8 (2-9, 2-10), 6-2 (6-3).
    /// </item>
    /// <item>
    /// B. Hofmann-Wellenhof, H. Moritz,
    /// <a href="https://doi.org/10.1007/978-3-211-33545-1">Physical Geodesy</a> (Second edition, Springer, 2006)
    /// </item>
    /// <item>
    /// H. Moritz,
    /// <a href="https://doi.org/10.1007/BF02521480">Geodetic Reference System 1980</a>,
    /// J. Geodesy 54(3), 395-405 (1980)
    /// </item>
    /// </list>
    /// For more information on normal gravity see <a href="https://geographiclib.sourceforge.io/html/normalgravity.html">Normal gravity</a>.
    /// </remarks>
    public class NormalGravity : IEllipsoid
    {
        private const int maxit_ = 20;
        private const double maxe_ = 1 - DBL_EPSILON;
        private static readonly double eps2_ = Sqrt(DBL_EPSILON) / 100;
        private static readonly double lg2eps_ = -Log2(DBL_EPSILON / 2);

        internal readonly double _a, _GM, _omega, _f, _J2, _omega2, _aomega2;
        private readonly double _e2, _ep2, _b, _E, _U0, _gammae, _gammap, _Q0, _k, _fstar;
        private readonly Geocentric _earth;

        /// <summary>
        /// Constructor for the normal gravity.
        /// </summary>
        /// <param name="a">equatorial radius (meters).</param>
        /// <param name="GM">mass constant of the ellipsoid (meters^3/seconds^2);
        /// this is the product of <i>G</i> the gravitational constant and <i>M</i> the mass of the earth
        /// (usually including the mass of the earth's atmosphere).</param>
        /// <param name="omega">the angular velocity (rad s^−1).</param>
        /// <param name="f_J2">either the flattening of the ellipsoid <i>f</i> or the the dynamical form factor <i>J2</i>.</param>
        /// <param name="geometricp">
        /// if <see langword="true"/> (the default), then <paramref name="f_J2"/> denotes the flattening,
        /// else it denotes the dynamical form factor <i>J2</i>.
        /// </param>
        /// <remarks>
        /// The shape of the ellipsoid can be given in one of two ways:
        /// <list type="bullet">
        /// <item>geometrically (<paramref name="geometricp"/> = <see langword="true"/>), the ellipsoid is defined by
        /// the flattening <i>f</i> = (<i>a</i> − <i>b</i>) / <i>a</i>,
        /// where <i>a</i> and <i>b</i> are the equatorial radius and the polar semi-axis.
        /// The parameters should obey <i>a</i> > 0, <i>f</i> &lt; 1.
        /// There are no restrictions on <paramref name="GM"/> or <paramref name="omega"/>, in particular, <paramref name="GM"/> need not be positive.</item>
        /// <item>physically (<paramref name="geometricp"/> = <see langword="false"/>),
        /// the ellipsoid is defined by the dynamical form factor <i>J2</i> = (<i>C</i> − <i>A</i>) / <i>Ma</i>^2,
        /// where <i>A</i> and <i>C</i> are the equatorial and polar moments of inertia and <i>M</i> is the mass of the earth.
        /// The parameters should obey <paramref name="a"/> > 0,
        /// <paramref name="GM"/> > 0 and <i>J2</i> &lt; 1/3 − (<paramref name="omega"/>^2 <paramref name="a"/>^3/<paramref name="GM"/>) 8/(45π).
        /// There is no restriction on omega.</item>
        /// </list>
        /// </remarks>
        public NormalGravity(double a, double GM, double omega, double f_J2,
                  bool geometricp = true)
        {
            _a = a;
            if (!(IsFinite(_a) && _a > 0))
                throw new GeographicException("Equatorial radius is not positive");
            _GM = GM;
            if (!IsFinite(_GM))
                throw new GeographicException("Gravitational constant is not finite");
            _omega = omega;
            _omega2 = Sq(_omega);
            _aomega2 = Sq(_omega * _a);
            if (!(IsFinite(_omega2) && IsFinite(_aomega2)))
                throw new GeographicException("Rotation velocity is not finite");
            _f = geometricp ? f_J2 : J2ToFlattening(_a, _GM, _omega, f_J2);
            _b = _a * (1 - _f);
            if (!(IsFinite(_b) && _b > 0))
                throw new GeographicException("Polar semi-axis is not positive");
            _J2 = geometricp ? FlatteningToJ2(_a, _GM, _omega, f_J2) : f_J2;
            _e2 = _f * (2 - _f);
            _ep2 = _e2 / (1 - _e2);
            var ex2 = _f < 0 ? -_e2 : _ep2;
            _Q0 = Qf(ex2, _f < 0);
            _earth = new Geocentric(_a, _f);
            _E = _a * Sqrt(Abs(_e2));   // H+M, Eq 2-54
                                        // H+M, Eq 2-61
            _U0 = _GM * Atanzz(ex2, _f < 0) / _b + _aomega2 / 3;
            var P = Hf(ex2, _f < 0) / (6 * _Q0);
            // H+M, Eq 2-73
            _gammae = _GM / (_a * _b) - (1 + P) * _a * _omega2;
            // H+M, Eq 2-74
            _gammap = _GM / (_a * _a) + 2 * P * _b * _omega2;
            // k = gammae * (b * gammap / (a * gammae) - 1)
            //   = (b * gammap - a * gammae) / a
            _k = -_e2 * _GM / (_a * _b) +
              _omega2 * (P * (_a + 2 * _b * (1 - _f)) + _a);
            // f* = (gammap - gammae) / gammae
            _fstar = (-_f * _GM / (_a * _b) + _omega2 * (P * (_a + 2 * _b) + _a)) /
              _gammae;
        }

        /// <summary>
        /// Evaluate the gravity on the surface of the ellipsoid.
        /// </summary>
        /// <param name="lat">the geographic latitude (degrees).</param>
        /// <returns>γ, the acceleration due to gravity, positive downwards (m s^−2).</returns>
        /// <remarks>
        /// Due to the axial symmetry of the ellipsoid, the result is independent of the value of the longitude.
        /// This acceleration is perpendicular to the surface of the ellipsoid. It includes the effects of the earth's rotation.
        /// </remarks>
        public double SurfaceGravity(double lat)
        {
            var sphi = Sind(LatFix(lat));
            // H+M, Eq 2-78
            return (_gammae + _k * Sq(sphi)) / Sqrt(1 - _e2 * Sq(sphi));
        }

        /// <summary>
        /// Evaluate the gravity at an arbitrary point above (or below) the ellipsoid.
        /// </summary>
        /// <param name="lat">the geographic latitude (degrees).</param>
        /// <param name="h">the height above the ellipsoid (meters).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>U</i>, the corresponding normal potential (m^2 s^−2).</item>
        /// <item><i>gammay</i>, the northerly component of the acceleration (m s^−2).</item>
        /// <item><i>gammaz</i>, the upward component of the acceleration (m s−^2); this is usually negative.</item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// Due to the axial symmetry of the ellipsoid, the result is independent of the value of the longitude and the easterly
        /// component of the acceleration vanishes, <i>gammax</i> = 0. The function includes the effects of the earth's rotation.
        /// When <i>h</i> = 0, this function gives <i>gammay</i> = 0 and the returned value matches that of
        /// <see cref="SurfaceGravity(double)"/>.
        /// </remarks>
        public (double U, double gammay, double gammaz) Gravity(double lat, double h)
        {
            Span<double> M = stackalloc double[Geocentric.dim2_];
            var (X, Y, Z) = _earth.IntForward(lat, 0, h, M);
            var (Ures, gammaX, gammaY, gammaZ) = U(X, Y, Z);
            // gammax = M[0] * gammaX + M[3] * gammaY + M[6] * gammaZ;
            var gammay = M[1] * gammaX + M[4] * gammaY + M[7] * gammaZ;
            var gammaz = M[2] * gammaX + M[5] * gammaY + M[8] * gammaZ;
            return (Ures, gammay, gammaz);
        }

        /// <summary>
        /// Evaluate the components of the acceleration due to gravity and the centrifugal acceleration in geocentric coordinates.
        /// </summary>
        /// <param name="X"><i>X</i> component of geocentric coordinate of point (meters).</param>
        /// <param name="Y"><i>Y</i> component of geocentric coordinate of point (meters).</param>
        /// <param name="Z"><i>Z</i> component of geocentric coordinate of point (meters).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>U</i> = <i>V</i>0 + Φ, the sum of the gravitational and centrifugal potentials (m^2 s^−2).</item>
        /// <item><i>gammaX</i>, the <i>X</i> component of the acceleration (m s^−2).</item>
        /// <item><i>gammaY</i>, the <i>Y</i> component of the acceleration (m s^−2).</item>
        /// <item><i>gammaZ</i>, the <i>Z</i> component of the acceleration (m s^−2).</item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// The acceleration given by <b>γ</b> = ∇<i>U</i> = ∇<i>V</i>0 + ∇Φ = <b>Γ</b> + <b>f</b>.
        /// </remarks>
        public (double U, double gammaX, double gammaY, double gammaZ) U(double X, double Y, double Z)
        {
            var (Ures, gammaX, gammaY, gammaZ) = V0(X, Y, Z);
            var (phi, fX, fY) = Phi(X, Y);
            Ures += phi;
            gammaX += fX;
            gammaY += fY;
            return (Ures, gammaX, gammaY, gammaZ);
        }

        /// <summary>
        /// Evaluate the components of the acceleration due to the gravitational force in geocentric coordinates.
        /// </summary>
        /// <param name="X"><i>X</i> component of geocentric coordinate of point (meters).</param>
        /// <param name="Y"><i>Y</i> component of geocentric coordinate of point (meters).</param>
        /// <param name="Z"><i>Z</i> component of geocentric coordinate of point (meters).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>V</i>0, the gravitational potential (m^2 s^−2).</item>
        /// <item><i>GammaX</i>, the <i>X</i> component of the acceleration due to the gravitational force (m s^−2).</item>
        /// <item><i>GammaY</i>, the <i>Y</i> component of the acceleration due to the gravitational force (m s^−2).</item>
        /// <item><i>GammaZ</i>, the <i>Z</i> component of the acceleration due to the gravitational force (m s^−2).</item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// This function excludes the centrifugal acceleration and is appropriate to use for space applications.
        /// In terrestrial applications, the function <see cref="U(double, double, double)"/>
        /// (which includes this effect) should usually be used.
        /// </remarks>
        public (double V0, double GammaX, double GammaY, double GammaZ) V0(double X, double Y, double Z)
        {
            // See H+M, Sec 6-2
            double
              p = Hypot(X, Y),
              clam = p != 0 ? X / p : 1,
              slam = p != 0 ? Y / p : 0,
              r = Hypot(p, Z);
            if (_f < 0) Swap(ref p, ref Z);
            double
              Q = Sq(r) - Sq(_E),
              t2 = Sq(2 * _E * Z),
              disc = Sqrt(Sq(Q) + t2),
              // This is H+M, Eq 6-8a, but generalized to deal with Q negative
              // accurately.
              u = Sqrt((Q >= 0 ? (Q + disc) : t2 / (disc - Q)) / 2),
              uE = Hypot(u, _E),
              // H+M, Eq 6-8b
              sbet = u != 0 ? Z * uE : CopySign(Sqrt(-Q), Z),
              cbet = u != 0 ? p * u : p,
              s = Hypot(cbet, sbet);
            sbet = s != 0 ? sbet / s : 1;
            cbet = s != 0 ? cbet / s : 0;
            double
              z = _E / u,
              z2 = Sq(z),
              den = Hypot(u, _E * sbet);
            if (_f < 0)
            {
                Swap(ref sbet, ref cbet);
                Swap(ref u, ref uE);
            }
            double
              invw = uE / den,          // H+M, Eq 2-63
              bu = _b / (u != 0 || _f < 0 ? u : _E),
              // Qf(z2->inf, false) = pi/(4*z^3)
              q = ((u != 0 || _f < 0 ? Qf(z2, _f < 0) : PI / 4) / _Q0) *
                bu * Sq(bu),
              qp = _b * Sq(bu) * (u != 0 || _f < 0 ? Hf(z2, _f < 0) : 2) / _Q0,
              ang = (Sq(sbet) - 1 / 3d) / 2,
              // H+M, Eqs 2-62 + 6-9, but omitting last (rotational) term.
              Vres = _GM * (u != 0 || _f < 0 ?
                            Atanzz(z2, _f < 0) / u :
                            PI / (2 * _E)) + _aomega2 * q * ang,
              // H+M, Eq 6-10
              gamu = -(_GM + (_aomega2 * qp * ang)) * invw / Sq(uE),
              gamb = _aomega2 * q * sbet * cbet * invw / uE,
              t = u * invw / uE,
              gamp = t * cbet * gamu - invw * sbet * gamb;
            // H+M, Eq 6-12
            var GammaX = gamp * clam;
            var GammaY = gamp * slam;
            var GammaZ = invw * sbet * gamu + t * cbet * gamb;
            return (Vres, GammaX, GammaY, GammaZ);
        }

        /// <summary>
        /// Evaluate the centrifugal acceleration in geocentric coordinates.
        /// </summary>
        /// <param name="X"><i>X</i> component of geocentric coordinate of point (meters).</param>
        /// <param name="Y"><i>Y</i> component of geocentric coordinate of point (meters).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item>Φ, the centrifugal potential (m^2 s^−2).</item>
        /// <item><i>fX</i>, the <i>X</i> component of the centrifugal acceleration (m s^−2).</item>
        /// <item><i>fY</i>, the <i>Y</i> component of the centrifugal acceleration (m s^−2).</item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// Φ is independent of <i>Z</i>, thus <i>fZ</i> = 0.
        /// This function <see cref="U(double, double, double)"/> sums the results of
        /// <see cref="V0(double, double, double)"/> and <see cref="Phi(double, double)"/>.
        /// </remarks>
        public (double phi, double fX, double fY) Phi(double X, double Y)
        {
            var fX = _omega2 * X;
            var fY = _omega2 * Y;
            // N.B. fZ = 0;
            return (_omega2 * (Sq(X) + Sq(Y)) / 2, fX, fY);
        }

        /// <summary>
        /// Gets the dynamical form factor of the ellipsoid.
        /// </summary>
        /// <param name="n"></param>
        /// <returns><i>Jn</i>, the dynamical form factors of the ellipsoid.</returns>
        /// <remarks>
        /// If <paramref name="n"/> = 2 (the default), this is the value of <i>J2</i> used in the constructor.
        /// Otherwise it is the zonal coefficient of the Legendre harmonic sum of the normal gravitational potential.
        /// Note that <i>Jn</i> = 0 if <paramref name="n"/> is odd. In most gravity applications,
        /// fully normalized Legendre functions are used and the corresponding coefficient is
        /// <i>Cn</i>0 = −<i>Jn</i> / sqrt(2 <paramref name="n"/> + 1).
        /// </remarks>
        public double DynamicalFormFactor(int n = 2) => n == 2 ? _J2 : Jn(n);

        /// <inheritdoc/>
        public double EquatorialRadius => _a;

        /// <inheritdoc/>
        public double Flattening => _f;

        /// <summary>
        /// Gets a value representing <i>GM</i>, the mass constant of the ellipsoid (m^3 s^−2).
        /// This is the value used in the constructor.
        /// </summary>
        public double MassConstant => _GM;

        /// <summary>
        /// Gets a value representing ω, the angular velocity of the ellipsoid (rad s^−1).
        /// This is the value used in the constructor.
        /// </summary>
        public double AngularVelocity => _omega;

        /// <summary>
        /// Gets a value representing γe, the normal gravity at equator (m s^−2).
        /// </summary>
        public double EquatorialGravity => _gammae;

        /// <summary>
        /// Gets a value representing γp, the normal gravity at poles (m s^−2).
        /// </summary>
        public double PolarGravity => _gammap;

        /// <summary>
        /// Gets a value representing <i>f</i> *, the gravity flattening (γp − γe) / γe.
        /// </summary>
        public double GravityFlattening => _fstar;

        /// <summary>
        /// Gets a value representing <i>U</i>0, the constant normal potential for the surface of the ellipsoid (m^2 s^−2).
        /// </summary>
        public double SurfacePotential => _U0;

        /// <summary>
        /// Gets a value representing the <see cref="Geocentric"/> object used by this instance.
        /// </summary>
        public Geocentric Earth => _earth;

        /// <summary>
        /// A global instantiation of <see cref="NormalGravity"/> for the WGS84 ellipsoid.
        /// </summary>
        public static NormalGravity WGS84 { get; } =
                new NormalGravity(
                    Constants.WGS84_a,
                    Constants.WGS84_GM,
                    Constants.WGS84_omega,
                    Constants.WGS84_f,
                    true
                );

        /// <summary>
        /// A global instantiation of <see cref="NormalGravity"/> for the GRS80 ellipsoid.
        /// </summary>
        public static NormalGravity GRS80 { get; } =
                new NormalGravity(
                    Constants.GRS80_a,
                    Constants.GRS80_GM,
                    Constants.GRS80_omega,
                    Constants.GRS80_J2,
                    false
                );

        /// <summary>
        /// Compute the flattening from the dynamical form factor.
        /// </summary>
        /// <param name="a">equatorial radius (meters).</param>
        /// <param name="GM">mass constant of the ellipsoid (meters^3/seconds^2);
        /// this is the product of <i>G</i> the gravitational constant and <i>M</i> the mass of the earth
        /// (usually including the mass of the earth's atmosphere).</param>
        /// <param name="omega">the angular velocity (rad s^−1).</param>
        /// <param name="J2">the dynamical form factor.</param>
        /// <returns><i>f</i>, the flattening of the ellipsoid.</returns>
        /// <remarks>
        /// This routine requires <i>a</i> > 0, <i>GM</i> > 0, <i>J2</i> &lt; 1/3 − <i>omega</i>^2 <i>a</i>^3/<i>GM</i> 8/(45π).
        /// A <see cref="double.NaN"/> is returned if these conditions do not hold. The restriction to positive <i>GM</i> is made
        /// because for negative <i>GM</i> two solutions are possible.
        /// </remarks>
        public static double J2ToFlattening(double a, double GM, double omega, double J2)
        {
            // Solve
            //   f = e^2 * (1 -  K * e/q0) - 3 * J2 = 0
            // for e^2 using Newton's method
            double
              K = 2 * Sq(a * omega) * a / (15 * GM),
              J0 = (1 - 4 * K / PI) / 3;
            if (!(GM > 0 && IsFinite(K) && K >= 0))
                return double.NaN;
            if (!(IsFinite(J2) && J2 <= J0)) return double.NaN;
            if (J2 == J0) return 1;
            // Solve e2 - f1 * f2 * K / Q0 - 3 * J2 = 0 for J2 close to J0;
            // subst e2 = ep2/(1+ep2), f2 = 1/(1+ep2), f1 = 1/sqrt(1+ep2), J2 = J0-dJ2,
            // Q0 = pi/(4*z^3) - 2/z^4 + (3*pi)/(4*z^5), z = sqrt(ep2), and balance two
            // leading terms to give
            double
              ep2 = Max(Sq(32 * K / (3 * Sq(PI) * (J0 - J2))),
                        -maxe_),
              e2 = Min(ep2 / (1 + ep2), maxe_);
            for (int j = 0; j < maxit_ || GEOGRAPHICLIB_PANIC; ++j)
            {
                double
                  e2a = e2, ep2a = ep2,
                  f2 = 1 - e2,            // (1 - f)^2
                  f1 = Sqrt(f2),          // (1 - f)
                  Q0 = Qf(e2 < 0 ? -e2 : ep2, e2 < 0),
                  h = e2 - f1 * f2 * K / Q0 - 3 * J2,
                  dh = 1 - 3 * f1 * K * QH3f(e2 < 0 ? -e2 : ep2, e2 < 0) /
                               (2 * Sq(Q0));
                e2 = Min(e2a - h / dh, maxe_);
                ep2 = Max(e2 / (1 - e2), -maxe_);
                if (Abs(h) < eps2_ || e2 == e2a || ep2 == ep2a)
                    break;
            }
            return e2 / (1 + Sqrt(1 - e2));
        }

        /// <summary>
        /// Compute the dynamical form factor from the flattening.
        /// </summary>
        /// <param name="a">equatorial radius (meters).</param>
        /// <param name="GM">mass constant of the ellipsoid (meters^3/seconds^2);
        /// this is the product of <i>G</i> the gravitational constant and <i>M</i> the mass of the earth
        /// (usually including the mass of the earth's atmosphere).</param>
        /// <param name="omega">the angular velocity (rad s^−1).</param>
        /// <param name="f">the flattening of the ellipsoid.</param>
        /// <returns><i>J2</i>, the dynamical form factor.</returns>
        /// <remarks>
        /// This routine requires <i>a</i> > 0, <i>GM</i> ≠ 0, <i>f</i> &lt; 1. The values of these parameters are not checked.
        /// </remarks>
        public static double FlatteningToJ2(double a, double GM, double omega, double f)
        {
            double
              K = 2 * Sq(a * omega) * a / (15 * GM),
              f1 = 1 - f,
              f2 = Sq(f1),
              e2 = f * (2 - f);
            // H+M, Eq 2-90 + 2-92'
            return (e2 - K * f1 * f2 / Qf(f < 0 ? -e2 : e2 / f2, f < 0)) / 3;
        }

        private static double Atanzz(double x, bool alt)
        {
            // This routine obeys the identity
            //   atanzz(x, alt) = atanzz(-x/(1+x), !alt)
            //
            // Require x >= -1.  Best to call with alt, s.t. x >= 0; this results in
            // a call to atan, instead of asin, or to asinh, instead of atanh.
            var z = Sqrt(Abs(x));
            return x == 0 ? 1 :
              (alt ?
               (!(x < 0) ? Asinh(z) : Asin(z)) / Sqrt(Abs(x) / (1 + x)) :
               (!(x < 0) ? Atan(z) : Atanh(z)) / z);
        }

        private static double Atan7Series(double x)
        {
            Frexp(x, out var e);
            e = Max(-e, 1);             // Here's where abs(x) < 1/2 is assumed
                                        // x = [0.5,1) * 2^(-e)
                                        // estimate n s.t. x^n/n < 1/7 * epsilon/2
                                        // a stronger condition is x^n < epsilon/2
                                        // taking log2 of both sides, a stronger condition is n*(-e) < -lg2eps;
                                        // or n*e > lg2eps or n > ceiling(lg2eps/e)
            int n = x == 0 ? 1 : (int)Ceiling(lg2eps_ / e);
            double v = 0;
            while (n-- != 0)                 // iterating from n-1 down to 0
                v = -x * v - 1 / (2 * n + 7);
            return v;
        }
        private static double Atan5Series(double x) =>
            // Compute Taylor series approximations to
            //   (atan(z)-(z-z^3/3))/z^5,
            // z = sqrt(x)
            // require abs(x) < 1/2, but better to restrict calls to abs(x) < 1/4
            1 / 5d + x * Atan7Series(x);

        private static double Qf(double x, bool alt)
        {
            // Compute
            //   Q(z) = (((1 + 3/z^2) * atan(z) - 3/z)/2) / z^3
            //        = q(z)/z^3 with q(z) defined by H+M, Eq 2-57 with z = E/u
            //   z = sqrt(x)
            var y = alt ? -x / (1 + x) : x;
            return !(4 * Abs(y) < 1) ?  // Backwards test to allow NaNs through
              ((1 + 3 / y) * Atanzz(x, alt) - 3 / y) / (2 * y) :
              (3 * (3 + y) * Atan5Series(y) - 1) / 6;
        }

        private static double Hf(double x, bool alt)
        {
            // z = sqrt(x)
            // Compute
            //   H(z) = (3*Q(z)+z*diff(Q(z),z))*(1+z^2)
            //        = (3 * (1 + 1/z^2) * (1 - atan(z)/z) - 1) / z^2
            //        = q'(z)/z^2, with q'(z) defined by H+M, Eq 2-67, with z = E/u
            var y = alt ? -x / (1 + x) : x;
            return !(4 * Abs(y) < 1) ?  // Backwards test to allow NaNs through
              (3 * (1 + 1 / y) * (1 - Atanzz(x, alt)) - 1) / y :
              1 - 3 * (1 + y) * Atan5Series(y);
        }

        private static double QH3f(double x, bool alt)
        {
            // z = sqrt(x)
            // (Q(z) - H(z)/3) / z^2
            //   = - (1+z^2)/(3*z) * d(Q(z))/dz - Q(z)
            //   = ((15+9*z^2)*atan(z)-4*z^3-15*z)/(6*z^7)
            //   = ((25+15*z^2)*atan7+3)/10
            var y = alt ? -x / (1 + x) : x;
            return !(4 * Abs(y) < 1) ? // Backwards test to allow NaNs through
              ((9 + 15 / y) * Atanzz(x, alt) - 4 - 15 / y) / (6 * Sq(y)) :
              ((25 + 15 * y) * Atan7Series(y) + 3) / 10;
        }

        internal double Jn(int n)
        {
            // Note Jn(0) = -1; Jn(2) = _J2; Jn(odd) = 0
            if ((n & 1) != 0 || n < 0)
                return 0;
            n /= 2;
            double e2n = 1;            // Perhaps this should just be e2n = pow(-_e2, n);
            for (int j = n; j-- != 0;)
                e2n *= -_e2;
            return                      // H+M, Eq 2-92
              -3 * e2n * ((1 - n) + 5 * n * _J2 / _e2) / ((2 * n + 1) * (2 * n + 3));
        }
    }
}
