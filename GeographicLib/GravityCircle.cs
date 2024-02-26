using GeographicLib.SphericalHarmonics;
using System;

using static GeographicLib.MathEx;

namespace GeographicLib
{
    /// <summary>
    /// Gravity on a circle of latitude.
    /// </summary>
    /// <remarks>
    /// Evaluate the earth's gravity field on a circle of constant height and latitude.
    /// This uses a <see cref="CircularEngine"/> to pre-evaluate the inner sum of the spherical harmonic sum,
    /// allowing the values of the field at several different longitudes to be evaluated rapidly.
    /// <para>
    /// Use <see cref="GravityModel.Circle"/> to create a <see cref="GravityCircle"/> object.
    /// (The constructor for this class is private.)
    /// </para>
    /// <para>
    /// See <a href="https://geographiclib.sourceforge.io/html/gravity.html#gravityparallel">
    /// Geoid heights on a multi-processor system</a> for an example of using
    /// <see cref="GravityCircle"/> (together with OpenMP) to speed up the computation of geoid heights.
    /// </para>
    /// </remarks>
    public class GravityCircle : IEllipsoid
    {
        private readonly GravityFlags _caps;
        private readonly double _a, _f, _lat, _h, _Z, _Px, _invR, _cpsi, _spsi,
          _cphi, _sphi, _amodel, _GMmodel, _dzonal0,
          _corrmult, _gamma0, _gamma, _frot;
        private readonly CircularEngine _gravitational, _disturbing, _correction;

        internal GravityCircle(GravityFlags caps, double a, double f, double lat, double h,
              double Z, double P, double cphi, double sphi,
              double amodel, double GMmodel,
              double dzonal0, double corrmult,
              double gamma0, double gamma, double frot,
              CircularEngine gravitational,
              CircularEngine disturbing,
              CircularEngine correction)
        {
            _caps = caps;
            _a = a;
            _f = f;
            _lat = LatFix(lat);
            _h = h;
            _Z = Z;
            _Px = P;
            _invR = 1 / Hypot(_Px, _Z);
            _cpsi = _Px * _invR;
            _spsi = _Z * _invR;
            _cphi = cphi;
            _sphi = sphi;
            _amodel = amodel;
            _GMmodel = GMmodel;
            _dzonal0 = dzonal0;
            _corrmult = corrmult;
            _gamma0 = gamma0;
            _gamma = gamma;
            _frot = frot;
            _gravitational = gravitational;
            _disturbing = disturbing;
            _correction = correction;
        }

        /// <inheritdoc/>
        public double EquatorialRadius => _a;

        /// <inheritdoc/>
        public double Flattening => _f;

        /// <summary>
        /// Gets a value representing the latitude of the circle (degrees).
        /// </summary>
        public double Latitude => _lat;

        /// <summary>
        /// Gets a value representing the height of the circle (meters).
        /// </summary>
        public double Height => _h;

        /// <summary>
        /// Gets a value representing the computational capabilities that this object was constructed with.
        /// </summary>
        public GravityFlags Capabilities => _caps;

        /// <summary>
        /// Check whether current <see cref="GravityCircle"/> has specified capabilities.
        /// </summary>
        /// <param name="testcaps">a set of bitor'ed <see cref="GravityFlags"/> values.</param>
        /// <returns><see langword="true"/> if the <see cref="GravityCircle"/> object has all these capabilities.</returns>
        public bool HasCapabilities(GravityFlags testcaps) => _caps.HasFlag(testcaps);

        /// <summary>
        /// Evaluate the gravity.
        /// </summary>
        /// <param name="lon">the geographic longitude (degrees).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>W</i>, the sum of the gravitational and centrifugal potentials (m^2 s^−2).</item>
        /// <item><i>gx</i>, the easterly component of the acceleration (m s^−2).</item>
        /// <item><i>gy</i>, the northerly component of the acceleration (m s^−2).</item>
        /// <item><i>gz</i>, the upward component of the acceleration (m s^−2); this is usually negative.</item>
        /// </list>
        /// </returns>
        /// <remarks>The function includes the effects of the earth's rotation.</remarks>
        public (double W, double gx, double gy, double gz) Gravity(double lon)
        {
            Span<double> M = stackalloc double[Geocentric.dim2_];
            SinCosd(lon, out var slam, out var clam);
            var (Wres, gx, gy, gz) = W(slam, clam);
            Geocentric.Rotation(_sphi, _cphi, slam, clam, M);
            Geocentric.Unrotate(M, gx, gy, gz, out gx, out gy, out gz);
            return (Wres, gx, gy, gz);
        }

        /// <summary>
        /// Evaluate the gravity disturbance vector.
        /// </summary>
        /// <param name="lon">the geographic longitude (degrees).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>T</i>, the corresponding disturbing potential (m2 s−2).</item>
        /// <item><i>deltax</i>, the easterly component of the disturbance vector (m s^−2).</item>
        /// <item><i>deltay</i>, the northerly component of the disturbance vector (m s^−2).</item>
        /// <item><i>deltaz</i>, the upward component of the disturbance vector (m s^−2).</item>
        /// </list>
        /// </returns>
        public (double T, double deltax, double deltay, double deltaz) Disturbance(double lon)
        {
            Span<double> M = stackalloc double[Geocentric.dim2_];
            SinCosd(lon, out var slam, out var clam);
            var Tres = InternalT(slam, clam, out var deltax, out var deltay, out var deltaz, true, true);
            Geocentric.Rotation(_sphi, _cphi, slam, clam, M);
            Geocentric.Unrotate(M, deltax, deltay, deltaz, out deltax, out deltay, out deltaz);
            return (Tres, deltax, deltay, deltaz);
        }

        /// <summary>
        /// Evaluate the geoid height.
        /// </summary>
        /// <param name="lon">the geographic longitude (degrees).</param>
        /// <returns><i>N</i>, the height of the geoid above the reference ellipsoid (meters).</returns>
        /// <remarks>
        /// Some approximations are made in computing the geoid height so that the results of the NGA codes are reproduced accurately.
        /// Details are given in
        /// <a href="https://geographiclib.sourceforge.io/html/gravity.html#gravitygeoid">Details of the geoid height and anomaly calculations</a>.
        /// </remarks>
        public double GeoidHeight(double lon)
        {
            if (!_caps.HasFlag(GravityFlags.GeoidHeight))
                return double.NaN;

            SinCosd(lon, out var slam, out var clam);
            var T = InternalT(slam, clam, out _, out _, out _, false, false);
            var correction = _corrmult * _correction.Evaluate(slam, clam);
            return T / _gamma0 + correction;
        }

        /// <summary>
        /// Evaluate the components of the gravity anomaly vector using the spherical approximation.
        /// </summary>
        /// <param name="lon">the geographic longitude (degrees).</param>   w
        /// <returns>
        /// <list type="bullet">
        /// <item><i>Dg01</i>, the gravity anomaly (m s^−2).</item>
        /// <item><i>xi</i>, the northerly component of the deflection of the vertical (degrees).</item>
        /// <item><i>eta</i>, the easterly component of the deflection of the vertical (degrees).</item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// The spherical approximation (see Heiskanen and Moritz, Sec 2-14) is used so that the results of the NGA codes are
        /// reproduced accurately. approximations used here.
        /// Details are given in
        /// <a href="https://geographiclib.sourceforge.io/html/gravity.html#gravitygeoid">Details of the geoid height and anomaly calculations</a>.
        /// </remarks>
        public (double Dg01, double xi, double eta) SphericalAnomaly(double lon)
        {
            if (!_caps.HasFlag(GravityFlags.SphericalAnomaly))
            {
                return (double.NaN, double.NaN, double.NaN);
            }

            SinCosd(lon, out var slam, out var clam);
            double
              deltax, deltay, deltaz,
              T = InternalT(slam, clam, out deltax, out deltay, out deltaz, true, false);
            // Rotate cartesian into spherical coordinates
            Span<double> MC = stackalloc double[Geocentric.dim2_];
            Geocentric.Rotation(_spsi, _cpsi, slam, clam, MC);
            Geocentric.Unrotate(MC, deltax, deltay, deltaz, out deltax, out deltay, out deltaz);
            // H+M, Eq 2-151c
            var Dg01 = -deltaz - 2 * T * _invR;
            var xi = -(deltay / _gamma) / Degree;
            var eta = -(deltax / _gamma) / Degree;

            return (Dg01, xi, eta);
        }

        /// <summary>
        /// Evaluate the components of the acceleration due to gravity and the centrifugal acceleration in geocentric coordinates.
        /// </summary>
        /// <param name="lon">the geographic longitude (degrees).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>W</i> = <i>V</i> + Φ, the sum of the gravitational and centrifugal potentials (m^2 s^−2).</item>
        /// <item><i>gX</i>, the <i>X</i> component of the acceleration (m s^−2).</item>
        /// <item><i>gY</i>, the <i>Y</i> component of the acceleration (m s^−2).</item>
        /// <item><i>gZ</i>, the <i>Z</i> component of the acceleration (m s^−2).</item>
        /// </list>
        /// </returns>
        public (double W, double gX, double gY, double gZ) W(double lon)
        {
            SinCosd(lon, out var slam, out var clam);
            return W(slam, clam);
        }

        /// <summary>
        /// Evaluate the components of the acceleration due to gravity in geocentric coordinates.
        /// </summary>
        /// <param name="lon">the geographic longitude (degrees).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>V</i> = <i>W</i> - Φ, the gravitational potential (m^2 s^−2).</item>
        /// <item><i>GX</i>, the <i>X</i> component of the acceleration (m s^−2).</item>
        /// <item><i>GY</i>, the <i>Y</i> component of the acceleration (m s^−2).</item>
        /// <item><i>GZ</i>, the <i>Z</i> component of the acceleration (m s^−2).</item>
        /// </list>
        /// </returns>
        public (double V, double GX, double GY, double GZ) V(double lon)
        {
            SinCosd(lon, out var slam, out var clam);
            return V(slam, clam);
        }

        /// <summary>
        /// Evaluate the components of the gravity disturbance in geocentric coordinates.
        /// </summary>
        /// <param name="lon">the geographic longitude (degrees).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>T</i> = <i>W</i> - <i>U</i>, the disturbing potential (also called the anomalous potential) (m^2 s^−2).</item>
        /// <item><i>deltaX</i>, the <i>X</i> component of the gravity disturbance (m s^−2).</item>
        /// <item><i>deltaY</i>, the <i>Y</i> component of the gravity disturbance (m s^−2).</item>
        /// <item><i>deltaZ</i>, the <i>Z</i> component of the gravity disturbance (m s^−2).</item>
        /// </list>
        /// </returns>
        public (double T, double deltaX, double deltaY, double deltaZ) Td(double lon)
        {
            SinCosd(lon, out var slam, out var clam);
            return (InternalT(slam, clam, out var deltaX, out var deltaY, out var deltaZ, true, true),
                    deltaX, deltaY, deltaZ);
        }

        /// <summary>
        /// Evaluate disturbing potential in geocentric coordinates.
        /// </summary>
        /// <param name="lon">the geographic longitude (degrees).</param>
        /// <returns>
        /// <i>T</i> = <i>W</i> - <i>U</i> the disturbing potential (also called the anomalous potential) (m^2 s^−2).
        /// </returns>
        public double T(double lon)
        {
            SinCosd(lon, out var slam, out var clam);
            return InternalT(slam, clam, out _, out _, out _, false, true);
        }

        private (double W, double gX, double gY, double gZ) W(double slam, double clam)
        {
            var (Wres, gX, gY, gZ) = V(slam, clam);
            Wres += _frot * _Px / 2;
            gX += _frot * clam;
            gY += _frot * slam;
            return (Wres, gX, gY, gZ);
        }

        private (double V, double GX, double GY, double GZ) V(double slam, double clam)
        {
            if (!_caps.HasFlag(GravityFlags.Gravity))
            {
                return (double.NaN, double.NaN, double.NaN, double.NaN);
            }
            double
              Vres = _gravitational.Evaluate(slam, clam, out var GX, out var GY, out var GZ),
              f = _GMmodel / _amodel;
            Vres *= f;
            GX *= f;
            GY *= f;
            GZ *= f;
            return (Vres, GX, GY, GZ);
        }

        private double InternalT(double slam, double clam,
                     out double deltaX, out double deltaY, out double deltaZ,
                     bool gradp, bool correct)
        {
            deltaX = deltaY = deltaZ = double.NaN;
            if (gradp)
            {
                if (!_caps.HasFlag(GravityFlags.Disturbance))
                {
                    return double.NaN;
                }
            }
            else
            {
                if (!_caps.HasFlag(GravityFlags.DisturbingPotential))
                    return double.NaN;
            }

            if (_dzonal0 == 0)
                correct = false;
            var T = (gradp
                      ? _disturbing.Evaluate(slam, clam, out deltaX, out deltaY, out deltaZ)
                      : _disturbing.Evaluate(slam, clam));
            T = (T / _amodel - (correct ? _dzonal0 : 0) * _invR) * _GMmodel;
            if (gradp)
            {
                var f = _GMmodel / _amodel;
                deltaX *= f;
                deltaY *= f;
                deltaZ *= f;
                if (correct)
                {
                    var r3 = _GMmodel * _dzonal0 * _invR * _invR * _invR;
                    deltaX += _Px * clam * r3;
                    deltaY += _Px * slam * r3;
                    deltaZ += _Z * r3;
                }
            }
            return T;
        }
    }
}
