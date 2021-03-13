using GeographicLib.SphericalHarmonics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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
        /// <param name="gx">the easterly component of the acceleration (m s^−2).</param>
        /// <param name="gy">the northerly component of the acceleration (m s^−2).</param>
        /// <param name="gz">the upward component of the acceleration (m s^−2); this is usually negative.</param>
        /// <returns>
        /// <i>W</i>, the sum of the gravitational and centrifugal potentials (m^2 s^−2).
        /// </returns>
        /// <remarks>The function includes the effects of the earth's rotation.</remarks>
        public double Gravity(double lon, out double gx, out double gy, out double gz)
        {
            Span<double> M=stackalloc double[Geocentric.dim2_];
            SinCosd(lon, out var slam, out var clam);
            var Wres = W(slam, clam, out gx, out gy, out gz);
            Geocentric.Rotation(_sphi, _cphi, slam, clam, M);
            Geocentric.Unrotate(M, gx, gy, gz,out gx, out gy, out gz);
            return Wres;
        }

        /// <summary>
        /// Evaluate the gravity disturbance vector.
        /// </summary>
        /// <param name="lon">the geographic longitude (degrees).</param>
        /// <param name="deltax">the easterly component of the disturbance vector (m s^−2).</param>
        /// <param name="deltay">the northerly component of the disturbance vector (m s^−2).</param>
        /// <param name="deltaz">the upward component of the disturbance vector (m s^−2).</param>
        /// <returns><i>T</i>, the corresponding disturbing potential (m2 s−2).</returns>
        public double Disturbance(double lon, out double deltax, out double deltay, out double deltaz)
        {
            Span<double> M = stackalloc double[Geocentric.dim2_];
            SinCosd(lon, out var slam, out var clam);
            var Tres = InternalT(slam, clam, out deltax, out deltay, out deltaz, true, true);
            Geocentric.Rotation(_sphi, _cphi, slam, clam, M);
            Geocentric.Unrotate(M, deltax, deltay, deltaz, out deltax, out deltay, out deltaz);
            return Tres;
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
            if (_caps .HasFlag( GravityFlags.GeoidHeight))
                return double.NaN;

            SinCosd(lon, out var slam, out var clam);
            var T = InternalT(slam, clam, out _, out _, out _, false, false);
            var correction = _corrmult * _correction.Evaluate(slam, clam);
            return T / _gamma0 + correction;
        }

        /// <summary>
        /// Evaluate the components of the gravity anomaly vector using the spherical approximation.
        /// </summary>
        /// <param name="lon">the geographic longitude (degrees).</param>
        /// <param name="Dg01">the gravity anomaly (m s^−2).</param>
        /// <param name="xi">the northerly component of the deflection of the vertical (degrees).</param>
        /// <param name="eta">the easterly component of the deflection of the vertical (degrees).</param>
        /// <remarks>
        /// The spherical approximation (see Heiskanen and Moritz, Sec 2-14) is used so that the results of the NGA codes are
        /// reproduced accurately. approximations used here.
        /// Details are given in
        /// <a href="https://geographiclib.sourceforge.io/html/gravity.html#gravitygeoid">Details of the geoid height and anomaly calculations</a>.
        /// </remarks>
        public void SphericalAnomaly(double lon, out double Dg01, out double xi, out double eta)
        {
            if (_caps.HasFlag(GravityFlags.SphericalAnomaly))
            {
                Dg01 = xi = eta = double.NaN;
                return;
            }

            SinCosd(lon, out var slam, out var clam);
            double
              deltax, deltay, deltaz,
              T = InternalT(slam, clam, out deltax, out deltay, out deltaz, true, false);
            // Rotate cartesian into spherical coordinates
            Span<double> MC=stackalloc double[Geocentric.dim2_];
            Geocentric.Rotation(_spsi, _cpsi, slam, clam, MC);
            Geocentric.Unrotate(MC, deltax, deltay, deltaz, out deltax, out deltay, out deltaz);
            // H+M, Eq 2-151c
            Dg01 = -deltaz - 2 * T * _invR;
            xi = -(deltay / _gamma) / Degree;
            eta = -(deltax / _gamma) / Degree;
        }

        /// <summary>
        /// Evaluate the components of the acceleration due to gravity and the centrifugal acceleration in geocentric coordinates.
        /// </summary>
        /// <param name="lon">the geographic longitude (degrees).</param>
        /// <param name="gX">the <i>X</i> component of the acceleration (m s^−2).</param>
        /// <param name="gY">the <i>Y</i> component of the acceleration (m s^−2).</param>
        /// <param name="gZ">the <i>Z</i> component of the acceleration (m s^−2).</param>
        /// <returns><i>W</i> = <i>V</i> + Φ, the sum of the gravitational and centrifugal potentials (m^2 s^−2).</returns>
        public double W(double lon, out double gX, out double gY, out double gZ)
        {
            SinCosd(lon, out var slam, out var clam);
            return W(slam, clam, out gX, out gY, out gZ);
        }

        /// <summary>
        /// Evaluate the components of the acceleration due to gravity in geocentric coordinates.
        /// </summary>
        /// <param name="lon">the geographic longitude (degrees).</param>
        /// <param name="GX">the <i>X</i> component of the acceleration (m s^−2).</param>
        /// <param name="GY">the <i>Y</i> component of the acceleration (m s^−2).</param>
        /// <param name="GZ">the <i>Z</i> component of the acceleration (m s^−2).</param>
        /// <returns><i>V</i> = <i>W</i> - Φ, the gravitational potential (m^2 s^−2).</returns>
        public double V(double lon, out double GX, out double GY, out double GZ)
        {
            SinCosd(lon, out var slam, out var clam);
            return V(slam, clam, out GX, out GY, out GZ);
        }

        /// <summary>
        /// Evaluate the components of the gravity disturbance in geocentric coordinates.
        /// </summary>
        /// <param name="lon">the geographic longitude (degrees).</param>
        /// <param name="deltaX">the <i>X</i> component of the gravity disturbance (m s^−2).</param>
        /// <param name="deltaY">the <i>Y</i> component of the gravity disturbance (m s^−2).</param>
        /// <param name="deltaZ">the <i>Z</i> component of the gravity disturbance (m s^−2).</param>
        /// <returns>
        /// <i>T</i> = <i>W</i> - <i>U</i> the disturbing potential (also called the anomalous potential) (m^2 s^−2).
        /// </returns>
        public double T(double lon, out double deltaX, out double deltaY, out double deltaZ)
        {
            SinCosd(lon, out var slam, out var clam);
            return InternalT(slam, clam, out deltaX, out deltaY, out deltaZ, true, true);
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

        private double W(double slam, double clam,
                     out double gX, out double gY, out double gZ)
        {
            var Wres = V(slam, clam, out gX, out gY, out gZ) + _frot * _Px / 2;
            gX += _frot * clam;
            gY += _frot * slam;
            return Wres;
        }

        private double V(double slam, double clam,
                     out double GX, out double GY, out double GZ)
        {
            if (_caps.HasFlag(GravityFlags.Gravity) )
            {
                GX = GY = GZ = double.NaN;
                return double.NaN;
            }
            double
              Vres = _gravitational.Evaluate(slam, clam, out GX, out GY, out GZ),
              f = _GMmodel / _amodel;
            Vres *= f;
            GX *= f;
            GY *= f;
            GZ *= f;
            return Vres;
        }

        private double InternalT(double slam, double clam,
                     out double deltaX, out double deltaY, out double deltaZ,
                     bool gradp, bool correct)
        {
            deltaX = deltaY = deltaZ = double.NaN;
            if (gradp)
            {
                if (_caps.HasFlag(GravityFlags.Disturbance))
                {
                    return double.NaN;
                }
            }
            else
            {
                if (_caps.HasFlag(GravityFlags.DisturbingPotential))
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
