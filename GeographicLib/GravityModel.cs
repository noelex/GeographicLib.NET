using GeographicLib.SphericalHarmonics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using static GeographicLib.Macros;
using static GeographicLib.MathEx;
using static System.Math;

namespace GeographicLib
{
    /// <summary>
    /// Model of the earth's gravity field.
    /// </summary>
    /// <remarks>
    /// Evaluate the earth's gravity field according to a model.
    /// The supported models treat only the gravitational field exterior to the mass of the earth.
    /// When computing the field at points near (but above) the surface of the earth a small correction can be applied to account
    /// for the mass of the atmosphere above the point in question; see The effect of the mass of the atmosphere.
    /// Determining the height of the geoid above the ellipsoid entails correcting for the mass of the earth above the geoid.
    /// The egm96 and egm2008 include separate correction terms to account for this mass.
    /// <para>
    /// Definitions and terminology (from Heiskanen and Moritz, Sec 2-13):
    /// </para>
    /// <list type="bullet">
    /// <item><i>V</i> = gravitational potential;</item>
    /// <item>Φ = rotational potential;</item>
    /// <item><i>W</i> = <i>V</i> + Φ = <i>T</i> + <i>U</i> = total potential;</item>
    /// <item><i>V</i>0 = normal gravitation potential;</item>
    /// <item><i>U</i> = <i>V</i>0 + Φ = total normal potential;</item>
    /// <item><i>T</i> = <i>W</i> − <i>U</i> = <i>V</i> − <i>V</i>0 = anomalous or disturbing potential;</item>
    /// <item><b>g</b> = ∇W = <b>γ</b> + <b>δ</b>;</item>
    /// <item><b>f</b> = ∇Φ;</item>
    /// <item><b>Γ</b> = ∇V0;</item>
    /// <item><b>γ</b> = ∇U;</item>
    /// <item><b>δ</b> = ∇<i>T</i> = gravity disturbance vector = <b>g</b><i>P</i> − <b>γ</b><i>P</i>;</item>
    /// <item>δ<i>g</i> = gravity disturbance = <i>gP</i> − γ<i>P</i>;</item>
    /// <item>Δ<b>g</b> = gravity anomaly vector = <b>g</b><i>P</i> − <b>γ</b><i>P</i>; here the line <i>PQ</i> is perpendicular to ellipsoid and the potential at <i>P</i> equals the normal potential at <i>Q</i>;</item>
    /// <item>Δ<i>g</i> = gravity anomaly = <i>gP</i> − γ<i>Q</i>;</item>
    /// <item>(ξ, η) deflection of the vertical, the difference in directions of <b>g</b><i>P</i> and <b>γ</b><i>Q</i>, ξ = NS, η = EW.</item>
    /// <item><i>X</i>, <i>Y</i>, <i>Z</i>, geocentric coordinates;</item>
    /// <item><i>x</i>, <i>y</i>, <i>z</i>, local cartesian coordinates used to denote the east, north and up directions.</item>
    /// </list>
    /// See <a href="https://geographiclib.sourceforge.io/html/gravity.html">Gravity models</a> for details of how to install the gravity models and the data format.
    /// <para>References:</para>
    /// <list type="bullet">
    /// <item>
    /// W. A. Heiskanen and H. Moritz, Physical Geodesy (Freeman, San Francisco, 1967).
    /// </item>
    /// </list>
    /// </remarks>
    public class GravityModel : IEllipsoid
    {
        private const int idlength_ = 8;

        private readonly DateTime? _date;
        private readonly string _name, _dir, _description, _filename, _id;
        private readonly double _amodel, _GMmodel, _zeta0, _corrmult;
        private readonly int _nmx, _mmx;
        private readonly Normalization _norm;
        private readonly NormalGravity _earth;
        private readonly double _dzonal0;              // A left over contribution to _zonal.
        private readonly SphericalHarmonic _gravitational;
        private readonly SphericalHarmonic1 _disturbing;
        private readonly SphericalHarmonic _correction;

        static GravityModel()
        {
            string GetDefaultPath()
            {
                var path = Environment.GetEnvironmentVariable("GEOGRAPHICLIB_MAGNETIC_PATH");
                if (!string.IsNullOrEmpty(path))
                    return path;

                var datapath = Environment.GetEnvironmentVariable("GEOGRAPHICLIB_DATA");
                return Path.Combine(!string.IsNullOrEmpty(datapath) ? datapath : GEOGRAPHICLIB_DATA, "gravity");
            }

            DefaultGravityPath = GetDefaultPath();
        }

        /// <summary>
        /// Construct a gravity model.
        /// </summary>
        /// <param name="name">the name of the model.</param>
        /// <param name="path">directory for data file.</param>
        /// <param name="Nmax">if non-negative, truncate the degree of the model this value.</param>
        /// <param name="Mmax">if non-negative, truncate the order of the model this value.</param>
        /// <remarks>
        /// A filename is formed by appending ".egm" (World Gravity Model) to the <paramref name="name"/>.
        /// If <paramref name="path"/> is specified (and is non-empty), then the file is loaded from directory, <paramref name="path"/>.
        /// Otherwise the <paramref name="path"/> is given by <see cref="DefaultGravityPath"/>.
        /// <para>
        /// This file contains the metadata which specifies the properties of the model.
        /// The coefficients for the spherical harmonic sums are obtained from a file obtained by appending ".cof"
        /// to metadata file (so the filename ends in ".egm.cof").
        /// </para>
        /// <para>
        /// If <paramref name="Nmax"/> ≥ 0 and <paramref name="Mmax"/> &lt; 0, then <paramref name="Mmax"/> is set to <paramref name="Nmax"/>.
        /// After the model is loaded, the maximum degree and order of the model can be found by the <see cref="Degree"/> and <see cref="Order"/> methods.
        /// </para>
        /// </remarks>
        public GravityModel(string name, string path = "", int Nmax = -1, int Mmax = -1)
            : this()
        {
            bool truncate = Nmax >= 0 || Mmax >= 0;
            if (truncate)
            {
                if (Nmax >= 0 && Mmax < 0) Mmax = Nmax;
                if (Nmax < 0) Nmax = int.MaxValue;
                if (Mmax < 0) Mmax = int.MaxValue;
            }

            _name = name;
            _dir = string.IsNullOrEmpty(path) ? DefaultGravityPath : path;
            _filename = _dir + "/" + name + ".egm";

            using (var metadata = File.OpenRead(_filename))
            {
                ReadMetadata(_filename, metadata, ref _name, ref _description, ref _date, ref _amodel, ref _GMmodel,
                    ref _zeta0, ref _corrmult, ref _norm, ref _id, ref _earth);
            }

            string coeff = _filename + ".cof";
            using (var coeffstr = File.OpenRead(coeff))
            {
                ReadCoefficients(
                    coeff, coeffstr, truncate, Nmax, Mmax,
                    ref _gravitational,
                    ref _correction,
                    ref _disturbing,
                    ref _nmx,
                    ref _mmx,
                    ref _dzonal0);
            }
        }

        /// <summary>
        /// Construct a gravity model form the given <paramref name="metadataStream"/> and <paramref name="coefficientsStream"/>.
        /// </summary>
        /// <param name="metadataStream">A <see cref="Stream"/> which contains the metadata of the gravity model.</param>
        /// <param name="coefficientsStream">A <see cref="Stream"/> which contains the coefficients of the gravity model.</param>
        /// <param name="Nmax">if non-negative, truncate the degree of the model this value.</param>
        /// <param name="Mmax">if non-negative, truncate the order of the model this value.</param>
        /// <param name="leaveOpen">
        /// <see langword="true"/> to leave the streams open after the constructor returns; otherwise, <see langword="false"/>.
        /// </param>
        /// <remarks>
        /// <para>
        /// If <paramref name="Nmax"/> ≥ 0 and <paramref name="Mmax"/> &lt; 0, then <paramref name="Mmax"/> is set to <paramref name="Nmax"/>.
        /// After the model is loaded, the maximum degree and order of the model can be found by the <see cref="Degree"/> and <see cref="Order"/> methods.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"></exception>
        public GravityModel(Stream metadataStream, Stream coefficientsStream, int Nmax = -1, int Mmax = -1, bool leaveOpen = false)
            : this()
        {
            if (metadataStream == null)
            {
                throw new ArgumentNullException(nameof(metadataStream));
            }

            if (coefficientsStream == null)
            {
                throw new ArgumentNullException(nameof(coefficientsStream));
            }

            bool truncate = Nmax >= 0 || Mmax >= 0;
            if (truncate)
            {
                if (Nmax >= 0 && Mmax < 0) Mmax = Nmax;
                if (Nmax < 0) Nmax = int.MaxValue;
                if (Mmax < 0) Mmax = int.MaxValue;
            }

            _dir = null;
            _filename = null;
            _name = null;

            try
            {
                ReadMetadata(null, metadataStream, ref _name, ref _description, ref _date, ref _amodel, ref _GMmodel,
                    ref _zeta0, ref _corrmult, ref _norm, ref _id, ref _earth);

                ReadCoefficients(
                    null, coefficientsStream, truncate, Nmax, Mmax,
                    ref _gravitational,
                    ref _correction,
                    ref _disturbing,
                    ref _nmx,
                    ref _mmx,
                    ref _dzonal0);
            }
            finally
            {
                if (!leaveOpen)
                {
                    metadataStream.Dispose();
                    coefficientsStream.Dispose();
                }
            }
        }

        /// <summary>
        /// Construct a gravity model form the given <paramref name="metadataBytes"/> and <paramref name="coefficientsBytes"/>.
        /// </summary>
        /// <param name="metadataBytes">A <see cref="Byte"/> array which contains the metadata of the gravity model.</param>
        /// <param name="coefficientsBytes">A <see cref="Byte"/> array which contains the coefficients of the gravity model.</param>
        /// <param name="Nmax">if non-negative, truncate the degree of the model this value.</param>
        /// <param name="Mmax">if non-negative, truncate the order of the model this value.</param>
        /// <remarks>
        /// <para>
        /// If <paramref name="Nmax"/> ≥ 0 and <paramref name="Mmax"/> &lt; 0, then <paramref name="Mmax"/> is set to <paramref name="Nmax"/>.
        /// After the model is loaded, the maximum degree and order of the model can be found by the <see cref="Degree"/> and <see cref="Order"/> methods.
        /// </para>
        /// </remarks>
        public GravityModel(byte[] metadataBytes, byte[] coefficientsBytes, int Nmax = -1, int Mmax = -1)
            : this(new MemoryStream(metadataBytes), new MemoryStream(coefficientsBytes), Nmax, Mmax, leaveOpen: false)
        {

        }

        private GravityModel()
        {
            _description = "NONE";
            _amodel = double.NaN;
            _GMmodel = double.NaN;
            _zeta0 = 0;
            _corrmult = 1;
            _nmx = -1;
            _mmx = -1;
            _norm = Normalization.Full;
        }

        /// <summary>
        /// Evaluate the gravity at an arbitrary point above (or below) the ellipsoid.
        /// </summary>
        /// <param name="lat">the geographic latitude (degrees).</param>
        /// <param name="lon">the geographic longitude (degrees).</param>
        /// <param name="h">the height above the ellipsoid (meters).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>W</i>, the sum of the gravitational and centrifugal potentials (m^2 s^−2).</item>
        /// <item><i>gx</i>, the easterly component of the acceleration (m s^−2).</item>
        /// <item><i>gy</i>, the northerly component of the acceleration (m s^−2).</item>
        /// <item><i>gz</i>, the upward component of the acceleration (m s^−2); this is usually negative.</item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// The function includes the effects of the earth's rotation.
        /// </remarks>
        public (double W, double gx, double gy, double gz) Gravity(double lat, double lon, double h)
        {
            Span<double> M = stackalloc double[Geocentric.dim2_];
            var (X, Y, Z) = _earth.Earth.IntForward(lat, lon, h, M);
            var (Wres, gx, gy, gz) = W(X, Y, Z);
            Geocentric.Unrotate(M, gx, gy, gz, out gx, out gy, out gz);
            return (Wres, gx, gy, gz);
        }

        /// <summary>
        /// Evaluate the gravity disturbance vector at an arbitrary point above (or below) the ellipsoid.
        /// </summary>
        /// <param name="lat">the geographic latitude (degrees).</param>
        /// <param name="lon">the geographic longitude (degrees).</param>
        /// <param name="h">the height above the ellipsoid (meters).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>T</i>, the corresponding disturbing potential (m2 s−2).</item>
        /// <item><i>deltax</i>, the easterly component of the disturbance vector (m s^−2).</item>
        /// <item><i>deltay</i>, the northerly component of the disturbance vector (m s^−2).</item>
        /// <item><i>deltaz</i>, the upward component of the disturbance vector (m s^−2).</item>
        /// </list>
        /// </returns>
        public (double T, double deltax, double deltay, double deltaz) Disturbance(double lat, double lon, double h)
        {
            Span<double> M = stackalloc double[Geocentric.dim2_];
            var (X, Y, Z) = _earth.Earth.IntForward(lat, lon, h, M);
            var Tres = InternalT(X, Y, Z, out var deltax, out var deltay, out var deltaz, true, true);
            Geocentric.Unrotate(M, deltax, deltay, deltaz, out deltax, out deltay, out deltaz);
            return (Tres, deltax, deltay, deltaz);
        }

        /// <summary>
        /// Evaluate the geoid height.
        /// </summary>
        /// <param name="lat">the geographic latitude (degrees).</param>
        /// <param name="lon">the geographic longitude (degrees).</param>
        /// <returns><i>N</i>, the height of the geoid above the <see cref="ReferenceEllipsoid"/> (meters).</returns>
        /// <remarks>
        /// This calls <see cref="NormalGravity.U(double, double, double)"/> for
        /// <see cref="ReferenceEllipsoid"/>. Some approximations are made in computing the geoid height so that the results of
        /// the NGA codes are reproduced accurately. Details are given in
        /// <a href="https://geographiclib.sourceforge.io/html/gravity.html#gravitygeoid">Details of the geoid height and anomaly calculations</a>.
        /// </remarks>
        public double GeoidHeight(double lat, double lon)
        {
            var (X, Y, Z) = _earth.Earth.IntForward(lat, lon, 0);
            double
              gamma0 = _earth.SurfaceGravity(lat),
              T = InternalT(X, Y, Z, out _, out _, out _, false, false),
              invR = 1 / Hypot(Hypot(X, Y), Z),
              correction = _corrmult * _correction.Evaluate(invR * X, invR * Y, invR * Z);
            // _zeta0 has been included in _correction
            return T / gamma0 + correction;
        }

        /// <summary>
        /// Evaluate the components of the gravity anomaly vector using the spherical approximation.
        /// </summary>
        /// <param name="lat">the geographic latitude (degrees).</param>
        /// <param name="lon">the geographic longitude (degrees).</param>
        /// <param name="h">the height above the ellipsoid (meters).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>Dg01</i>, the gravity anomaly (m s^−2).</item>
        /// <item><i>xi</i>, the northerly component of the deflection of the vertical (degrees).</item>
        /// <item><i>eta</i>, the easterly component of the deflection of the vertical (degrees).</item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// The spherical approximation (see Heiskanen and Moritz, Sec 2-14) is used so that the results of
        /// the NGA codes are reproduced accurately. approximations used here.
        /// Details are given in
        /// <a href="https://geographiclib.sourceforge.io/html/gravity.html#gravitygeoid">Details of the geoid height and anomaly calculations</a>.
        /// </remarks>
        public (double Dg01, double xi, double eta) SphericalAnomaly(double lat, double lon, double h)
        {
            Span<double> M = stackalloc double[Geocentric.dim2_];
            var (X, Y, Z) = _earth.Earth.IntForward(lat, lon, h, M);
            double
              T = InternalT(X, Y, Z, out var deltax, out var deltay, out var deltaz, true, false),
              clam = M[3], slam = -M[0],
              P = Hypot(X, Y),
              R = Hypot(P, Z),
              // theta is geocentric latitude
              ctheta = R != 0 ? P / R : M[7],
              stheta = R != 0 ? Z / R : M[8];
            // Rotate cartesian into spherical coordinates
            Span<double> MC = stackalloc double[Geocentric.dim2_];
            Geocentric.Rotation(stheta, ctheta, slam, clam, MC);
            Geocentric.Unrotate(MC, deltax, deltay, deltaz, out deltax, out deltay, out deltaz);
            // H+M, Eq 2-151c
            var Dg01 = -deltaz - 2 * T / R;
            var (_, gammaX, gammaY, gammaZ) = _earth.U(X, Y, Z);
            var gamma = Hypot(Hypot(gammaX, gammaY), gammaZ);
            var xi = -(deltay / gamma) / Degree;
            var eta = -(deltax / gamma) / Degree;

            return (Dg01, xi, eta);
        }

        /// <summary>
        /// Evaluate the components of the acceleration due to gravity and the centrifugal acceleration in geocentric coordinates.
        /// </summary>
        /// <param name="X">the <i>X</i> component of geocentric coordinate of point (meters).</param>
        /// <param name="Y">the <i>Y</i> component of geocentric coordinate of point (meters).</param>
        /// <param name="Z">the <i>Z</i> component of geocentric coordinate of point (meters).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>W</i> = <i>V</i> + Φ, the sum of the gravitational and centrifugal potentials (m^2 s^−2).</item>
        /// <item><i>gX</i>, the <i>X</i> component of the acceleration (m s^−2).</item>
        /// <item><i>gY</i>, the <i>Y</i> component of the acceleration (m s^−2).</item>
        /// <item><i>gZ</i>, the <i>Z</i> component of the acceleration (m s^−2).</item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// This calls <see cref="NormalGravity.U(double, double, double)"/>
        /// for <see cref="ReferenceEllipsoid"/>.
        /// </remarks>
        public (double W, double gX, double gY, double gZ) W(double X, double Y, double Z)
        {
            var (Wres, gX, gY, gZ) = V(X, Y, Z);
            var (phi, fX, fY) = _earth.Phi(X, Y);
            Wres += phi;
            gX += fX;
            gY += fY;
            return (Wres, gX, gY, gZ);
        }

        /// <summary>
        /// Evaluate the components of the acceleration due to gravity in geocentric coordinates.
        /// </summary>
        /// <param name="X">the <i>X</i> component of geocentric coordinate of point (meters).</param>
        /// <param name="Y">the <i>Y</i> component of geocentric coordinate of point (meters).</param>
        /// <param name="Z">the <i>Z</i> component of geocentric coordinate of point (meters).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>V</i> = <i>W</i> - Φ, the gravitational potential (m^2 s^−2).</item>
        /// <item><i>GX</i>, the <i>X</i> component of the acceleration (m s^−2).</item>
        /// <item><i>GY</i>, the <i>Y</i> component of the acceleration (m s^−2).</item>
        /// <item><i>GZ</i>, the <i>Z</i> component of the acceleration (m s^−2).</item>
        /// </list>
        /// </returns>
        public (double V, double GX, double GY, double GZ) V(double X, double Y, double Z)
        {
            double
              Vres = _gravitational.Evaluate(X, Y, Z, out var GX, out var GY, out var GZ),
              f = _GMmodel / _amodel;
            Vres *= f;
            GX *= f;
            GY *= f;
            GZ *= f;
            return (Vres, GX, GY, GZ);
        }

        /// <summary>
        /// Evaluate the components of the gravity disturbance in geocentric coordinates.
        /// </summary>
        /// <param name="X">the <i>X</i> component of geocentric coordinate of point (meters).</param>
        /// <param name="Y">the <i>Y</i> component of geocentric coordinate of point (meters).</param>
        /// <param name="Z">the <i>Z</i> component of geocentric coordinate of point (meters).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>T</i> = <i>W</i> - <i>U</i>, the disturbing potential (also called the anomalous potential) (m^2 s^−2).</item>
        /// <item><i>deltaX</i>, the <i>X</i> component of the gravity disturbance (m s^−2).</item>
        /// <item><i>deltaY</i>, the <i>Y</i> component of the gravity disturbance (m s^−2).</item>
        /// <item><i>deltaZ</i>, the <i>Z</i> component of the gravity disturbance (m s^−2).</item>
        /// </list>
        /// </returns>
        public (double T, double deltaX, double deltaY, double deltaZ) Td(double X, double Y, double Z)
            => (InternalT(X, Y, Z, out var deltaX, out var deltaY, out var deltaZ, true, true), deltaX, deltaY, deltaZ);

        /// <summary>
        /// Evaluate disturbing potential in geocentric coordinates.
        /// </summary>
        /// <param name="X">the <i>X</i> component of geocentric coordinate of point (meters).</param>
        /// <param name="Y">the <i>Y</i> component of geocentric coordinate of point (meters).</param>
        /// <param name="Z">the <i>Z</i> component of geocentric coordinate of point (meters).</param>
        /// <returns>
        /// <i>T</i> = <i>W</i> - <i>U</i>, the disturbing potential (also called the anomalous potential) (m^2 s^−2).
        /// </returns>
        public double T(double X, double Y, double Z)
            => InternalT(X, Y, Z, out _, out _, out _, false, true);

        /// <summary>
        /// Evaluate the components of the acceleration due to normal gravity and the centrifugal acceleration in geocentric coordinates.
        /// </summary>
        /// <param name="X">the <i>X</i> component of geocentric coordinate of point (meters).</param>
        /// <param name="Y">the <i>Y</i> component of geocentric coordinate of point (meters).</param>
        /// <param name="Z">the <i>Z</i> component of geocentric coordinate of point (meters).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>U</i> = <i>V</i>0 + Φ, the sum of the normal gravitational and centrifugal potentials (m^2 s^−2).</item>
        /// <item><i>gammaX</i>, the <i>X</i> component of the normal acceleration (m s^−2).</item>
        /// <item><i>gammaY</i>, the <i>Y</i> component of the normal acceleration (m s^−2).</item>
        /// <item><i>gammaZ</i>, the <i>Z</i> component of the normal acceleration (m s^−2).</item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// This calls <see cref="NormalGravity.U(double, double, double)"/>
        /// for <see cref="ReferenceEllipsoid"/>.
        /// </remarks>
        public (double U, double gammaX, double gammaY, double gammaZ) U(double X, double Y, double Z)
            => _earth.U(X, Y, Z);

        /// <summary>
        /// Evaluate the centrifugal acceleration in geocentric coordinates.
        /// </summary>
        /// <param name="X">the <i>X</i> component of geocentric coordinate of point (meters).</param>
        /// <param name="Y">the <i>Y</i> component of geocentric coordinate of point (meters).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item>Φ, the centrifugal potential (m^2 s^−2).</item>
        /// <item><i>fX</i>, the <i>X</i> component of the centrifugal acceleration (m s^−2).</item>
        /// <item><i>fY</i>, the <i>Y</i> component of the centrifugal acceleration (m s^−2).</item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// This calls <see cref="NormalGravity.Phi(double, double)"/> for <see cref="ReferenceEllipsoid"/>.
        /// </remarks>
        public (double phi, double fX, double fY) Phi(double X, double Y)
            => _earth.Phi(X, Y);

        /// <summary>
        /// Create a <see cref="GravityCircle"/> object to allow the gravity field at many points with constant
        /// <paramref name="lat"/> and <paramref name="h"/> and varying <i>lon</i> to be computed efficiently.
        /// </summary>
        /// <param name="lat">latitude of the point (degrees).</param>
        /// <param name="h">the height of the point above the ellipsoid (meters).</param>
        /// <param name="caps">bitor'ed combination of <see cref="GravityFlags"/> values specifying the capabilities of the resulting
        /// <see cref="GravityCircle"/> object.</param>
        /// <returns>a <see cref="GravityCircle"/> object whose member functions computes the gravitational field at a particular values
        /// of <i>lon</i>.</returns>
        /// <remarks>
        /// The <see cref="GravityFlags"/> values are
        /// <list type="bullet">
        /// <item><paramref name="caps"/> |= <see cref="GravityFlags.Gravity"/></item>
        /// <item><paramref name="caps"/> |= <see cref="GravityFlags.Disturbance"/></item>
        /// <item><paramref name="caps"/> |= <see cref="GravityFlags.DisturbingPotential"/></item>
        /// <item><paramref name="caps"/> |= <see cref="GravityFlags.SphericalAnomaly"/></item>
        /// <item><paramref name="caps"/> |= <see cref="GravityFlags.GeoidHeight"/></item>
        /// </list>
        /// The default value of <paramref name="caps"/> is <see cref="GravityFlags.All"/> which turns on all the capabilities.
        /// If an unsupported function is invoked, it will return <see cref="double.NaN"/>.
        /// Note that <see cref="GravityFlags.GeoidHeight"/> will only be honored if <paramref name="h"/> = 0.
        /// <para>
        /// If the field at several points on a circle of latitude need to be calculated then creating a <see cref="GravityCircle"/>
        /// object and using its member functions will be substantially faster, especially for high-degree models.
        /// See <a href="https://geographiclib.sourceforge.io/html/gravity.html#gravityparallel">Geoid heights on a multi-processor
        /// system</a> for an example of using <see cref="GravityCircle"/> (together with OpenMP) to speed up the computation of geoid heights.
        /// </para>
        /// </remarks>
        public GravityCircle Circle(double lat, double h, GravityFlags caps = GravityFlags.All)
        {
            if (h != 0)
                // Disallow invoking GeoidHeight unless h is zero.
                caps &= ~(GravityFlags)(GravityCapability.Gamma0 | GravityCapability.C);
            Span<double> M = stackalloc double[Geocentric.dim2_];
            var (X, Y, Z) = _earth.Earth.IntForward(lat, 0, h, M);
            // Y = 0, cphi = M[7], sphi = M[8];
            double
              invR = 1 / Hypot(X, Z),
              gamma0 = ((int)caps & (int)GravityCapability.Gamma0) != 0 ? _earth.SurfaceGravity(lat)
                        : double.NaN,
              fx, fy, gamma;
            if (((int)caps & (int)GravityCapability.Gamma) != 0)
            {
                double fz;
                (_, fx, fy, fz) = _earth.U(X, Y, Z); // fy = 0
                gamma = Hypot(fx, fz);
            }
            else
                gamma = double.NaN;
            (_, fx, fy) = _earth.Phi(X, Y);
            return new GravityCircle(caps,
                                 _earth._a, _earth._f, lat, h, Z, X, M[7], M[8],
                                 _amodel, _GMmodel, _dzonal0, _corrmult,
                                 gamma0, gamma, fx,
                                 ((int)caps & (int)GravityCapability.G) != 0 ?
                                 _gravitational.Circle(X, Z, true) :
                                 new CircularEngine(),
                                 // N.B. If CAP_DELTA is set then CAP_T should be too.
                                 ((int)caps & (int)GravityCapability.T) != 0 ?
                                 _disturbing.Circle(-1, X, Z, ((int)caps & (int)GravityCapability.Delta) != 0) :
                                 new CircularEngine(),
                                 ((int)caps & (int)GravityCapability.C) != 0 ?
                                 _correction.Circle(invR * X, invR * Z, false) :
                                 new CircularEngine());
        }

        /// <summary>
        /// Gets a value representing the <see cref="NormalGravity"/> object for the reference ellipsoid.
        /// </summary>
        public NormalGravity ReferenceEllipsoid => _earth;

        /// <summary>
        /// Gets a value representing the description of the gravity model, if available, in the data file; if absent, return "NONE".
        /// </summary>
        public string Description => _description;

        /// <summary>
        /// Gets a value representing date of the model; if absent, return <see langword="null"/>.
        /// </summary>
        public DateTime? DateTime => _date;

        /// <summary>
        /// Gets a value representing full file name used to load the gravity model.
        /// </summary>
        /// <remarks>
        /// This property returns <see langword="null"/> if the <see cref="GravityModel"/> object
        /// is constructed from <see cref="Stream"/> or <see cref="Byte"/> array.
        /// </remarks>
        public string GravityFile => _filename;

        /// <summary>
        /// Gets a value representing "name" used to load the gravity model
        /// (from the first argument of the constructor, but this may be overridden by the model file).
        /// </summary>
        /// <remarks>
        /// This property returns <see langword="null"/> if the <see cref="GravityModel"/> object
        /// is constructed from <see cref="Stream"/> or <see cref="Byte"/> array, and the model file doesn't
        /// define <c>Name</c> attribute.
        /// </remarks>
        public string GravityModelName => _name;

        /// <summary>
        /// Gets a value representing directory used to load the gravity model.
        /// </summary>
        ///  <remarks>
        /// This property returns <see langword="null"/> if the <see cref="GravityModel"/> object
        /// is constructed from <see cref="Stream"/> or <see cref="Byte"/> array.
        /// </remarks>
        public string GravityModelDirectory => _dir;

        /// <inheritdoc/>
        public double EquatorialRadius => _earth.EquatorialRadius;

        /// <inheritdoc/>
        public double Flattening => _earth.Flattening;

        /// <summary>
        /// Gets a value representing <i>GM</i>, the mass constant of the model (m^3 s^−2);
        /// this is the product of <i>G</i> the gravitational constant and <i>M</i> the mass of the earth
        /// (usually including the mass of the earth's atmosphere).
        /// </summary>
        public double MassConstant => _GMmodel;

        /// <summary>
        /// Gets a value representing <i>GM</i>, the mass constant of the <see cref="ReferenceEllipsoid"/> (m^3 s^−2).
        /// </summary>
        public double ReferenceMassConstant => _earth.MassConstant;

        /// <summary>
        /// Gets a value representing ω, the angular velocity of the model and the <see cref="ReferenceEllipsoid"/> (rad s^−1).
        /// </summary>
        public double AngularVelocity => _earth.AngularVelocity;

        /// <summary>
        /// Gets a value representing <i>Nmax</i>, the maximum degree of the components of the model.
        /// </summary>
        public int Degree => _nmx;

        /// <summary>
        /// Gets a value representing <i>Mmax</i>, the maximum order of the components of the model.
        /// </summary>
        public int Order => _mmx;

        /// <summary>
        /// Gets a value representing the default path for gravity model data files.
        /// </summary>
        /// <remarks>
        /// This is the value of the environment variable GEOGRAPHICLIB_GRAVITY_PATH, if set;
        /// otherwise, it is $GEOGRAPHICLIB_DATA/gravity if the environment variable GEOGRAPHICLIB_DATA is set;
        /// otherwise, it is a compile-time default (/usr/local/share/GeographicLib/gravity on non-Windows systems
        /// and C:/ProgramData/GeographicLib/gravity on Windows systems).
        /// </remarks>
        public static string DefaultGravityPath { get; }

        /// <summary>
        /// Gets a value representing the default name for the gravity model.
        /// </summary>
        /// <remarks>
        /// This is the value of the environment variable GEOGRAPHICLIB_GRAVITY_NAME, if set;
        /// otherwise, it is "egm96". The <see cref="GravityModel"/> class does not use this function;
        /// it is just provided as a convenience for a calling program when constructing a <see cref="GravityModel"/> object.
        /// </remarks>
        public static string DefaultGravityName { get; } = Environment.GetEnvironmentVariable("GEOGRAPHICLIB_GRAVITY_NAME") ?? "egm96";

        private void ReadMetadata(string name, Stream stream, ref string _name, ref string _description,
            ref DateTime? _date, ref double _amodel, ref double _GMmodel, ref double _zeta0, ref double _corrmult,
            ref Normalization _norm, ref string _id, ref NormalGravity _earth)
        {
            name = name == null ? "metadata stream" : name;
            using (var metastr = new StreamReader(stream, Encoding.UTF8, true, bufferSize: 1024, leaveOpen: true))
            {
                var line = metastr.ReadLine();
                if (!(line.Length >= 6 && line.StartsWith("EGMF-")))
                    throw new GeographicException(name + " does not contain EGMF-n signature");

                var parts = line.TrimEnd().Split(new[] { "EGMF-" }, StringSplitOptions.RemoveEmptyEntries);
                if (parts.Length == 0 || parts[0] != "1")
                    throw new GeographicException("Unknown version in " + name + ": " + parts[0]);
                double a = double.NaN, GM = a, omega = a, f = a, J2 = a;
                while ((line = metastr.ReadLine()) != null)
                {
                    if (!Utility.ParseLine(line.AsSpan(), out var key, out var val))
                    {
                        continue;
                    }
                    // Process key words
                    if (key == "Name")
                        _name = val;
                    else if (key == "Description")
                        _description = val;
                    else if (key == "ReleaseDate")
                        _date = val.ParseDateTime();
                    else if (key == "ModelRadius")
                        _amodel = val.ParseDoubleInfNan();
                    else if (key == "ModelMass")
                        _GMmodel = val.ParseDoubleInfNan();
                    else if (key == "AngularVelocity")
                        omega = val.ParseDoubleInfNan();
                    else if (key == "ReferenceRadius")
                        a = val.ParseDoubleInfNan();
                    else if (key == "ReferenceMass")
                        GM = val.ParseDoubleInfNan();
                    else if (key == "Flattening")
                        f = val.ParseFract();
                    else if (key == "DynamicalFormFactor")
                        J2 = val.ParseFract();
                    else if (key == "HeightOffset")
                        _zeta0 = val.ParseFract();
                    else if (key == "CorrectionMultiplier")
                        _corrmult = val.ParseFract();
                    else if (key == "Normalization")
                    {
                        if (val == "FULL" || val == "Full" || val == "full")
                            _norm = Normalization.Full;
                        else if (val == "SCHMIDT" || val == "Schmidt" || val == "schmidt")
                            _norm = Normalization.Schmidt;
                        else
                            throw new GeographicException("Unknown normalization " + val);
                    }
                    else if (key == "ByteOrder")
                    {
                        if (val == "Big" || val == "big")
                            throw new GeographicException("Only little-endian ordering is supported");
                        else if (!(val == "Little" || val == "little"))
                            throw new GeographicException("Unknown byte ordering " + val);
                    }
                    else if (key == "ID")
                        _id = val;
                    // else unrecognized keywords are skipped
                }
                // Check values
                if (!(IsFinite(_amodel) && _amodel > 0))
                    throw new GeographicException("Model radius must be positive");
                if (!(IsFinite(_GMmodel) && _GMmodel > 0))
                    throw new GeographicException("Model mass constant must be positive");
                if (!(IsFinite(_corrmult) && _corrmult > 0))
                    throw new GeographicException("Correction multiplier must be positive");
                if (!(IsFinite(_zeta0)))
                    throw new GeographicException("Height offset must be finite");
                if (_id.Length != idlength_)
                    throw new GeographicException("Invalid ID");
                if (IsFinite(f) && IsFinite(J2))
                    throw new GeographicException("Cannot specify both f and J2");
                _earth = new NormalGravity(a, GM, omega,
                                       IsFinite(f) ? f : J2, IsFinite(f));
            }
        }

        private void ReadCoefficients(
            string name,
            Stream coeffstr,
            bool truncate, int Nmax, int Mmax,
            ref SphericalHarmonic _gravitational,
            ref SphericalHarmonic _correction,
            ref SphericalHarmonic1 _disturbing,
            ref int _nmx,
            ref int _mmx,
            ref double _dzonal0)
        {
            name = name == null ? "coefficients stream" : name;

            double[] _Cx, _Sx;
            Span<byte> id = stackalloc byte[idlength_];
            if (coeffstr.Read(id) != idlength_)
                throw new GeographicException("No header in " + name);
            if (MemoryMarshal.Cast<char, byte>(_id.AsSpan()).SequenceEqual(id))
                throw new GeographicException($"ID mismatch: {_id} vs {Encoding.ASCII.GetString(id.ToArray())}");
            int N = 0, M = 0;
            if (truncate) { N = Nmax; M = Mmax; }

            var scoeff = SphericalEngine.Coeff.FromStream(coeffstr, ref N, ref M, truncate);
            _Cx = new double[scoeff.Cnm.Length];
            _Sx = new double[scoeff.Snm.Length];
            scoeff.Cnm.CopyTo(_Cx);
            scoeff.Snm.CopyTo(_Sx);

            if (!(N >= 0 && M >= 0))
                throw new GeographicException("Degree and order must be at least 0");
            if (_Cx[0] != 0)
                throw new GeographicException("The degree 0 term should be zero");

            _Cx[0] = 1;               // Include the 1/r term in the sum
            _gravitational = new SphericalHarmonic(_Cx, _Sx, N, N, M, _amodel, _norm);
            if (truncate) { N = Nmax; M = Mmax; }


            scoeff = SphericalEngine.Coeff.FromStream(coeffstr, ref N, ref M, truncate);
            double[] _CC, _CS;
            if (N < 0)
            {
                N = M = 0;
                _CC = new[] { 0d };
            }
            else
            {
                _CC = new double[scoeff.Cnm.Length];
                scoeff.Cnm.CopyTo(_CC);
            }

            _CS = new double[scoeff.Snm.Length];
            scoeff.Snm.CopyTo(_CS);

            _CC[0] += _zeta0 / _corrmult;
            _correction = new SphericalHarmonic(_CC, _CS, N, N, M, 1, _norm);
            var pos = (int)coeffstr.Position;
            coeffstr.Seek(0, SeekOrigin.End);
            if (pos != coeffstr.Position)
                throw new GeographicException("Extra data in " + name);

            int nmx = _gravitational.Coefficients.Nmx;
            _nmx = Max(nmx, _correction.Coefficients.Nmx);
            _mmx = Max(_gravitational.Coefficients.Mmx,
                       _correction.Coefficients.Mmx);
            // Adjust the normalization of the normal potential to match the model.
            var mult = _earth._GM / _GMmodel;
            var amult = Sq(_earth._a / _amodel);
            // The 0th term in _zonal should be is 1 + _dzonal0.  Instead set it to 1
            // to give exact cancellation with the (0,0) term in the model and account
            // for _dzonal0 separately.
            var _zonal = new List<double> { 1 };
            _dzonal0 = (_earth.MassConstant - _GMmodel) / _GMmodel;
            for (int n = 2; n <= nmx; n += 2)
            {
                // Only include as many normal zonal terms as matter.  Figuring the limit
                // in this way works because the coefficients of the normal potential
                // (which is smooth) decay much more rapidly that the corresponding
                // coefficient of the model potential (which is bumpy).  Typically this
                // goes out to n = 18.
                mult *= amult;
                double
                  r = _Cx[n],                                        // the model term
                  s = -mult * _earth.Jn(n) / Sqrt(2 * n + 1), // the normal term
                  t = r - s;                                         // the difference
                if (t == r)               // the normal term is negligible
                    break;
                _zonal.Add(0);      // index = n - 1; the odd terms are 0
                _zonal.Add(s);
            }
            int nmx1 = _zonal.Count - 1;
            var za = _zonal.ToArray();
            _disturbing = new SphericalHarmonic1(_Cx, _Sx,
                                             _gravitational.Coefficients.N,
                                             nmx, _gravitational.Coefficients.Mmx,
                                             za,
                                             za, // This is not accessed!
                                             nmx1, nmx1, 0,
                                             _amodel,
                                             _norm);
        }

        private double InternalT(double X, double Y, double Z,
                             out double deltaX, out double deltaY, out double deltaZ,
                             bool gradp, bool correct)
        {
            deltaX = deltaY = deltaZ = double.NaN;
            // If correct, then produce the correct T = W - U.  Otherwise, neglect the
            // n = 0 term (which is proportial to the difference in the model and
            // reference values of GM).
            if (_dzonal0 == 0)
                // No need to do the correction
                correct = false;
            double T, invR = correct ? 1 / Hypot(Hypot(X, Y), Z) : 1;
            if (gradp)
            {
                T = _disturbing.Evaluate(-1, X, Y, Z, out deltaX, out deltaY, out deltaZ);
                var f = _GMmodel / _amodel;
                deltaX *= f;
                deltaY *= f;
                deltaZ *= f;
                if (correct)
                {
                    invR = _GMmodel * _dzonal0 * invR * invR * invR;
                    deltaX += X * invR;
                    deltaY += Y * invR;
                    deltaZ += Z * invR;
                }
            }
            else
                T = _disturbing.Evaluate(-1, X, Y, Z);
            T = (T / _amodel - (correct ? _dzonal0 : 0) * invR) * _GMmodel;
            return T;
        }
    }
}
