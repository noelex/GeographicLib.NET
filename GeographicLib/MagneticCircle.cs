using System;
using System.Collections.Generic;
using System.Text;
using GeographicLib.SphericalHarmonics;

using static GeographicLib.MathEx;

namespace GeographicLib
{
    /// <summary>
    /// Geomagnetic field on a circle of latitude.
    /// </summary>
    /// <remarks>
    /// Evaluate the earth's magnetic field on a circle of constant height and latitude.
    /// This uses a <see cref="CircularEngine"/> to pre-evaluate the inner sum of the spherical harmonic sum,
    /// allowing the values of the field at several different longitudes to be evaluated rapidly.
    /// <para>
    /// Use <see cref="MagneticModel.Circle"/>  to create a <see cref="MagneticCircle"/> object. (The constructor for this class is private.)
    /// </para>
    /// </remarks>
    public class MagneticCircle : IEllipsoid
    {
        private readonly double _a, _f, _lat, _h, _t, _cphi, _sphi, _t1, _dt0;
        private readonly bool _interpolate, _constterm;
        private readonly CircularEngine _circ0, _circ1, _circ2;

        internal MagneticCircle(double a, double f, double lat, double h, double t,
                                double cphi, double sphi, double t1, double dt0,
                                bool interpolate,
                                CircularEngine circ0, CircularEngine circ1, CircularEngine circ2 = null)
        {
            if (circ0 is null || circ1 is null)
            {
                throw new ArgumentNullException("circ0 and circ1 cannot be null.");
            }

            _a = a;
            _f = f;
            _lat = LatFix(lat);
            _h = h;
            _t = t;
            _cphi = cphi;
            _sphi = sphi;
            _t1 = t1;
            _dt0 = dt0;
            _interpolate = interpolate;
            _constterm = circ2 != null;
            _circ0 = circ0;
            _circ1 = circ1;
            _circ2 = circ2;
        }

        internal MagneticCircle(IEllipsoid ellipsoid, double lat, double h, double t,
                        double cphi, double sphi, double t1, double dt0,
                        bool interpolate,
                        CircularEngine circ0, CircularEngine circ1, CircularEngine circ2 = null)
            : this(ellipsoid.EquatorialRadius, ellipsoid.Flattening, lat, h, t, cphi, sphi, t1, dt0, interpolate, circ0, circ1, circ2) { }

        private (double Bx, double By, double Bz) Field(double lon, bool diffp,
                           out double Bxt, out double Byt, out double Bzt)
        {
            Bxt = Byt = Bzt = default;

            SinCosd(lon, out var slam, out var clam);
            Span<double> M = stackalloc double[Geocentric.dim2_];
            Geocentric.Rotation(_sphi, _cphi, slam, clam, M);

            var (BX, BY, BZ) = FieldGeocentric(slam, clam, out var BXt, out var BYt, out var BZt);
            if (diffp)
                Geocentric.Unrotate(M, BXt, BYt, BZt, out Bxt, out Byt, out Bzt);
            Geocentric.Unrotate(M, BX, BY, BZ, out var Bx, out var By, out var Bz);

            return (Bx, By, Bz);
        }

        private (double BX, double BY, double BZ) FieldGeocentric(double slam, double clam,
                         out double BXt, out double BYt, out double BZt)
        {
            double BXc = 0, BYc = 0, BZc = 0;
            _circ0.Evaluate(slam, clam, out var BX, out var BY, out var BZ);
            _circ1.Evaluate(slam, clam, out BXt, out BYt, out BZt);
            if (_constterm)
                _circ2.Evaluate(slam, clam, out BXc, out BYc, out BZc);
            if (_interpolate)
            {
                BXt = (BXt - BX) / _dt0;
                BYt = (BYt - BY) / _dt0;
                BZt = (BZt - BZ) / _dt0;
            }
            BX += _t1 * BXt + BXc;
            BY += _t1 * BYt + BYc;
            BZ += _t1 * BZt + BZc;

            BXt *= -_a;
            BYt *= -_a;
            BZt *= -_a;

            BX *= -_a;
            BY *= -_a;
            BZ *= -_a;

            return (BX, BY, BZ);
        }

        /// <summary>
        /// Evaluate the components of the geomagnetic field at a particular longitude.
        /// </summary>
        /// <param name="lon">longitude of the point (degrees).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>Bx</i>, the easterly component of the magnetic field (nanotesla).</item>
        /// <item><i>By</i>, the northerly component of the magnetic field (nanotesla).</item>
        /// <item><i>Bz</i>, the vertical (up) component of the magnetic field (nanotesla).</item>
        /// </list>
        /// </returns>
        public (double Bx, double By, double Bz) Evaluate(double lon)
            => Field(lon, false, out _, out _, out _);

        /// <summary>
        /// Evaluate the components of the geomagnetic field and their time derivatives at a particular longitude.
        /// </summary>
        /// <param name="lon">longitude of the point (degrees).</param>
        /// <param name="Bxt">the rate of change of <i>Bx</i> (nT/yr).</param>
        /// <param name="Byt">the rate of change of <i>By</i> (nT/yr).</param>
        /// <param name="Bzt">the rate of change of <i>Bz</i> (nT/yr).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>Bx</i>, the easterly component of the magnetic field (nanotesla).</item>
        /// <item><i>By</i>, the northerly component of the magnetic field (nanotesla).</item>
        /// <item><i>Bz</i>, the vertical (up) component of the magnetic field (nanotesla).</item>
        /// </list>
        /// </returns>
        public (double Bx, double By, double Bz) Evaluate(double lon,
                             out double Bxt, out double Byt, out double Bzt)
            => Field(lon, false, out Bxt, out Byt, out Bzt);

        /// <summary>
        /// Evaluate the components of the geomagnetic field and their time derivatives at a particular longitude.
        /// </summary>
        /// <param name="lon">longitude of the point (degrees).</param>
        /// <param name="BXt">the rate of change of <i>BX</i> (nT/yr).</param>
        /// <param name="BYt">the rate of change of <i>BY</i> (nT/yr).</param>
        /// <param name="BZt">the rate of change of <i>BZ</i> (nT/yr).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>BX</i>, the <i>X</i> component of the magnetic field (nanotesla).</item>
        /// <item><i>BY</i>, the <i>Y</i> component of the magnetic field (nanotesla).</item>
        /// <item><i>BZ</i>, the <i>Z</i> component of the magnetic field (nanotesla).</item>
        /// </list>
        /// </returns>
        public (double BX, double BY, double BZ) FieldGeocentric(double lon, out double BXt, out double BYt, out double BZt)
        {
            SinCosd(lon, out var slam, out var clam);
            return FieldGeocentric(slam, clam, out BXt, out BYt, out BZt);
        }

        /// <summary>
        /// Gets a value representing the equatorial radius (<i>a</i>) of the ellipsoid.
        /// </summary>
        public double EquatorialRadius => _a;

        /// <summary>
        /// Gets a value representing the flattening (<i>f</i>) of the ellipsoid.
        /// </summary>
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
        /// Gets a value representing the time (fractional years).
        /// </summary>
        public double Time => _t;
    }
}
