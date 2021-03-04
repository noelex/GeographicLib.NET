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

        private void Field(double lon, bool diffp,
                           out double Bx, out double By, out double Bz,
                           out double Bxt, out double Byt, out double Bzt)
        {
            SinCosd(lon, out var slam, out var clam);
            Span<double> M = stackalloc double[Geocentric.dim2_];
            Geocentric.Rotation(_sphi, _cphi, slam, clam, M);

            double BXc = 0, BYc = 0, BZc = 0;
            _circ0.Evaluate(slam, clam, out var BX0, out var BY0, out var BZ0);
            _circ1.Evaluate(slam, clam, out var BX1, out var BY1, out var BZ1);

            if (_constterm)
                _circ2.Evaluate(slam, clam, out BXc, out BYc, out BZc);

            if (_interpolate)
            {
                BX1 = (BX1 - BX0) / _dt0;
                BY1 = (BY1 - BY0) / _dt0;
                BZ1 = (BZ1 - BZ0) / _dt0;
            }

            BX0 += _t1 * BX1 + BXc;
            BY0 += _t1 * BY1 + BYc;
            BZ0 += _t1 * BZ1 + BZc;

            if (diffp)
            {
                Geocentric.Unrotate(M, BX1, BY1, BZ1, out Bxt, out Byt, out Bzt);
                Bxt *= -_a;
                Byt *= -_a;
                Bzt *= -_a;
            }
            else
            {
                Bxt = Byt = Bzt = double.NaN;
            }

            Geocentric.Unrotate(M, BX0, BY0, BZ0, out Bx, out By, out Bz);
            Bx *= -_a;
            By *= -_a;
            Bz *= -_a;
        }

        /// <summary>
        /// Evaluate the components of the geomagnetic field at a particular longitude.
        /// </summary>
        /// <param name="lon">longitude of the point (degrees).</param>
        /// <param name="Bx">the easterly component of the magnetic field (nanotesla).</param>
        /// <param name="By">the northerly component of the magnetic field (nanotesla).</param>
        /// <param name="Bz">the vertical (up) component of the magnetic field (nanotesla).</param>
        public void Evaluate(double lon, out double Bx, out double By, out double Bz)
            => Field(lon, false, out Bx, out By, out Bz, out _, out _, out _);

        /// <summary>
        /// Evaluate the components of the geomagnetic field and their time derivatives at a particular longitude.
        /// </summary>
        /// <param name="lon">longitude of the point (degrees).</param>
        /// <param name="Bx">the easterly component of the magnetic field (nanotesla).</param>
        /// <param name="By">the northerly component of the magnetic field (nanotesla).</param>
        /// <param name="Bz">the vertical (up) component of the magnetic field (nanotesla).</param>
        /// <param name="Bxt">the rate of change of <i>Bx</i> (nT/yr).</param>
        /// <param name="Byt">the rate of change of <i>By</i> (nT/yr).</param>
        /// <param name="Bzt">the rate of change of <i>Bz</i> (nT/yr).</param>
        public void Evaluate(double lon,
                             out double Bx, out double By, out double Bz,
                             out double Bxt, out double Byt, out double Bzt)
            => Field(lon, false, out Bx, out By, out Bz, out Bxt, out Byt, out Bzt);

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
