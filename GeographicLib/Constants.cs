using System;
using System.Collections.Generic;
using System.Text;

namespace GeographicLib
{
    /// <summary>
    /// Constants need by GeographicLib.
    /// </summary>
    /// <remarks>
    /// Define constants specifying the WGS84 ellipsoid, the UTM and UPS
    /// projections, and various unit conversions.
    /// </remarks>
    public static class Constants
    {
        /// <summary>
        /// A synonym for <see cref="MathEx.Degree"/>.
        /// </summary>
        public const double Degree = MathEx.Degree;

        /// <summary>
        /// The number of radians in an arcminute.
        /// </summary>
        public const double Arcminute = Degree / MathEx.dm;

        /// <summary>
        /// The number of radians in an arcsecond.
        /// </summary>
        public const double Arcsecond = Degree / MathEx.ds;

        #region Ellipoid parameters

        /// <summary>
        /// The equatorial radius of WGS84 ellipsoid (6378137 m).
        /// </summary>
        public const double WGS84_a = 6378137 * Meter;

        /// <summary>
        /// The flattening of WGS84 ellipsoid (1/298.257223563).
        /// </summary>
        // Evaluating this as 1000000000 / T(298257223563LL) reduces the
        // round-off error by about 10%.  However, expressing the flattening as
        // 1/298.257223563 is well ingrained.
        public const double WGS84_f = 1 / (298257223563d / 1000000000);

        /// <summary>
        /// The gravitational constant of the WGS84 ellipsoid, <i>GM</i> , in m^3/s^2
        /// </summary>
        public const double WGS84_GM = 3986004d * 100000000 + 41800000;

        /// <summary>
        /// The angular velocity of the WGS84 ellipsoid, <i>ω</i>, in rad/s
        /// </summary>
        public const double WGS84_omega = 7292115 / (1000000d * 100000);

        /// <summary>
        /// The equatorial radius of GRS80 ellipsoid (6378137 m).
        /// </summary>
        public const double GRS80_a = 6378137 * Meter;

        /// <summary>
        /// The gravitational constant of the GRS80 ellipsoid, <i>GM</i> , in m^3/s^2
        /// </summary>
        public const double GRS80_GM = 3986005d * 100000000;

        /// <summary>
        /// The angular velocity of the GRS80 ellipsoid, <i>ω</i>, in rad/s.
        /// </summary>
        /// <remarks>
        /// This is about 2π * 366.25 / (365.25 * 24 * 3600) rad/s.
        /// 365.25 is the number of days in a Julian year and
        /// 365.35/366.25 converts from solar days to sidereal days.Using the
        /// number of days in a Gregorian year(365.2425) results in a worse
        /// approximation(because the Gregorian year includes the precession of the
        /// earth's axis).
        /// </remarks>
        public const double GRS80_omega = 7292115 / (1000000d * 100000);

        /// <summary>
        /// Dynamical form factor of the GRS80 ellipsoid, in <i>J^2</i>.
        /// </summary>
        public const double GRS80_J2 = 108263d / 100000000;

        /// <summary>
        /// The central scale factor for UTM (0.9996).
        /// </summary>
        public const double UTM_k0 = 9996d / 10000;

        /// <summary>
        /// The central scale factor for UPS (0.994).
        /// </summary>
        public const double UPS_k0 = 994d / 1000;

        #endregion

        #region SI Units

        /// <summary>
        /// The number of meters in a meter.
        /// </summary>
        public const double Meter = 1;

        /// <summary>
        /// The number of meters in a kilometer.
        /// </summary>
        public const double Kilometer = 1000 * Meter;

        /// <summary>
        /// The number of meters in a nautical mile (approximately 1 arcminute).
        /// </summary>
        public const double NauticalMile = 1852 * Meter;

        /// <summary>
        /// The number of square meters in a square meter
        /// </summary>
        public const double SquareMeter = Meter * Meter;

        /// <summary>
        /// The number of square meters in a hectare.
        /// </summary>
        public const double Hectare = 10000 * SquareMeter;

        /// <summary>
        /// The number of square meters in a square kilometer.
        /// </summary>
        public const double SquareKilometer = Kilometer * Kilometer;

        /// <summary>
        /// The number of square meters in a square nautical mile.
        /// </summary>
        public const double SquareNauticleMile = NauticalMile * NauticalMile;

        #endregion

        #region Anachronistic British units

        /// <summary>
        /// The number of meters in an international foot.
        /// </summary>
        public const double Foot = (254d * 12) / 10000 * Meter;

        /// <summary>
        /// The number of meters in a yard.
        /// </summary>
        public const double Yard = 3 * Foot;

        /// <summary>
        /// The number of meters in a fathom.
        /// </summary>
        public const double Fathom = 2 * Yard;

        /// <summary>
        /// The number of meters in a chain.
        /// </summary>
        public const double Chain = 22 * Yard;

        /// <summary>
        /// The number of meters in a furlong.
        /// </summary>
        public const double Furlong = 10 * Chain;

        /// <summary>
        /// The number of meters in a statute mile.
        /// </summary>
        public const double Mile = 8 * Furlong;

        /// <summary>
        /// The number of square meters in an acre.
        /// </summary>
        public const double Acre= Chain * Furlong;

        /// <summary>
        /// The number of square meters in a square statute mile.
        /// </summary>
        public const double SquareMile = Mile * Mile;

        #endregion

        #region Anachronistic US units

        /// <summary>
        /// The number of meters in a US survey foot.
        /// </summary>
        public const double SurveyFoot = 1200d / 3937 * Meter;

        #endregion
    }
}
