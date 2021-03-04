using System;
using System.Collections.Generic;
using System.Text;
using GeographicLib;

using static System.Math;
using static GeographicLib.MathEx;
using static GeographicLib.Macros;

namespace GeographicLib.Projections
{
    /// <summary>
    /// Gnomonic projection.
    /// </summary>
    /// <remarks>
    /// 
    ///  </remarks>
    public class Gnomonic : IEllipsoid
    {
        // numit_ increased from 10 to 20 to fix convergence failure with high
        // precision (e.g., GEOGRAPHICLIB_DIGITS=2000) calculations.  Reverse uses
        // Newton's method which converges quadratically and so numit_ = 10 would
        // normally be big enough.  However, since the Geodesic class is based on a
        // series it is of limited accuracy; in particular, the derivative rules
        // used by Reverse only hold approximately.  Consequently, after a few
        // iterations, the convergence in the Reverse falls back to improvements in
        // each step by a constant (albeit small) factor.
        private const int numit_ = 20;

        private static readonly double eps0_ = DBL_EPSILON, eps_ = 0.01 * Sqrt(eps0_);

        private readonly IGeodesic _earth;
        private readonly double _a, _f;

        /// <summary>
        /// Initialize a new <see cref="Gnomonic"/> instance with specified <see cref="Geodesic"/> instance.
        /// </summary>
        /// <param name="earth">the <see cref="IGeodesic"/> object to use for geodesic calculations.</param>
        public Gnomonic(IGeodesic earth)
        {
            _earth = earth;
            _a = _earth.EquatorialRadius;
            _f = _earth.Flattening;
        }

        /// <summary>
        /// Initialize a new <see cref="Gnomonic"/> instance with <see cref="Geodesic.WGS84"/>.
        /// </summary>
        public Gnomonic() : this(Geodesic.WGS84) { }

        /// <inheritdoc/>
        public double EquatorialRadius => _a;

        /// <inheritdoc/>
        public double Flattening => _f;

        /// <summary>
        /// Forward projection, from geographic to gnomonic.
        /// </summary>
        /// <param name="lat0">latitude of center point of projection (degrees).</param>
        /// <param name="lon0">longitude of center point of projection (degrees).</param>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <param name="azi">azimuth of geodesic at point (degrees).</param>
        /// <param name="rk">reciprocal of azimuthal scale at point.</param>
        /// <returns>
        /// <i>x</i>, easting of point (meters) and <i>y</i>, northing of point (meters).
        /// </returns>
        /// <remarks>
        /// <paramref name="lat0"/> and <paramref name="lat"/> should be in the range [−90°, 90°].
        /// The scale of the projection is 1/<paramref name="rk"/>^2 in the "radial" direction, <paramref name="azi"/> clockwise from true north,
        /// and is 1/<paramref name="rk"/> in the direction perpendicular to this. If the point lies "over the horizon", i.e.,
        /// if <paramref name="rk"/> ≤ 0, then <see cref="double.NaN"/>s are returned for <i>x</i> and <i>y</i> 
        /// (the correct values are returned for <paramref name="azi"/> and <paramref name="rk"/>).
        /// A call to <see cref="Forward(double, double, double, double, out double, out double)"/> followed by a call to
        /// <see cref="Reverse(double, double, double, double, out double, out double)"/> will return the 
        /// original (<paramref name="lat"/>, <paramref name="lon"/>) (to within roundoff) provided the point in not over the horizon.
        /// </remarks>
        public (double x, double y) Forward(double lat0, double lon0, double lat, double lon, out double azi, out double rk)
        {
            double x, y;

            _earth.GenInverse(lat0, lon0, lat, lon,
                              GeodesicFlags.Azimuth | GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale,
                              out _, out var azi0, out azi, out var m, out var M, out _, out _);
            rk = M;
            if (M <= 0)
                x = y = double.NaN;
            else
            {
                var rho = m / M;
                SinCosd(azi0, out x, out y);
                x *= rho; y *= rho;
            }

            return (x, y);
        }

        /// <summary>
        /// <see cref="Forward(double, double, double, double, out double, out double)"/> without returning the azimuth and scale.
        /// </summary>
        /// <param name="lat0">latitude of center point of projection (degrees).</param>
        /// <param name="lon0">longitude of center point of projection (degrees).</param>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <returns>
        /// <i>x</i>, easting of point (meters) and <i>y</i>, northing of point (meters).
        /// </returns>
        public (double x, double y) Forward(double lat0, double lon0, double lat, double lon)
            => Forward(lat0, lon0, lat, lon, out _, out _);

        /// <summary>
        /// Reverse projection, from gnomonic to geographic.
        /// </summary>
        /// <param name="lat0">latitude of center point of projection (degrees).</param>
        /// <param name="lon0">longitude of center point of projection (degrees).</param>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <param name="azi">azimuth of geodesic at point (degrees).</param>
        /// <param name="rk">reciprocal of azimuthal scale at point.</param>
        /// <returns>
        /// <i>lat</i>, latitude of point (degrees) and <i>lon</i>, longitude of point (degrees).
        /// </returns>
        /// <remarks>
        /// <paramref name="lat0"/> should be in the range [−90°, 90°].
        /// <i>lat</i> will be in the range [−90°, 90°] and <i>lon</i> will be in the range [−180°, 180°].
        /// The scale of the projection is 1/<paramref name="rk"/>^2 in the "radial" direction, <paramref name="azi"/> clockwise from true north,
        /// and is 1/<paramref name="rk"/> in the direction perpendicular to this. Even though all inputs should return a valid <i>lat</i> and <i>lon</i>,
        /// it's possible that the procedure fails to converge for very large <paramref name="x"/> or <paramref name="y"/>;
        /// in this case <see cref="double.NaN"/>s are returned for all the output arguments.
        /// A call to <see cref="Reverse(double, double, double, double, out double, out double)"/> followed by a call to
        /// <see cref="Forward(double, double, double, double, out double, out double)"/> 
        /// will return the original (<paramref name="x"/>, <paramref name="y"/>) (to roundoff).
        /// </remarks>
        public (double lat, double lon) Reverse(double lat0, double lon0, double x, double y, out double azi, out double rk)
        {
            double
                lat, lon,
                azi0 = Atan2d(x, y),
                rho = Hypot(x, y),
                s = _a * Atan(rho / _a);

            bool little = rho <= _a;
            if (!little)
                rho = 1 / rho;

            var line=_earth.Line(lat0, lon0, azi0,
                                          GeodesicFlags.Latitude | GeodesicFlags.Longitude |
                                          GeodesicFlags.Azimuth | GeodesicFlags.DistanceIn |
                                          GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale);
            int count = numit_, trip = 0;
            double lat1=0, lon1 = 0, azi1 = 0, M = 0;

            while (count--!=0 || GEOGRAPHICLIB_PANIC)
            {
                line.Position(s, out lat1, out lon1, out azi1, out var m, out M, out _);
                if (trip!=0)
                    break;
                // If little, solve rho(s) = rho with drho(s)/ds = 1/M^2
                // else solve 1/rho(s) = 1/rho with d(1/rho(s))/ds = -1/m^2
                var ds = little ? (m - rho * M) * M : (rho * m - M) * m;
                s -= ds;
                // Reversed test to allow escape with NaNs
                if (!(Abs(ds) >= eps_ * _a))
                    ++trip;
            }

            if (trip!=0)
            {
                (lat, lon, azi, rk) = (lat1, lon1, azi1, M);
            }
            else
                lat = lon = azi = rk = double.NaN;

            return (lat, lon);
        }

        /// <summary>
        /// <see cref="Reverse(double, double, double, double, out double, out double)"/> without returning the azimuth and scale.
        /// </summary>
        /// <param name="lat0">latitude of center point of projection (degrees).</param>
        /// <param name="lon0">longitude of center point of projection (degrees).</param>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <returns>
        /// <i>lat</i>, latitude of point (degrees) and <i>lon</i>, longitude of point (degrees).
        /// </returns>
        public (double lat, double lon) Reverse(double lat0, double lon0, double x, double y)
            => Reverse(lat0, lon0, x, y, out _, out _);
    }
}
