using System;
using System.Collections.Generic;
using System.Text;
using GeographicLib;

using static System.Math;
using static GeographicLib.MathEx;

namespace GeographicLib.Projections
{
    /// <summary>
    /// Azimuthal equidistant projection.
    /// </summary>
    /// <remarks>
    /// Azimuthal equidistant projection centered at an arbitrary position on the ellipsoid.
    /// For a point in projected space (<i>x</i>, <i>y</i>), the geodesic distance from the center position is hypot(<i>x</i>, <i>y</i>) 
    /// and the azimuth of the geodesic from the center point is atan2(<i>x</i>, <i>y</i>).
    /// The <see cref="Forward(double, double, double, double, out double, out double)"/> and 
    /// <see cref="Reverse(double, double, double, double, out double, out double)"/> methods also return the azimuth <i>azi</i> of the geodesic at (<i>x</i>, <i>y</i>)
    /// and reciprocal scale <i>rk</i> in the azimuthal direction which, together with the basic properties of the projection, serve to specify
    /// completely the local affine transformation between geographic and projected coordinates.
    /// <para>
    /// The conversions all take place using a <see cref="Geodesic"/> object (by default <see cref="Geodesic.WGS84"/>).
    /// For more information on geodesics see <a href="https://geographiclib.sourceforge.io/html/geodesic.html">Geodesics on an ellipsoid of revolution</a>.
    /// </para>
    /// </remarks>
    public class AzimuthalEquidistant : IEllipsoid
    {
        private readonly double eps_ = 0.01 * Sqrt(DBL_MIN);
        private readonly IGeodesic _earth;

        /// <summary>
        /// Constructor for <see cref="AzimuthalEquidistant"/>.
        /// </summary>
        /// <param name="earth">the <see cref="IGeodesic"/> object to use for geodesic calculations.</param>
        public AzimuthalEquidistant(IGeodesic earth) => _earth = earth;

        /// <summary>
        /// Initialize a <see cref="AzimuthalEquidistant"/> with <see cref="Geodesic.WGS84"/>.
        /// </summary>
        public AzimuthalEquidistant() : this(Geodesic.WGS84) { }

        /// <summary>
        /// Gets a value representing the equatorial radius (<i>a</i>) of the ellipsoid.
        /// </summary>
        public double EquatorialRadius => _earth.EquatorialRadius;

        /// <summary>
        /// Gets a value representing the flattening (<i>f</i>) of the ellipsoid.
        /// </summary>
        public double Flattening => _earth.Flattening;

        /// <summary>
        /// Forward projection, from geographic to azimuthal equidistant.
        /// </summary>
        /// <param name="lat0">latitude of center point of projection (degrees).</param>
        /// <param name="lon0">longitude of center point of projection (degrees).</param>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <param name="azi">azimuth of easting direction at point (degrees).</param>
        /// <param name="rk">reciprocal of azimuthal northing scale at point.</param>
        /// <returns>
        /// <i>x</i>, easting of point and <i>y</i>, northing of point, in meters.
        /// </returns>
        /// <remarks>
        /// <paramref name="lat0"/> and <paramref name="lat"/>  should be in the range [−90°, 90°].
        /// The scale of the projection is 1 in the "radial" direction, <paramref name="azi"/> clockwise from true north,
        /// and is 1/<paramref name="rk"/> in the direction perpendicular to this.
        /// A call to <see cref="Forward(double, double, double, double, out double, out double)"/> 
        /// followed by a call to <see cref="Reverse(double, double, double, double, out double, out double)"/> 
        /// will return the original (<paramref name="lat"/>, <paramref name="lon"/>) (to within roundoff).
        /// </remarks>
        public (double x, double y) Forward(double lat0, double lon0, double lat, double lon, out double azi, out double rk)
        {
            var sig = _earth.Inverse(lat0, lon0, lat, lon, out var s, out var azi0, out azi, out var m);
            SinCosd(azi0, out var x, out var y);
            x *= s; y *= s;
            rk = !(sig <= eps_) ? m / s : 1;

            return (x, y);
        }

        /// <summary>
        /// Forward without returning the azimuth and scale.
        /// </summary>
        /// <param name="lat0">latitude of center point of projection (degrees).</param>
        /// <param name="lon0">longitude of center point of projection (degrees).</param>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <returns>
        /// <i>x</i>, easting of point and <i>y</i>, northing of point, in meters.
        /// </returns>
        /// <remarks>
        /// <paramref name="lat0"/> and <paramref name="lat"/>  should be in the range [−90°, 90°].
        /// The scale of the projection is 1 in the "radial" direction, <i>azi</i> clockwise from true north,
        /// and is 1/<i>rk</i> in the direction perpendicular to this.
        /// A call to <see cref="Forward(double, double, double, double)"/> followed by a call to
        /// <see cref="Reverse(double, double, double, double)"/> 
        /// will return the original (<paramref name="lat"/>, <paramref name="lon"/>) (to within roundoff).
        /// </remarks>
        public (double x, double y) Forward(double lat0, double lon0, double lat, double lon)
            => Forward(lat0, lon0, lat, lon, out _, out _);

        /// <summary>
        /// Reverse projection, from azimuthal equidistant to geographic.
        /// </summary>
        /// <param name="lat0">latitude of center point of projection (degrees).</param>
        /// <param name="lon0">longitude of center point of projection (degrees).</param>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <param name="azi">azimuth of easting direction at point (degrees).</param>
        /// <param name="rk">reciprocal of azimuthal northing scale at point.</param>
        /// <returns>
        /// <i>lat</i>, latitude of point and <i>lon</i>, longitude of point, in degress.
        /// </returns>
        /// <remarks>
        /// <paramref name="lat0"/> should be in the range [−90°, 90°].
        /// <i>lat</i> will be in the range [−90°, 90°] and <i>lon</i> will be in the range [−180°, 180°].
        /// The scale of the projection is 1 in the "radial" direction, <paramref name="azi"/> clockwise from true north,
        /// and is 1/<paramref name="rk"/> in the direction perpendicular to this.
        /// A call to <see cref="Reverse(double, double, double, double, out double, out double)"/> 
        /// followed by a call to <see cref="Forward(double, double, double, double, out double, out double)"/> 
        /// will return the original (<paramref name="x"/>, <paramref name="y"/>) (to within roundoff) only if the geodesic to 
        /// (<paramref name="x"/>, <paramref name="y"/>) is a shortes path.
        /// </remarks>
        public (double lat, double lon) Reverse(double lat0, double lon0, double x, double y, out double azi, out double rk)
        {
            double
              azi0 = Atan2d(x, y),
              s = Hypot(x, y);

            var sig = _earth.Direct(lat0, lon0, azi0, s, out var lat, out var lon, out azi, out var m);
            rk = !(sig <= eps_) ? m / s : 1;

            return (lat, lon);
        }

        /// <summary>
        /// Reverse without returning the azimuth and scale.
        /// </summary>
        /// <param name="lat0">latitude of center point of projection (degrees).</param>
        /// <param name="lon0">longitude of center point of projection (degrees).</param>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <returns>
        /// <i>lat</i>, latitude of point and <i>lon</i>, longitude of point, in degress.
        /// </returns>
        /// <remarks>
        /// <paramref name="lat0"/> should be in the range [−90°, 90°].
        /// <i>lat</i> will be in the range [−90°, 90°] and <i>lon</i> will be in the range [−180°, 180°].
        /// The scale of the projection is 1 in the "radial" direction, <i>azi</i> clockwise from true north,
        /// and is 1/<i>rk</i> in the direction perpendicular to this.
        /// A call to <see cref="Reverse(double, double, double, double)"/> 
        /// followed by a call to <see cref="Forward(double, double, double, double)"/> 
        /// will return the original (<paramref name="x"/>, <paramref name="y"/>) (to within roundoff) only if the geodesic to 
        /// (<paramref name="x"/>, <paramref name="y"/>) is a shortes path.
        /// </remarks>
        public (double lat, double lon) Reverse(double lat0, double lon0, double x, double y)
            => Reverse(lat0, lon0, x, y, out _, out _);

    }
}
