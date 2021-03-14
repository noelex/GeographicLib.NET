using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib
{
    /// <summary>
    /// Defines properties and methods of a geodesic-like object.
    /// </summary>
    public interface IGeodesicLike : IEllipsoid
    {
        /// <summary>
        /// Gets a value representing the total area of ellipsoid in meters^2.
        /// The area of a polygon encircling a pole can be found by adding <see cref="EllipsoidArea"/>/2 to the sum of <i>S12</i> for each side of the polygon.
        /// </summary>
        double EllipsoidArea { get; }

        #region General version of direct and inverse geodesic solution
        /// <summary>
        /// The general direct geodesic problem.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="arcmode">boolean flag determining the meaning of the <paramref name="s12_a12"/>.</param>
        /// <param name="s12_a12">
        /// if <paramref name="arcmode"/> is <see langword="false"/>, this is the distance between
        /// point 1 and point 2 (meters); otherwise it is the arc length between
        /// point 1 and point 2 (degrees); it can be negative.</param>
        /// <param name="outmask">
        /// a bitor'ed combination of <see cref="GeodesicFlags"/> values specifying
        /// which of the following parameters should be set.
        /// </param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="azi2">(forward) azimuth at point 2 (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters).</param>
        /// <param name="m12">reduced length of geodesic (meters).</param>
        /// <param name="M12">geodesic scale of point 2 relative to point 1 (dimensionless).</param>
        /// <param name="M21">geodesic scale of point 1 relative to point 2 (dimensionless).</param>
        /// <param name="S12">area under the geodesic (meters^2).</param>
        /// <returns><i>a12</i> arc length of between point 1 and point 2 (degrees).</returns>
        /// <remarks>
        /// The <see cref="GeodesicFlags"/> values possible for outmask are
        /// <list type="bullet">
        /// <item><i>outmask</i> |= <see cref="GeodesicFlags.Latitude"/> for the latitude <i>lat2</i>;</item>
        /// <item><i>outmask</i> |= <see cref="GeodesicFlags.Longitude"/> for the latitude <i>lon2</i>;</item>
        /// <item><i>outmask</i> |= <see cref="GeodesicFlags.Azimuth"/> for the latitude <i>azi2</i>;</item>
        /// <item><i>outmask</i> |= <see cref="GeodesicFlags.Distance"/> for the distance <i>s12</i>;</item>
        /// <item><i>outmask</i> |= <see cref="GeodesicFlags.ReducedLength"/> for the reduced length <i>m12</i>;</item>
        /// <item><i>outmask</i> |= <see cref="GeodesicFlags.GeodesicScale"/> for the geodesic scales <i>M12</i> and <i>M21</i>;</item>
        /// <item><i>outmask</i> |= <see cref="GeodesicFlags.Area"/> for the area <i>S12</i>;</item>
        /// <item><i>outmask</i> |= <see cref="GeodesicFlags.All"/> for all of the above;</item>
        /// <item><i>outmask</i> |= <see cref="GeodesicFlags.LongUnroll"/> to unroll <i>lon2</i> instead of wrapping it into the range [−180°, 180°].</item>
        /// </list>
        /// <para>
        /// The function value <i>a12</i> is always computed and returned and this equals <i>s12_a12</i> is arcmode is <see langword="true"/>.
        /// If outmask includes <see cref="GeodesicFlags.Distance"/> and arcmode is <see langword="false"/>, 
        /// then <i>s12</i> = <i>s12_a12</i>. It is not necessary to include <see cref="GeodesicFlags.DistanceIn"/> in outmask;
        /// this is automatically included is arcmode is <see langword="false"/>.
        /// </para>
        /// <para>
        /// With the <see cref="GeodesicFlags.LongUnroll"/> bit set, the quantity <i>lon2</i> − <i>lon1</i> indicates how many times 
        /// and in what sense the geodesic encircles the ellipsoid.
        /// </para>
        /// </remarks>
        double GenDirect(double lat1, double lon1, double azi1, bool arcmode, double s12_a12, GeodesicFlags outmask, out double lat2, out double lon2, out double azi2, out double s12, out double m12, out double M12, out double M21, out double S12);

        /// <summary>
        /// The general inverse geodesic calculation. 
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="outmask">
        /// a bitor'ed combination of <see cref="GeodesicFlags"/> values specifying
        /// which of the following parameters should be set.
        /// </param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="azi2">(forward) azimuth at point 2 (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters).</param>
        /// <param name="m12">reduced length of geodesic (meters).</param>
        /// <param name="M12">geodesic scale of point 2 relative to point 1 (dimensionless).</param>
        /// <param name="M21">geodesic scale of point 1 relative to point 2 (dimensionless).</param>
        /// <param name="S12">area under the geodesic (meters^2).</param>
        /// <returns><i>a12</i>, arc length of between point 1 and point 2 (degrees).</returns>
        double GenInverse(double lat1, double lon1, double lat2, double lon2, GeodesicFlags outmask, out double s12, out double azi1, out double azi2, out double m12, out double M12, out double M21, out double S12);

        #endregion
    }
}
