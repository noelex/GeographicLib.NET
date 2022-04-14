using System.Collections.Generic;

namespace GeographicLib
{
    /// <summary>
    /// Computes polygon areas.
    /// </summary>
    public interface IPolygonArea
    {
        /// <summary>
        /// Report the number of points currently in the polygon or polyline.
        /// </summary>
        /// <remarks>
        /// If no points have been added, then 0 is returned.
        /// </remarks>
        int Count { get; }

        /// <summary>
        /// Report whether the current object is a polygon or a polyline.
        /// </summary>
        bool IsPolyline { get; }

        /// <summary>
        /// Gets a value representing the previous vertex added to the polygon or polyline.
        /// </summary>
        /// <remarks>
        /// If no points have been added, then <see cref="double.NaN"/>s are returned. Otherwise, <i>lon</i> will be in the range [−180°, 180°].
        /// </remarks>
        (double lat, double lon) CurrentPoint { get; }

        /// <summary>
        /// Add an edge to the polygon or polyline.
        /// </summary>
        /// <param name="azi">azimuth at current point (degrees).</param>
        /// <param name="s">distance from current point to next point (meters).</param>
        /// <remarks>
        /// This does nothing if no points have been added yet.
        /// Use <see cref="CurrentPoint"/> to determine the position of the new vertex.
        /// </remarks>
        void AddEdge(double azi, double s);

        /// <summary>
        /// Add a point to the polygon or polyline.
        /// </summary>
        /// <param name="lat">the latitude of the point (degrees).</param>
        /// <param name="lon">the longitude of the point (degrees).</param>
        /// <remarks>
        /// <paramref name="lat"/> should be in the range [−90°, 90°].
        /// </remarks>
        void AddPoint(double lat, double lon);

        /// <summary>
        /// Add a point to the polygon or polyline.
        /// </summary>
        /// <param name="coords">point to add.</param>
        void AddPoint(GeoCoords coords);

        /// <summary>
        /// Add points to the polygon or polyline.
        /// </summary>
        /// <param name="coords">points to add.</param>
        void AddPoints(IEnumerable<GeoCoords> coords);

        /// <summary>
        /// Clear current instance of <see cref="PolygonArea{T}"/>, allowing a new polygon to be started.
        /// </summary>
        void Clear();

        /// <summary>
        /// Return the results so far.
        /// </summary>
        /// <param name="reverse">if <see langword="true"/> then clockwise (instead of counter-clockwise) traversal counts as a positive area.</param>
        /// <param name="sign">if <see langword="true"/> then return a signed result for the area if the polygon is traversed in the "wrong" direction
        /// instead of returning the area for the rest of the earth.</param>
        /// <returns>
        /// <list type="bullet">
        /// <item>
        /// <i>points</i>, the number of points.
        /// </item>
        /// <item>
        /// <i>perimeter</i>, the perimeter of the polygon or length of the polyline (meters).
        /// </item>
        /// <item>
        /// <i>area</i>, the area of the polygon (meters^2); only set if polyline is <see langword="false"/> in the constructor.
        /// </item>
        /// </list>
        /// </returns>
        (int points, double perimeter, double area) Compute(bool reverse, bool sign);

        /// <summary>
        /// Return the results assuming a tentative final test point is added via an azimuth and distance;
        /// however, the data for the test point is not saved.
        /// This lets you report a running result for the perimeter and area as the user moves the mouse cursor.
        /// Ordinary floating point arithmetic is used to accumulate the data for the test point;
        /// thus the area and perimeter returned are less accurate than if
        /// <see cref="AddEdge(double, double)"/> and <see cref="Compute(bool, bool)"/> are used.
        /// </summary>
        /// <param name="azi">azimuth at current point (degrees).</param>
        /// <param name="s">distance from current point to final test point (meters).</param>
        /// <param name="reverse">if <see langword="true"/> then clockwise (instead of counter-clockwise) traversal counts as a positive area.</param>
        /// <param name="sign">if <see langword="true"/> then return a signed result for the area if the polygon is traversed in the "wrong" direction
        /// instead of returning the area for the rest of the earth.</param>
        /// <returns>
        /// <list type="bullet">
        /// <item>
        /// <i>points</i>, the number of points.
        /// </item>
        /// <item>
        /// <i>perimeter</i>, the perimeter of the polygon or length of the polyline (meters).
        /// </item>
        /// <item>
        /// <i>area</i>, the area of the polygon (meters^2); only set if polyline is <see langword="false"/> in the constructor.
        /// </item>
        /// </list>
        /// </returns>
        (int points, double perimeter, double area) TestEdge(double azi, double s, bool reverse, bool sign);

        /// <summary>
        /// Return the results assuming a tentative final test point is added; however, the data for the test point is not saved.
        /// This lets you report a running result for the perimeter and area as the user moves the mouse cursor.
        /// Ordinary floating point arithmetic is used to accumulate the data for the test point;
        /// thus the area and perimeter returned are less accurate than if
        /// <see cref="AddPoint(double, double)"/> and <see cref="Compute(bool, bool)"/> are used.
        /// </summary>
        /// <param name="lat">the latitude of the test point (degrees).</param>
        /// <param name="lon">the longitude of the test point (degrees).</param>
        /// <param name="reverse">if <see langword="true"/> then clockwise (instead of counter-clockwise) traversal counts as a positive area.</param>
        /// <param name="sign">if <see langword="true"/> then return a signed result for the area if the polygon is traversed in the "wrong" direction
        /// instead of returning the area for the rest of the earth.</param>
        /// <returns>
        /// <list type="bullet">
        /// <item>
        /// <i>points</i>, the number of points.
        /// </item>
        /// <item>
        /// <i>perimeter</i>, the perimeter of the polygon or length of the polyline (meters).
        /// </item>
        /// <item>
        /// <i>area</i>, the area of the polygon (meters^2); only set if polyline is <see langword="false"/> in the constructor.
        /// </item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// <paramref name="lat"/> should be in the range [−90°, 90°].
        /// </remarks>
        (int points, double perimeter, double area) TestPoint(double lat, double lon, bool reverse, bool sign);

        /// <summary>
        /// Return the results assuming a tentative final test point is added; however, the data for the test point is not saved.
        /// This lets you report a running result for the perimeter and area as the user moves the mouse cursor.
        /// Ordinary floating point arithmetic is used to accumulate the data for the test point;
        /// thus the area and perimeter returned are less accurate than if
        /// <see cref="AddPoint(GeoCoords)"/> and <see cref="Compute(bool, bool)"/> are used.
        /// </summary>
        /// <param name="coords">point to test.</param>
        /// <param name="reverse">if <see langword="true"/> then clockwise (instead of counter-clockwise) traversal counts as a positive area.</param>
        /// <param name="sign">if <see langword="true"/> then return a signed result for the area if the polygon is traversed in the "wrong" direction
        /// instead of returning the area for the rest of the earth.</param>
        /// <returns>
        /// <list type="bullet">
        /// <item>
        /// <i>points</i>, the number of points.
        /// </item>
        /// <item>
        /// <i>perimeter</i>, the perimeter of the polygon or length of the polyline (meters).
        /// </item>
        /// <item>
        /// <i>area</i>, the area of the polygon (meters^2); only set if polyline is <see langword="false"/> in the constructor.
        /// </item>
        /// </list>
        /// </returns>
        (int points, double perimeter, double area) TestPoint(GeoCoords coords, bool reverse, bool sign);
    }
}