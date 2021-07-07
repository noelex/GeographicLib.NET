using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using static System.Math;
using static GeographicLib.MathEx;

namespace GeographicLib
{
    /// <summary>
    /// Computes polygon areas.
    /// </summary>
    /// <typeparam name="T">the geodesic class to use.</typeparam>
    /// <remarks>
    /// This computes the area of a polygon whose edges are geodesics using the method given in Section 6 of
    /// <list type="bullet">
    /// <item>
    /// C. F. F. Karney,
    /// <a href="https://doi.org/10.1007/s00190-012-0578-z">Algorithms for geodesics</a>,
    /// J. Geodesy 87, 43–55 (2013);
    /// DOI: 10.1007/s00190-012-0578-z;
    /// addenda: geod-addenda.html.
    /// </item>
    /// </list>
    /// Arbitrarily complex polygons are allowed.
    /// In the case self-intersecting of polygons the area is accumulated "algebraically", e.g., the areas of the 2 loops in a figure-8 polygon
    /// will partially cancel.
    /// <para>
    /// This class lets you add vertices and edges one at a time to the polygon.
    /// The sequence must start with a vertex and thereafter vertices and edges can be added in any order.
    /// Any vertex after the first creates a new edge which is the shortest geodesic from the previous vertex.
    /// In some cases there may be two or many such shortest geodesics and the area is then not uniquely defined.
    /// In this case, either add an intermediate vertex or add the edge as an edge (by defining its direction and length).
    /// </para>
    /// <para>
    /// The area and perimeter are accumulated at two times the standard floating point precision to guard against the loss of accuracy
    /// with many-sided polygons. At any point you can ask for the perimeter and area so far. There's an option to treat the points as defining
    /// a polyline instead of a polygon; in that case, only the perimeter is computed.
    /// </para>
    /// <para>
    /// This is a templated class to allow it to be used with <see cref="Geodesic"/>, <see cref="GeodesicExact"/>, and <see cref="Rhumb"/>.
    /// <see cref="PolygonArea"/>, <see cref="PolygonAreaExact"/>, and <see cref="PolygonAreaRhumb"/> are typedefs for these cases.
    /// </para>
    /// <para>
    /// For <see cref="PolygonArea"/> (edges defined by <see cref="Geodesic"/>), an upper bound on the error is about 0.1 m^2 per vertex.
    /// However this is a wildly pessimistic estimate in most cases.
    /// A more realistic estimate of the error is given by a test involving 10^7 approximately regular polygons on the WGS84 ellipsoid.
    /// The centers and the orientations of the polygons were uniformly distributed, the number of vertices was log-uniformly distributed in [3, 300],
    /// and the center to vertex distance log-uniformly distributed in [0.1 m, 9000 km].
    /// </para>
    /// <para>
    /// Using double precision (the standard precision for GeographicLib and GeographicLib.NET), the maximum error in the perimeter was 200 nm,
    /// and the maximum error in the area was
    /// <code>
    /// <para>0.0013 m^2 for perimeter &lt; 10 km</para>
    /// <para>0.0070 m^2 for perimeter &lt; 100 km</para>
    /// <para>0.070 m^2 for perimeter &lt; 1000 km</para>
    /// <para>0.11 m^2 for all perimeters</para>
    /// </code>
    /// </para>
    /// <para>
    /// The errors are given in terms of the perimeter, because it is expected that the errors depend mainly on the number of edges and the edge lengths.
    /// </para>
    /// <para>
    /// Using long doubles (GEOGRPAHICLIB_PRECISION = 3, which is not supported by GeographicLib.NET), the maximum error in the perimeter was 200 pm,
    /// and the maximum error in the area was
    /// </para>
    /// <code>
    /// <para>0.7 mm^2 for perimeter &lt; 10 km</para>
    /// <para>3.2 mm^2 for perimeter &lt; 100 km</para>
    /// <para>21 mm^2 for perimeter &lt; 1000 km</para>
    /// <para>45 mm^2 for all perimeters</para>
    /// </code>
    /// </remarks>
    public class PolygonArea<T> where T : IGeodesicLike
    {
        private readonly T _earth;

        private double _area0;
        private bool _polyline;
        private GeodesicFlags _mask;
        private int _num;
        private int _crossings;

        private Accumulator _areasum, _perimetersum;

        private double _lat0, _lon0, _lat1, _lon1;

        /// <summary>
        /// Constructor for <see cref="PolygonArea{T}"/>.
        /// </summary>
        /// <param name="earth">the <see cref="IGeodesicLike"/> object to use for geodesic calculations.</param>
        /// <param name="polyline">
        /// if <see langword="true"/> that treat the points as defining a polyline instead of a polygon (default = <see langword="false"/>).
        /// </param>
        public PolygonArea(T earth, bool polyline = false)
        {
            _earth = earth;
            _area0 = _earth.EllipsoidArea;
            _polyline = polyline;
            _mask = GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Distance |
                (_polyline ? GeodesicFlags.None : GeodesicFlags.Area | GeodesicFlags.LongUnroll);

            Clear();
        }

        private static int Transit(double lon1, double lon2)
        {
            // Return 1 or -1 if crossing prime meridian in east or west direction.
            // Otherwise return zero.
            // Compute lon12 the same way as Geodesic::Inverse.
            lon1 = AngNormalize(lon1);
            lon2 = AngNormalize(lon2);
            var lon12 = AngDiff(lon1, lon2);
            // Treat 0 as negative in these tests.  This balances +/- 180 being
            // treated as positive, i.e., +180.
            int cross =
              lon1 <= 0 && lon2 > 0 && lon12 > 0 ? 1 :
              (lon2 <= 0 && lon1 > 0 && lon12 < 0 ? -1 : 0);
            return cross;
        }

        // an alternate version of transit to deal with longitudes in the direct
        // problem.
        private static int TransitDirect(double lon1, double lon2)
        {
            // Compute exactly the parity of
            //   int(ceil(lon2 / 360)) - int(ceil(lon1 / 360))
            lon1 = IEEERemainder(lon1, 720);
            lon2 = IEEERemainder(lon2, 720);
            return ((lon2 <= 0 && lon2 > -360 ? 1 : 0) -
                     (lon1 <= 0 && lon1 > -360 ? 1 : 0));
        }

        private void Remainder(ref Accumulator a) => a.Remainder(_area0);

        private void Remainder(ref double a) => a = IEEERemainder(a, _area0);

        private void AreaReduce(ref Accumulator area, int crossings, bool reverse, bool sign)
        {
            Remainder(ref area);
            if ((crossings & 1) !=0) area += (area < 0 ? 1 : -1) * _area0 / 2;
            // area is with the clockwise sense.  If !reverse convert to
            // counter-clockwise convention.
            if (!reverse) area *= -1;
            // If sign put area in (-_area0/2, _area0/2], else put area in [0, _area0)
            if (sign)
            {
                if (area > _area0 / 2)
                    area -= _area0;
                else if (area <= -_area0 / 2)
                    area += _area0;
            }
            else
            {
                if (area >= _area0)
                    area -= _area0;
                else if (area < 0)
                    area += _area0;
            }
        }

        private void AreaReduce(ref double area, int crossings, bool reverse, bool sign)
        {
            Remainder(ref area);
            if ((crossings & 1) != 0) area += (area < 0 ? 1 : -1) * _area0 / 2;
            // area is with the clockwise sense.  If !reverse convert to
            // counter-clockwise convention.
            if (!reverse) area *= -1;
            // If sign put area in (-_area0/2, _area0/2], else put area in [0, _area0)
            if (sign)
            {
                if (area > _area0 / 2)
                    area -= _area0;
                else if (area <= -_area0 / 2)
                    area += _area0;
            }
            else
            {
                if (area >= _area0)
                    area -= _area0;
                else if (area < 0)
                    area += _area0;
            }
        }

        /// <summary>
        /// Clear current instance of <see cref="PolygonArea{T}"/>, allowing a new polygon to be started.
        /// </summary>
        public void Clear()
        {
            _num = 0;
            _crossings = 0;
            _areasum = 0;
            _perimetersum = 0;
            _lat0 = _lon0 = _lat1 = _lon1 = double.NaN;
        }

        /// <summary>
        /// Add a point to the polygon or polyline.
        /// </summary>
        /// <param name="lat">the latitude of the point (degrees).</param>
        /// <param name="lon">the longitude of the point (degrees).</param>
        /// <remarks>
        /// <paramref name="lat"/> should be in the range [−90°, 90°].
        /// </remarks>
        public void AddPoint(double lat, double lon)
        {
            lat = LatFix(lat);
            lon = AngNormalize(lon);
            if (_num == 0)
            {
                _lat0 = _lat1 = lat;
                _lon0 = _lon1 = lon;
            }
            else
            {
                _earth.GenInverse(_lat1, _lon1, lat, lon, _mask,
                                  out var s12, out _, out _, out _, out _, out _, out var S12);
                _perimetersum += s12;
                if (!_polyline)
                {
                    _areasum += S12;
                    _crossings += Transit(_lon1, lon);
                }
                _lat1 = lat; _lon1 = lon;
            }
            ++_num;
        }

        /// <summary>
        /// Add a point to the polygon or polyline.
        /// </summary>
        /// <param name="coords">point to add.</param>
        public void AddPoint(GeoCoords coords) => AddPoint(coords.Latitude, coords.Longitude);

        /// <summary>
        /// Add points to the polygon or polyline.
        /// </summary>
        /// <param name="coords">points to add.</param>
        public void AddPoints(IEnumerable<GeoCoords> coords)
        {
            foreach(var coord in coords)
            {
                AddPoint(coord);
            }
        }

        /// <summary>
        /// Add an edge to the polygon or polyline.
        /// </summary>
        /// <param name="azi">azimuth at current point (degrees).</param>
        /// <param name="s">distance from current point to next point (meters).</param>
        /// <remarks>
        /// This does nothing if no points have been added yet.
        /// Use <see cref="CurrentPoint"/> to determine the position of the new vertex.
        /// </remarks>
        public void AddEdge(double azi, double s)
        {
            if (_num!=0)
            {                 // Do nothing if _num is zero
                _earth.GenDirect(_lat1, _lon1, azi, false, s, _mask,
                                 out var lat, out var lon, out _, out _, out _, out _, out _, out var S12);
                _perimetersum += s;
                if (!_polyline)
                {
                    _areasum += S12;
                    _crossings += TransitDirect(_lon1, lon);
                    lon = AngNormalize(lon);
                }
                _lat1 = lat; _lon1 = lon;
                ++_num;
            }
        }

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
        public (int points, double perimeter, double area) Compute(bool reverse, bool sign)
        {
            double perimeter, area=0;
            if (_num < 2)
            {
                perimeter = 0;
                if (!_polyline)
                    area = 0;
                return (_num, perimeter, area);
            }
            if (_polyline)
            {
                perimeter = _perimetersum;
                return (_num, perimeter, area);
            }
            _earth.GenInverse(_lat1, _lon1, _lat0, _lon0, _mask,
                              out var s12, out _, out _, out _, out _, out _, out var S12);
            perimeter = _perimetersum.Sum(s12);
            var tempsum = _areasum;
            tempsum += S12;
            int crossings = _crossings + Transit(_lon1, _lon0);
            AreaReduce(ref tempsum, crossings, reverse, sign);
            area = 0 + tempsum;
            return (_num, perimeter, area);
        }

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
        public (int points, double perimeter, double area) TestPoint(double lat, double lon, bool reverse, bool sign)
        {
            double perimeter, area = 0;
            if (_num == 0)
            {
                perimeter = 0;
                if (!_polyline)
                    area = 0;
                return (1, perimeter, area);
            }
            perimeter = _perimetersum;
            double tempsum = _polyline ? 0 : _areasum;
            int crossings = _crossings;
            var num = _num + 1;
            for (int i = 0; i < (_polyline ? 1 : 2); ++i)
            {
                _earth.GenInverse(i == 0 ? _lat1 : lat, i == 0 ? _lon1 : lon,
                                  i != 0 ? _lat0 : lat, i != 0 ? _lon0 : lon,
                                  _mask, out var s12, out _, out _, out _, out _, out _, out var S12);
                perimeter += s12;
                if (!_polyline)
                {
                    tempsum += S12;
                    crossings += Transit(i == 0 ? _lon1 : lon,
                                         i != 0 ? _lon0 : lon);
                }
            }

            if (_polyline)
                return (num, perimeter, area);

            AreaReduce(ref tempsum, crossings, reverse, sign);
            area = 0 + tempsum;
            return (num, perimeter, area);
        }

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
        public (int points, double perimeter, double area) TestPoint(GeoCoords coords, bool reverse, bool sign)
            => TestPoint(coords.Latitude, coords.Longitude, reverse, sign);

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
        public (int points, double perimeter, double area) TestEdge(double azi, double s, bool reverse, bool sign)
        {
            double perimeter, area = 0;
            if (_num == 0)
            {            // we don't have a starting point!
                perimeter = double.NaN;
                if (!_polyline)
                    area = double.NaN;
                return (0, perimeter, area);
            }
            var num = _num + 1;
            perimeter = _perimetersum + s;
            if (_polyline)
                return (num, perimeter, area);

            double tempsum = _areasum;
            int crossings = _crossings;
            {
                _earth.GenDirect(_lat1, _lon1, azi, false, s, _mask,
                                 out var lat, out var lon, out _, out _, out _, out _, out _, out var S12);
                tempsum += S12;
                crossings += TransitDirect(_lon1, lon);
                lon = AngNormalize(lon);
                _earth.GenInverse(lat, lon, _lat0, _lon0, _mask,
                                  out var s12, out _, out _, out _, out _, out _, out S12);
                perimeter += s12;
                tempsum += S12;
                crossings += Transit(lon, _lon0);
            }

            AreaReduce(ref tempsum, crossings, reverse, sign);
            area = 0 + tempsum;
            return (num, perimeter, area);
        }

        /// <summary>
        /// Gets a value representing the previous vertex added to the polygon or polyline.
        /// </summary>
        /// <remarks>
        /// If no points have been added, then <see cref="double.NaN"/>s are returned. Otherwise, <i>lon</i> will be in the range [−180°, 180°].
        /// </remarks>
        public (double lat, double lon) CurrentPoint => (_lat1, _lon1);
    }

    /// <summary>
    /// Polygon areas using <see cref="Geodesic"/>. This should be used if the flattening is small.
    /// </summary>
    public class PolygonArea : PolygonArea<Geodesic>
    {
        /// <summary>
        /// Constructor for <see cref="PolygonArea"/>.
        /// </summary>
        /// <param name="earth">the <see cref="Geodesic"/> object to use for geodesic calculations.</param>
        /// <param name="polyline">
        /// if <see langword="true"/> that treat the points as defining a polyline instead of a polygon (default = <see langword="false"/>).
        /// </param>
        public PolygonArea(Geodesic earth, bool polyline = false):base(earth, polyline) { }
    }

    /// <summary>
    /// Polygon areas using <see cref="GeodesicExact"/>.
    /// (But note that the implementation of areas in <see cref="GeodesicExact"/> uses a high order series and this is only accurate for modest flattenings.)
    /// </summary>
    public class PolygonAreaExact : PolygonArea<GeodesicExact>
    {
        /// <summary>
        /// Constructor for <see cref="PolygonAreaExact"/>.
        /// </summary>
        /// <param name="earth">the <see cref="GeodesicExact"/> object to use for geodesic calculations.</param>
        /// <param name="polyline">
        /// if <see langword="true"/> that treat the points as defining a polyline instead of a polygon (default = <see langword="false"/>).
        /// </param>
        public PolygonAreaExact(GeodesicExact earth, bool polyline = false) : base(earth, polyline) { }
    }

    /// <summary>
    /// Polygon areas using <see cref="Rhumb"/>.
    /// </summary>
    public class PolygonAreaRhumb : PolygonArea<Rhumb>
    {
        /// <summary>
        /// Constructor for <see cref="PolygonAreaRhumb"/>.
        /// </summary>
        /// <param name="earth">the <see cref="Rhumb"/> object to use for geodesic calculations.</param>
        /// <param name="polyline">
        /// if <see langword="true"/> that treat the points as defining a polyline instead of a polygon (default = <see langword="false"/>).
        /// </param>
        public PolygonAreaRhumb(Rhumb earth, bool polyline = false) : base(earth, polyline) { }
    }
}
