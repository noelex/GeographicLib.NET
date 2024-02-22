using System;
using System.Collections.Generic;
using static GeographicLib.MathEx;
using static System.Math;

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
    /// addenda: <a href="https://geographiclib.sourceforge.io/geod-addenda.html">geod-addenda.html</a>.
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
    public class PolygonArea<T> : IPolygonArea where T : IGeodesicLike
    {
        private readonly T _earth;

        private readonly double _area0;
        private readonly bool _polyline;
        private readonly GeodesicFlags _mask;

        private int _num, _crossings;
        private double _lat0, _lon0, _lat1, _lon1;
        private Accumulator _areasum, _perimetersum;

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

        /// <inheritdoc/>
        public int Count => _num;

        /// <inheritdoc/>
        public bool IsPolyline => _polyline;

        private static int Transit(double lon1, double lon2)
        {
            // Return 1 or -1 if crossing prime meridian in east or west direction.
            // Otherwise return zero.  longitude = +/-0 considered to be positive.
            // This is (should be?) compatible with transitdirect which computes
            // exactly the parity of
            //   int(floor((lon1 + lon12) / 360)) - int(floor(lon1 / 360)))
            var lon12 = AngDiff(lon1, lon2);
            lon1 = AngNormalize(lon1);
            lon2 = AngNormalize(lon2);
            // N.B. lon12 == 0 gives cross = 0
            return
                // edge case lon1 = 180, lon2 = 360->0, lon12 = 180 to give 1
                lon12 > 0 && ((lon1 < 0 && lon2 >= 0) ||
                              // lon12 > 0 && lon1 > 0 && lon2 == 0 implies lon1 == 180
                              (lon1 > 0 && lon2 == 0)) ? 1 :
                // non edge case lon1 = -180, lon2 = -360->-0, lon12 = -180
                (lon12 < 0 && lon1 >= 0 && lon2 < 0 ? -1 : 0);
            // This was the old method (treating +/- 0 as negative).  However, with the
            // new scheme for handling longitude differences this fails on:
            // lon1 = -180, lon2 = -360->-0, lon12 = -180 gives 0 not -1.
            //    return
            //      lon1 <= 0 && lon2 > 0 && lon12 > 0 ? 1 :
            //      (lon2 <= 0 && lon1 > 0 && lon12 < 0 ? -1 : 0);
        }

        // an alternate version of transit to deal with longitudes in the direct
        // problem.
        private static int TransitDirect(double lon1, double lon2)
        {
            // Compute exactly the parity of
            //   int(floor(lon2 / 360)) - int(floor(lon1 / 360))

            // C++ C remainder -> [-360,360]
            // Java % -> (-720, 720) switch to IEEEremainder?
            // JS % -> (-720, 720)
            // Python fmod -> (-720, 720)
            // Fortran, Octave skip
            // If mod function gives result in [-360, 360]
            // [0, 360) -> 0; [-360, 0) or 360 -> 1
            // If mod function gives result in (-720, 720)
            // [0, 360) or [-inf, -360) -> 0; [-360, 0) or [360, inf) -> 1
            lon1 = IEEERemainder(lon1, 2 * TD);
            lon2 = IEEERemainder(lon2, 2 * TD);
            return ((lon2 >= 0 && lon2 < TD ? 0 : 1) -
                     (lon1 >= 0 && lon1 < TD ? 0 : 1));
        }

        private void Remainder(ref Accumulator a) => a.Remainder(_area0);

        private void Remainder(ref double a) => a = IEEERemainder(a, _area0);

        private void AreaReduce(ref Accumulator area, int crossings, bool reverse, bool sign)
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

        /// <inheritdoc/>
        public void Clear()
        {
            _num = 0;
            _crossings = 0;
            _areasum = 0;
            _perimetersum = 0;
            _lat0 = _lon0 = _lat1 = _lon1 = double.NaN;
        }

        /// <inheritdoc/>
        public void AddPoint(double lat, double lon)
        {
            lat = LatFix(lat);
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

        /// <inheritdoc/>
        public void AddPoint(GeoCoords coords) => AddPoint(coords.Latitude, coords.Longitude);

        /// <inheritdoc/>
        public void AddPoints(IEnumerable<GeoCoords> coords)
        {
            foreach (var coord in coords)
            {
                AddPoint(coord);
            }
        }

        /// <inheritdoc/>
        public void AddEdge(double azi, double s)
        {
            if (_num != 0)
            {                 // Do nothing if _num is zero
                _earth.GenDirect(_lat1, _lon1, azi, false, s, _mask,
                                 out var lat, out var lon, out _, out _, out _, out _, out _, out var S12);
                _perimetersum += s;
                if (!_polyline)
                {
                    _areasum += S12;
                    _crossings += TransitDirect(_lon1, lon);
                }
                _lat1 = lat; _lon1 = lon;
                ++_num;
            }
        }

        /// <inheritdoc/>
        public (int points, double perimeter, double area) Compute(bool reverse, bool sign)
        {
            double perimeter, area = 0;
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
            area = 0d + tempsum;
            return (_num, perimeter, area);
        }

        /// <inheritdoc/>
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
            area = 0d + tempsum;
            return (num, perimeter, area);
        }

        /// <inheritdoc/>
        public (int points, double perimeter, double area) TestPoint(GeoCoords coords, bool reverse, bool sign)
            => TestPoint(coords.Latitude, coords.Longitude, reverse, sign);

        /// <inheritdoc/>
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

        /// <inheritdoc/>
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
        public PolygonArea(Geodesic earth, bool polyline = false) : base(earth, polyline) { }
    }

    /// <summary>
    /// Polygon areas using <see cref="GeodesicExact"/>.
    /// (But note that the implementation of areas in <see cref="GeodesicExact"/> uses a high order series and this is only accurate for modest flattenings.)
    /// </summary>
    [Obsolete("Use PolygonArea with a Geodesic object specified with exact = true.")]
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
