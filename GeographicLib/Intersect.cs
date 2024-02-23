using System;
using System.Collections.Generic;
using System.Linq;
using static GeographicLib.Macros;
using static GeographicLib.MathEx;
using static System.Math;

namespace GeographicLib
{
    /// <summary>
    /// Geodesic intersections
    /// </summary>
    /// <remarks>
    /// Find the intersections of two geodesics <i>X</i> and <i>Y</i>. Four calling sequences are supported.
    /// <list type="bullet">
    /// <item>
    /// The geodesics are defined by a position (latitude and longitude) and an azimuth.
    /// In this case the <b>closest</b> intersection is found.
    /// </item>
    /// <item>
    /// The geodesics are defined by two endpoints.
    /// The intersection of the two segments is found.
    /// It they don't intersect, the the closest intersection is returned.
    /// </item>
    /// <item>
    /// The geodesics are defined as an intersection point, a single position and two azimuths.
    /// In this case, the next closest intersection is found.
    /// </item>
    /// <item>The geodesics are defined as in the first case and all intersection within a specified distance are returned.</item>
    /// </list>
    /// In all cases the position of the intersection is given by
    /// the signed displacements <i>x</i> and <i>y</i> along the geodesics from the starting point
    /// (the first point in the case of a geodesic segment).
    /// The closest itersection is defined as the one that minimizes the L1 distance, Intersect.Dist([<i>x</i>, <i>y</i>) = |<i>x</i>| + |<i>y</i>|.
    /// <para>
    /// The routines also optionally return a coincidence indicator <i>c</i>.
    /// This is typically 0. However if the geodesics lie on top of one another at the point of intersection,
    /// then <i>c</i> is set to +1, if they are parallel, and −1, if they are antiparallel.
    /// </para>
    /// This solution for intersections is described in
    /// <list type="bullet">
    /// <item>
    /// C. F. F. Karney,
    /// <a href="https://arxiv.org/abs/2308.00495">Geodesic intersections</a>,
    /// Technical Report, SRI International(2023).
    /// <a href="https://arxiv.org/abs/2308.00495">arxiv:2308.00495</a>
    /// </item>
    /// </list>
    /// It is based on the work of
    /// <list type="bullet">
    /// <item>
    /// S. Baseldga and J. C. Martinez-Llario,
    /// <a href="https://doi.org/10.1007/s11200-017-1020-z">Intersection and point-to-line solutions for geodesics on the ellipsoid</a>,
    /// Stud. Geophys. Geod. 62, 353–363 (2018);
    /// DOI: <a href="https://doi.org/10.1007/s11200-017-1020-z">10.1007/s11200-017-1020-z</a>.
    /// </item>
    /// </list>
    /// </remarks>
    public class Intersect
    {
        /// <summary>
        /// The minimum capabilities for <see cref="GeodesicLine"/> objects which are passed to this class.
        /// </summary>
        public const GeodesicFlags LineCaps =
            GeodesicFlags.Latitude |
            GeodesicFlags.Longitude |
            GeodesicFlags.Azimuth |
            GeodesicFlags.ReducedLength |
            GeodesicFlags.GeodesicScale |
            GeodesicFlags.DistanceIn;

        private const int numit_ = 100;

        private readonly Geodesic _geod;
        private readonly double
          _a, _f,                   // equatorial radius, flattening
          _R,                       // authalic radius
          _d,                       // pi*_R
          _eps,                     // criterion for intersection + coincidence
          _tol,                     // convergence for Newton in Solve1
          _delta,                   // for equality tests, safety margin for tiling
          _t1,                      // min distance between intersections
          _t2,                      // furthest dist to closest intersection
          _t3,                      // 1/2 furthest min dist to next intersection
          _t4,                      // capture radius for spherical sol in Solve0
          _t5,                      // longest shortest geodesic
          _d1,                      // tile spacing for Closest
          _d2,                      // tile spacing for Next
          _d3;                      // tile spacing for All
        private readonly SetComp _comp;

        private long _cnt0, _cnt1, _cnt2, _cnt3, _cnt4;

        private static readonly RankPoint s_defaultRankPoint = new RankPoint(Point.Zero);

        /// <summary>
        /// Initializes a new instance of the <see cref="AuxLatitude"/> class
        /// with the given <see cref="Geodesic"/> object. 
        /// </summary>
        /// <param name="geod">A <see cref="Geodesic"/> object. This sets the parameters a and f for the ellipsoid.</param>
        /// <remarks>
        /// This class has been validated for -1/4 ≤ <i>f</i> ≤ 1/5.
        /// It may give satisfactory results slightly outside this range;
        /// however sufficient far outside the range, some internal checks will fail and an exception thrown.
        /// <para>
        /// If |<i>f</i>| > 1/50, then the <see cref="Geodesic"/> object should be constructed with <i>exact</i> = <see langword="true"/>.
        /// </para>
        /// </remarks>
        public unsafe Intersect(Geodesic geod)
        {
            _geod = geod;
            _a = _geod.EquatorialRadius;
            _f = _geod.Flattening;
            _R = Sqrt(_geod.EllipsoidArea / (4 * PI));
            _d = _R * PI;       // Used to normalize intersection points
            _eps = 3 * DBL_EPSILON;
            _tol = _d * Pow(DBL_EPSILON, 3 / 4.0);
            _delta = _d * Pow(DBL_EPSILON, 1 / 5.0);
            _comp = new SetComp(_delta);
            _cnt0 = 0;
            _cnt1 = 0;
            _cnt2 = 0;
            _cnt3 = 0;
            _cnt4 = 0;

            _t1 = _t4 = _a * (1 - _f) * PI;
            _t2 = 2 * distpolar(90);
            _geod.Inverse(0, 0, 90, 0, out _t5); _t5 *= 2;
            if (_f > 0)
            {
                _t3 = distoblique();
                _t4 = _t1;
            }
            else
            {
                _t3 = _t5;
                _t4 = polarb();
                Swap(ref _t1, ref _t2);
            }
            _d1 = _t2 / 2;
            _d2 = 2 * _t3 / 3;
            _d3 = _t4 - _delta;
            if (!(_d1 < _d3 && _d2 < _d3 && _d2 < 2 * _t1))
                throw new GeographicException("Ellipsoid too eccentric for Closest");
        }

        /// <summary>
        /// Find the closest intersection point, with each geodesic specified by position and azimuth.
        /// </summary>
        /// <param name="latX">Latitude of starting point for geodesic <i>X</i> (degrees).</param>
        /// <param name="lonX">Longitude of starting point for geodesic <i>X</i> (degrees).</param>
        /// <param name="aziX">Azimuth of starting point for geodesic <i>X</i> (degrees).</param>
        /// <param name="latY">Latitude of starting point for geodesic <i>Y</i> (degrees).</param>
        /// <param name="lonY">Longitude of starting point for geodesic <i>Y</i> (degrees).</param>
        /// <param name="aziY">Azimuth of starting point for geodesic <i>Y</i> (degrees).</param>
        /// <param name="p0">Offset for the starting points (meters).</param>
        /// <param name="c">Coincidence indicator.</param>
        /// <returns>
        /// <i>p</i>, the intersection point closest to <i>p0</i>.
        /// </returns>
        /// <remarks>
        /// The returned intersection minimizes Intersect.Dist(<i>p</i>, <i>p0</i>).
        /// </remarks>
        public unsafe Point Closest(
                double latX, double lonX, double aziX,
                double latY, double lonY, double aziY,
                Point p0, out int c)
        {
            c = 0;
            fixed (int* pc = &c)
            {
                return Closest(
                            _geod.Line(latX, lonX, aziX),
                            _geod.Line(latY, lonY, aziY),
                            p0, pc);
            }
        }

        /// <summary>
        /// Find the closest intersection point, with each geodesic specified by position and azimuth.
        /// </summary>
        /// <param name="latX">Latitude of starting point for geodesic <i>X</i> (degrees).</param>
        /// <param name="lonX">Longitude of starting point for geodesic <i>X</i> (degrees).</param>
        /// <param name="aziX">Azimuth of starting point for geodesic <i>X</i> (degrees).</param>
        /// <param name="latY">Latitude of starting point for geodesic <i>Y</i> (degrees).</param>
        /// <param name="lonY">Longitude of starting point for geodesic <i>Y</i> (degrees).</param>
        /// <param name="aziY">Azimuth of starting point for geodesic <i>Y</i> (degrees).</param>
        /// <param name="c">Coincidence indicator.</param>
        /// <returns>
        /// <i>p</i>, the intersection point closest to <i>p0</i>.
        /// </returns>
        /// <remarks>
        /// The returned intersection minimizes Intersect.Dist(<i>p</i>, <i>p0</i>).
        /// </remarks>
        public unsafe Point Closest(
                double latX, double lonX, double aziX,
                double latY, double lonY, double aziY,
                out int c)
        {
            c = 0;
            fixed (int* pc = &c)
            {
                return Closest(
                        _geod.Line(latX, lonX, aziX),
                        _geod.Line(latY, lonY, aziY),
                        Point.Zero, pc);
            }
        }

        /// <summary>
        /// Find the closest intersection point, with each geodesic specified by position and azimuth.
        /// </summary>
        /// <param name="latX">Latitude of starting point for geodesic <i>X</i> (degrees).</param>
        /// <param name="lonX">Longitude of starting point for geodesic <i>X</i> (degrees).</param>
        /// <param name="aziX">Azimuth of starting point for geodesic <i>X</i> (degrees).</param>
        /// <param name="latY">Latitude of starting point for geodesic <i>Y</i> (degrees).</param>
        /// <param name="lonY">Longitude of starting point for geodesic <i>Y</i> (degrees).</param>
        /// <param name="aziY">Azimuth of starting point for geodesic <i>Y</i> (degrees).</param>
        /// <param name="p0">Offset for the starting points (meters).</param>
        /// <returns>
        /// <i>p</i>, the intersection point closest to <i>p0</i>.
        /// </returns>
        /// <remarks>
        /// The returned intersection minimizes Intersect.Dist(<i>p</i>, <i>p0</i>).
        /// </remarks>
        public unsafe Point Closest(
            double latX, double lonX, double aziX,
            double latY, double lonY, double aziY,
            Point p0)
        {
            return Closest(
                        _geod.Line(latX, lonX, aziX),
                        _geod.Line(latY, lonY, aziY),
                        p0, null);
        }

        /// <summary>
        /// Find the closest intersection point, with each geodesic specified by position and azimuth.
        /// </summary>
        /// <param name="lineX">Geodesic <i>X</i>.</param>
        /// <param name="lineY">Geodesic <i>Y</i></param>
        /// <param name="p0">Offset for the starting points (meters).</param>
        /// <param name="c">Coincidence indicator.</param>
        /// <returns>
        /// <i>p</i>, the intersection point closest to <i>p0</i>.
        /// </returns>
        /// <remarks>
        /// <paramref name="lineX"/> and <paramref name="lineY"/> should be created with
        /// minimum capabilities <see cref="LineCaps"/>.
        /// The methods for creating a <see cref="IGeodesicLine"/> include all these capabilities by default.
        /// </remarks>
        public unsafe Point Closest(
                IGeodesicLine lineX, IGeodesicLine lineY,
                Point p0, out int c)
        {
            c = 0;
            fixed (int* pc = &c)
            {
                return Closest(lineX, lineY, p0, pc);
            }
        }

        /// <summary>
        /// Find the closest intersection point, with each geodesic specified by position and azimuth.
        /// </summary>
        /// <param name="lineX">Geodesic <i>X</i>.</param>
        /// <param name="lineY">Geodesic <i>Y</i></param>
        /// <param name="p0">Offset for the starting points (meters).</param>
        /// <returns>
        /// <i>p</i>, the intersection point closest to <i>p0</i>.
        /// </returns>
        /// <remarks>
        /// <paramref name="lineX"/> and <paramref name="lineY"/> should be created with
        /// minimum capabilities <see cref="LineCaps"/>.
        /// The methods for creating a <see cref="IGeodesicLine"/> include all these capabilities by default.
        /// </remarks>
        public unsafe Point Closest(
                IGeodesicLine lineX, IGeodesicLine lineY,
                Point p0)
        {
            return Closest(lineX, lineY, p0, null);
        }

        /// <summary>
        /// Find the closest intersection point, with each geodesic specified by position and azimuth.
        /// </summary>
        /// <param name="lineX">Geodesic <i>X</i>.</param>
        /// <param name="lineY">Geodesic <i>Y</i></param>
        /// <param name="c">Coincidence indicator.</param>
        /// <returns>
        /// <i>p</i>, the intersection point closest to <i>p0</i>.
        /// </returns>
        /// <remarks>
        /// <paramref name="lineX"/> and <paramref name="lineY"/> should be created with
        /// minimum capabilities <see cref="Intersect.LineCaps"/>.
        /// The methods for creating a <see cref="IGeodesicLine"/> include all these capabilities by default.
        /// </remarks>
        public unsafe Point Closest(
                IGeodesicLine lineX, IGeodesicLine lineY,
                out int c)
        {
            c = 0;
            fixed (int* pc = &c)
            {
                return Closest(lineX, lineY, Point.Zero, pc);
            }
        }

        /// <summary>
        /// Find the intersection of two geodesic segments defined by their endpoints.
        /// </summary>
        /// <param name="latX1">Latitude of first point for segment <i>X</i> (degrees).</param>
        /// <param name="lonX1">Longitude of first point for segment <i>X</i> (degrees).</param>
        /// <param name="latX2">Latitude of second point for segment <i>X</i> (degrees).</param>
        /// <param name="lonX2">Longitude of second point for segment <i>X</i> (degrees).</param>
        /// <param name="latY1">Latitude of first point for segment <i>Y</i> (degrees).</param>
        /// <param name="lonY1">Longitude of first point for segment <i>Y</i> (degrees).</param>
        /// <param name="latY2">Latitude of second point for segment <i>Y</i> (degrees).</param>
        /// <param name="lonY2">Longitude of second point for segment <i>Y</i> (degrees).</param>
        /// <param name="segmode">An indicator equal to zero if the segments intersect.</param>
        /// <param name="c">Coincidence indicator.</param>
        /// <returns>
        /// <i>p</i>, the intersection point if the segments intersect, otherwise the intersection point closest to the midpoints of the two segments.
        /// </returns>
        /// <remarks>
        /// The results are only well defined if there's a <i>unique</i> shortest geodesic between the endpoints of the two segments.
        /// <para>
        /// <paramref name="segmode"/> codes up information about the closest intersection in the case where the segments intersect.
        /// Let <i>x12</i> be the length of the segment <i>X</i> and <i>x</i> = <i>p</i>.X, the position of the intersection on segment <i>X</i>. Define
        /// <list type="bullet">
        /// <item><i>kx</i> = −1, if <i>x</i> &lt; 0,</item>
        /// <item><i>kx</i> = 0, if 0 ≤ <i>x</i> ≤ <i>x12</i>,</item>
        /// <item><i>kx</i> = 1, if <i>x12</i> &lt; <i>x</i>.</item>
        /// </list>
        /// and similarly for segment <i>Y</i>. Then <paramref name="segmode"/> = 3 <i>kx</i> + <i>ky</i>.
        /// </para>
        /// </remarks>
        public unsafe Point Segment(
            double latX1, double lonX1,
            double latX2, double lonX2,
            double latY1, double lonY1,
            double latY2, double lonY2,
            out int segmode, out int c)
        {
            c = 0;
            fixed (int* pc = &c)
            {
                return Segment(
                    _geod.InverseLine(latX1, lonX1, latX2, lonX2, LineCaps),
                    _geod.InverseLine(latY1, lonY1, latY2, lonY2, LineCaps),
                    out segmode, pc);
            }
        }

        /// <summary>
        /// Find the intersection of two geodesic segments defined by their endpoints.
        /// </summary>
        /// <param name="latX1">Latitude of first point for segment <i>X</i> (degrees).</param>
        /// <param name="lonX1">Longitude of first point for segment <i>X</i> (degrees).</param>
        /// <param name="latX2">Latitude of second point for segment <i>X</i> (degrees).</param>
        /// <param name="lonX2">Longitude of second point for segment <i>X</i> (degrees).</param>
        /// <param name="latY1">Latitude of first point for segment <i>Y</i> (degrees).</param>
        /// <param name="lonY1">Longitude of first point for segment <i>Y</i> (degrees).</param>
        /// <param name="latY2">Latitude of second point for segment <i>Y</i> (degrees).</param>
        /// <param name="lonY2">Longitude of second point for segment <i>Y</i> (degrees).</param>
        /// <param name="segmode">An indicator equal to zero if the segments intersect.</param>
        /// <returns>
        /// <i>p</i>, the intersection point if the segments intersect, otherwise the intersection point closest to the midpoints of the two segments.
        /// </returns>
        /// <remarks>
        /// The results are only well defined if there's a <i>unique</i> shortest geodesic between the endpoints of the two segments.
        /// <para>
        /// <paramref name="segmode"/> codes up information about the closest intersection in the case where the segments intersect.
        /// Let <i>x12</i> be the length of the segment <i>X</i> and <i>x</i> = <i>p</i>.X, the position of the intersection on segment <i>X</i>. Define
        /// <list type="bullet">
        /// <item><i>kx</i> = −1, if <i>x</i> &lt; 0,</item>
        /// <item><i>kx</i> = 0, if 0 ≤ <i>x</i> ≤ <i>x12</i>,</item>
        /// <item><i>kx</i> = 1, if <i>x12</i> &lt; <i>x</i>.</item>
        /// </list>
        /// and similarly for segment <i>Y</i>. Then <paramref name="segmode"/> = 3 <i>kx</i> + <i>ky</i>.
        /// </para>
        /// </remarks>
        public unsafe Point Segment(
            double latX1, double lonX1,
            double latX2, double lonX2,
            double latY1, double lonY1,
            double latY2, double lonY2,
            out int segmode)
        {
            return Segment(
                _geod.InverseLine(latX1, lonX1, latX2, lonX2, LineCaps),
                _geod.InverseLine(latY1, lonY1, latY2, lonY2, LineCaps),
                out segmode, null);
        }

        /// <summary>
        /// Find the intersection of two geodesic segments each defined by a <see cref="GeodesicLine"/>.
        /// </summary>
        /// <param name="lineX">Segment <i>X</i>.</param>
        /// <param name="lineY">Segment <i>Y</i></param>
        /// <param name="segmode">An indicator equal to zero if the segments intersect.</param>
        /// <param name="c">Coincidence indicator.</param>
        /// <returns>
        /// <i>p</i>, the intersection point if the segments intersect, otherwise the intersection point closest to the midpoints of the two segments.
        /// </returns>
        /// <remarks>
        /// <paramref name="lineX"/> and <paramref name="lineY"/> must represent shortest geodesics,
        /// e.g., they can be created by <see cref="Geodesic.InverseLine(double, double, double, double, GeodesicFlags)"/>. 
        /// The results are only well defined if there's a <i>unique</i> shortest geodesic between the endpoints of the two segments.
        /// <para>
        /// <paramref name="lineX"/> and <paramref name="lineY"/> should be created with minimum capabilities <see cref="LineCaps"/>.
        /// The methods for creating a <see cref="GeodesicLine"/> include all these capabilities by default.
        /// </para>
        /// <para>
        /// <paramref name="segmode"/> codes up information about the closest intersection in the case where the segments intersect.
        /// Let <i>x12</i> be the length of the segment <i>X</i> and <i>x</i> = <i>p</i>.X, the position of the intersection on segment <i>X</i>. Define
        /// <list type="bullet">
        /// <item><i>kx</i> = −1, if <i>x</i> &lt; 0,</item>
        /// <item><i>kx</i> = 0, if 0 ≤ <i>x</i> ≤ <i>x12</i>,</item>
        /// <item><i>kx</i> = 1, if <i>x12</i> &lt; <i>x</i>.</item>
        /// </list>
        /// and similarly for segment <i>Y</i>. Then <paramref name="segmode"/> = 3 <i>kx</i> + <i>ky</i>.
        /// </para>
        /// </remarks>
        public unsafe Point Segment(
            GeodesicLine lineX, GeodesicLine lineY,
            out int segmode, out int c)
        {
            c = 0;
            fixed (int* pc = &c)
            {
                return Segment(
                    lineX, lineY,
                    out segmode, pc);
            }
        }

        /// <summary>
        /// Find the next closest intersection point to a given intersection, specified by position and two azimuths.
        /// </summary>
        /// <param name="latX">Latitude of starting points for geodesics <i>X</i> and <i>Y</i> (degrees).</param>
        /// <param name="lonX">Longitude of starting points for geodesics <i>X</i> and <i>Y</i> (degrees).</param>
        /// <param name="aziX">Azimuth at starting point for geodesic <i>X</i> (degrees).</param>
        /// <param name="aziY">Azimuth at starting point for geodesic <i>Y</i> (degrees).</param>
        /// <param name="c">Coincidence indicator.</param>
        /// <returns>
        /// <i>p</i>, the next closest intersection point.
        /// </returns>
        /// <remarks>
        /// The returned intersection minimizes Intersect.Dist(<i>p</i>) (excluding <i>p</i> = [0,0]).
        /// <para>
        /// Equidistant closest intersections are surprisingly common.
        /// If this may be a problem, use <see cref="All(double, double, double, double, double, double, double, out int[], Point)"/> with a sufficiently large <i>maxdist</i> to capture close intersections.
        /// </para>
        /// </remarks>
        public unsafe Point Next(
            double latX, double lonX,
            double aziX, double aziY,
            out int c)
        {
            c = 0;
            fixed (int* pc = &c)
            {
                return Next(
                    _geod.Line(latX, lonX, aziX, LineCaps),
                    _geod.Line(latX, lonX, aziY, LineCaps),
                    pc);
            }
        }

        /// <summary>
        /// Find the next closest intersection point to a given intersection, specified by position and two azimuths.
        /// </summary>
        /// <param name="latX">Latitude of starting points for geodesics <i>X</i> and <i>Y</i> (degrees).</param>
        /// <param name="lonX">Longitude of starting points for geodesics <i>X</i> and <i>Y</i> (degrees).</param>
        /// <param name="aziX">Azimuth at starting point for geodesic <i>X</i> (degrees).</param>
        /// <param name="aziY">Azimuth at starting point for geodesic <i>Y</i> (degrees).</param>
        /// <returns>
        /// <i>p</i>, the next closest intersection point.
        /// </returns>
        /// <remarks>
        /// The returned intersection minimizes Intersect.Dist(<i>p</i>) (excluding <i>p</i> = [0,0]).
        /// <para>
        /// Equidistant closest intersections are surprisingly common.
        /// If this may be a problem, use <see cref="All(double, double, double, double, double, double, double, Point)"/> with a sufficiently large <i>maxdist</i> to capture close intersections.
        /// </para>
        /// </remarks>
        public unsafe Point Next(
            double latX, double lonX,
            double aziX, double aziY)
        {
            return Next(
                _geod.Line(latX, lonX, aziX, LineCaps),
                _geod.Line(latX, lonX, aziY, LineCaps),
                null);
        }

        /// <summary>
        /// Find the next closest intersection point to a given intersection, with each geodesic specified a <see cref="GeodesicLine"/>.
        /// </summary>
        /// <param name="lineX">Geodesic <i>X</i>.</param>
        /// <param name="lineY">Geodesic <i>Y</i>.</param>
        /// <param name="c">Coincidence indicator.</param>
        /// <returns>
        /// <i>p</i>, the next closest intersection point.
        /// </returns>
        /// <remarks>
        /// <paramref name="lineX"/> and <paramref name="lineY"/> must both have the same starting point,
        /// i.e., the distance between [<paramref name="lineX"/>.Latitude, <paramref name="lineX"/>.Longitude]
        /// and [<paramref name="lineY"/>.Latitude, <paramref name="lineY"/>.Longitude] must be zero.
        /// <para>
        /// <paramref name="lineX"/> and <paramref name="lineY"/> should be created with minimum capabilities <see cref="LineCaps"/>.
        /// The methods for creating a <see cref="GeodesicLine"/> include all these capabilities by default.
        /// </para>
        /// <para>
        /// Equidistant closest intersections are surprisingly common.
        /// If this may be a problem, use <see cref="All(GeodesicLine, GeodesicLine, double, out int[], Point)"/> with a sufficiently large <i>maxdist</i> to capture close intersections.
        /// </para>
        /// </remarks>
        public unsafe Point Next(
            GeodesicLine lineX, GeodesicLine lineY, out int c)
        {
            c = 0;
            fixed (int* pc = &c)
            {
                return Next(
                    lineX, lineY,
                    pc);
            }
        }

        /// <summary>
        /// Find the next closest intersection point to a given intersection, with each geodesic specified a <see cref="GeodesicLine"/>.
        /// </summary>
        /// <param name="lineX">Geodesic <i>X</i>.</param>
        /// <param name="lineY">Geodesic <i>Y</i>.</param>
        /// <returns>
        /// <i>p</i>, the next closest intersection point.
        /// </returns>
        /// <remarks>
        /// <paramref name="lineX"/> and <paramref name="lineY"/> must both have the same starting point,
        /// i.e., the distance between [<paramref name="lineX"/>.Latitude, <paramref name="lineX"/>.Longitude]
        /// and [<paramref name="lineY"/>.Latitude, <paramref name="lineY"/>.Longitude] must be zero.
        /// <para>
        /// <paramref name="lineX"/> and <paramref name="lineY"/> should be created with minimum capabilities <see cref="LineCaps"/>.
        /// The methods for creating a <see cref="GeodesicLine"/> include all these capabilities by default.
        /// </para>
        /// <para>
        /// Equidistant closest intersections are surprisingly common.
        /// If this may be a problem, use <see cref="All(GeodesicLine, GeodesicLine, double, Point)"/> with a sufficiently large <i>maxdist</i> to capture close intersections.
        /// </para>
        /// </remarks>
        public unsafe Point Next(
            GeodesicLine lineX, GeodesicLine lineY)
        {
            return Next(
                lineX, lineY,
                null);
        }

        /// <summary>
        /// Find all intersections within a certain distance, with each geodesic specified by position and azimuth.
        /// </summary>
        /// <param name="latX">Latitude of starting point for geodesic <i>X</i> (degrees).</param>
        /// <param name="lonX">Longitude of starting point for geodesic <i>X</i> (degrees).</param>
        /// <param name="aziX">Azimuth of starting point for geodesic <i>X</i> (degrees).</param>
        /// <param name="latY">Latitude of starting point for geodesic <i>Y</i> (degrees).</param>
        /// <param name="lonY">Longitude of starting point for geodesic <i>Y</i> (degrees).</param>
        /// <param name="aziY">Azimuth of starting point for geodesic <i>Y</i> (degrees).</param>
        /// <param name="maxdist">The maximum distance for the returned intersections (meters).</param>
        /// <param name="c">Vector of coincidences.</param>
        /// <param name="p0">Offset for the starting points (meters), default = [0,0].</param>
        /// <returns>A vector for the intersections closest to <i>p0</i>.</returns>
        /// <remarks>
        /// Each intersection point satisfies Intersect.Dist(<i>p</i>, <i>p0</i>) ≤ <i>maxdist</i>.
        /// The vector of returned intersections is sorted on the distance from <i>p0</i>.
        /// </remarks>
        public Point[] All(
                  double latX, double lonX, double aziX,
                  double latY, double lonY, double aziY,
                  double maxdist, out int[] c, Point p0 = default)
        {
            return AllInternal(_geod.Line(latX, lonX, aziX, LineCaps),
                               _geod.Line(latY, lonY, aziY, LineCaps),
                               maxdist, p0, out c, true);
        }

        /// <summary>
        /// Find all intersections within a certain distance, with each geodesic specified by position and azimuth.
        /// Don't return vector of coincidences.
        /// </summary>
        /// <param name="latX">Latitude of starting point for geodesic <i>X</i> (degrees).</param>
        /// <param name="lonX">Longitude of starting point for geodesic <i>X</i> (degrees).</param>
        /// <param name="aziX">Azimuth of starting point for geodesic <i>X</i> (degrees).</param>
        /// <param name="latY">Latitude of starting point for geodesic <i>Y</i> (degrees).</param>
        /// <param name="lonY">Longitude of starting point for geodesic <i>Y</i> (degrees).</param>
        /// <param name="aziY">Azimuth of starting point for geodesic <i>Y</i> (degrees).</param>
        /// <param name="maxdist">The maximum distance for the returned intersections (meters).</param>
        /// <param name="p0">Offset for the starting points (meters), default = [0,0].</param>
        /// <returns>A vector for the intersections closest to <i>p0</i>.</returns>
        /// <remarks>
        /// Each intersection point satisfies Intersect.Dist(<i>p</i>, <i>p0</i>) ≤ <i>maxdist</i>.
        /// The vector of returned intersections is sorted on the distance from <i>p0</i>.
        /// </remarks>
        public Point[] All(
                  double latX, double lonX, double aziX,
                  double latY, double lonY, double aziY,
                  double maxdist, Point p0 = default)
        {
            return AllInternal(_geod.Line(latX, lonX, aziX, LineCaps),
                               _geod.Line(latY, lonY, aziY, LineCaps),
                               maxdist, p0, out _, false);
        }

        /// <summary>
        /// Find all intersections within a certain distance, with each geodesic specified by position and azimuth.
        /// </summary>
        /// <param name="lineX">Geodesic <i>X</i>.</param>
        /// <param name="lineY">Geodesic <i>Y</i>.</param>
        /// <param name="maxdist">The maximum distance for the returned intersections (meters).</param>
        /// <param name="c">Vector of coincidences.</param>
        /// <param name="p0">Offset for the starting points (meters), default = [0,0].</param>
        /// <returns>A vector for the intersections closest to <i>p0</i>.</returns>
        /// <remarks>
        /// Each intersection point satisfies Intersect.Dist(<i>p</i>, <i>p0</i>) ≤ <i>maxdist</i>.
        /// The vector of returned intersections is sorted on the distance from <i>p0</i>.
        /// <para>
        /// <paramref name="lineX"/> and <paramref name="lineY"/> should be created with minimum capabilities <see cref="LineCaps"/>.
        /// The methods for creating a <see cref="GeodesicLine"/> include all these capabilities by default.
        /// </para>
        /// </remarks>
        public Point[] All(
                  GeodesicLine lineX, GeodesicLine lineY,
                  double maxdist, out int[] c, Point p0 = default)
        {
            return AllInternal(lineX, lineY, maxdist, p0, out c, true);
        }

        /// <summary>
        /// Find all intersections within a certain distance, with each geodesic specified by position and azimuth.
        /// Don't return vector of coincidences.
        /// </summary>
        /// <param name="lineX">Geodesic <i>X</i>.</param>
        /// <param name="lineY">Geodesic <i>Y</i>.</param>
        /// <param name="maxdist">The maximum distance for the returned intersections (meters).</param>
        /// <param name="p0">Offset for the starting points (meters), default = [0,0].</param>
        /// <returns>A vector for the intersections closest to <i>p0</i>.</returns>
        /// <remarks>
        /// Each intersection point satisfies Intersect.Dist(<i>p</i>, <i>p0</i>) ≤ <i>maxdist</i>.
        /// The vector of returned intersections is sorted on the distance from <i>p0</i>.
        /// <para>
        /// <paramref name="lineX"/> and <paramref name="lineY"/> should be created with minimum capabilities <see cref="LineCaps"/>.
        /// The methods for creating a <see cref="GeodesicLine"/> include all these capabilities by default.
        /// </para>
        /// </remarks>
        public Point[] All(
                  GeodesicLine lineX, GeodesicLine lineY,
                  double maxdist, Point p0 = default)
        {
            return AllInternal(lineX, lineY, maxdist, p0, out _, false);
        }

        /// <summary>
        /// The L1 distance.
        /// </summary>
        /// <param name="p">The position along geodesics <i>X</i> and <i>Y</i>.</param>
        /// <param name="p0">The reference position, default = [0, 0].</param>
        /// <returns>
        /// The L1 distance of <i>p</i> from <i>p0</i>, i.e., |<i>px</i> − <i>p0x</i>| + |<i>py</i> − <i>p0y</i>|.
        /// </returns>
        public static double Dist(Point p, Point p0 = default)
        {
            return Abs(p.X - p0.X) + Abs(p.Y - p0.Y);
        }

        /// <summary>
        /// The <see cref="Geodesic"/> object used in the constructor.
        /// </summary>
        /// <remarks>
        /// This can be used to query <see cref="Geodesic.EquatorialRadius"/>,
        /// <see cref="Geodesic.Flattening"/>, <see cref="Geodesic.IsExact"/>,
        /// and <see cref="Geodesic.EllipsoidArea"/>.
        /// </remarks>
        public Geodesic Geodesic => _geod;

        /// <summary>
        /// The cumulative number of invocations of <b>h</b>.
        /// </summary>
        /// <remarks>
        /// This is a count of the number of times the spherical triangle needs to be solved.
        /// Each involves a call to Geodesic.GenInverse
        /// and this is a good metric for the overall cost. This counter is set to zero by the constructor.
        /// <para>
        /// The counter is a mutable variable and so is not thread safe.
        /// </para>
        /// </remarks>
        public long NumInverse => _cnt0;

        /// <summary>
        /// The cumulative number of invocations of <b>b</b>.
        /// </summary>
        /// <remarks>
        /// This is a count of the number of invocations of the basic algorithm, which is used by all the intersection methods.
        /// This counter is set to zero by the constructor.
        /// <para>
        /// The counter is a mutable variable and so is not thread safe.
        /// </para>
        /// </remarks>
        public long NumBasic => _cnt1;

        /// <summary>
        /// The number of times intersection point was changed in Intersect.Closest and Intersect.Next.
        /// </summary>
        /// <remarks>
        /// If this counter is incremented by just 1 in Intersect.Closest,
        /// then the initial result of the basic algorithm was eventually accepted.
        /// This counter is set to zero by the constructor.
        /// <para>
        /// This counter is also incremented by Intersect.Segment, which calls Intersect.Closest.
        /// </para>
        /// <para>
        /// The counter is a mutable variable and so is not thread safe.
        /// </para>
        /// </remarks>
        public long NumChange => _cnt2;

        /// <summary>
        /// The number of times a corner point is checked in Intersect.Segment.
        /// </summary>
        /// <remarks>
        /// This counter is set to zero by the constructor.
        /// <para>
        /// The counter is a mutable variable and so is not thread safe.
        /// </para>
        /// </remarks>
        public long NumCorner => _cnt3;

        /// <summary>
        /// The number of times a corner point is returned by Intersect.Segment.
        /// </summary>
        /// <remarks>
        /// This counter is set to zero by the constructor.
        /// <para>
        /// A conjecture is that a corner point never results in an intersection that
        /// overrides the intersection closest to the midpoints of the segments; i.e., <see cref="NumCorner"/> always returns 0.
        /// </para>
        /// <para>
        /// The counter is a mutable variable and so is not thread safe.
        /// </para>
        /// </remarks>
        public long NumOverride => _cnt4;

        private unsafe Point Closest(
                  IGeodesicLine lineX, IGeodesicLine lineY,
                  Point p0, int* c)
        {
            XPoint p = ClosestInt(lineX, lineY, new XPoint(p0));
            if (c != null) *c = p.C;
            return p.ToPoint();
        }

        private unsafe Point Segment(
            IGeodesicLine lineX, IGeodesicLine lineY,
            out int segmode, int* c)
        {

            XPoint p = SegmentInt(lineX, lineY, out segmode);
            if (c != null) *c = p.C;
            return p.ToPoint();
        }

        private unsafe Point Next(
            IGeodesicLine lineX, IGeodesicLine lineY,
            int* c)
        {
            XPoint p = NextInt(lineX, lineY);
            if (c != null) *c = p.C;
            return p.ToPoint();
        }

        private static double d1(double x, double y) => Abs(x) + Abs(y);

        /// <summary>
        /// The spherical solution
        /// </summary>
        /// <param name="lineX"></param>
        /// <param name="lineY"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        private XPoint Spherical(IGeodesicLine lineX, IGeodesicLine lineY,
                     XPoint p)
        {
            // threshold for coincident geodesics and intersections; this corresponds
            // to about 4.3 nm on WGS84.
            double latX, lonX, aziX, latY, lonY, aziY;
            lineX.Position(p.X, out latX, out lonX, out aziX);
            lineY.Position(p.Y, out latY, out lonY, out aziY);
            double z, aziXa, aziYa;
            _geod.Inverse(latX, lonX, latY, lonY, out z, out aziXa, out aziYa);
            double sinz = Sin(z / _R), cosz = Cos(z / _R);
            // X = interior angle at X, Y = exterior angle at Y
            double dX, dY, dXY,
              X = AngDiff(aziX, aziXa, out dX), Y = AngDiff(aziY, aziYa, out dY),
              XY = AngDiff(X, Y, out dXY);
            double s = CopySign(1, XY + (dXY + dY - dX)); // inverted triangle
                                                          // For z small, sinz -> z, cosz -> 1
                                                          // ( sinY*cosX*cosz - cosY*sinX) =
                                                          // (-sinX*cosY*cosz + cosX*sinY) -> sin(Y-X)
                                                          // for z = pi, sinz -> 0, cosz -> -1
                                                          // ( sinY*cosX*cosz - cosY*sinX) -> -sin(Y+X)
                                                          // (-sinX*cosY*cosz + cosX*sinY) ->  sin(Y+X)
            double sinX, cosX; SinCosde(s * X, s * dX, out sinX, out cosX);
            double sinY, cosY; SinCosde(s * Y, s * dY, out sinY, out cosY);
            double sX, sY;
            int c;
            if (z <= _eps * _R)
            {
                sX = sY = 0;              // Already at intersection
                                          // Determine whether lineX and lineY are parallel or antiparallel
                if (Abs(sinX - sinY) <= _eps && Abs(cosX - cosY) <= _eps)
                    c = 1;
                else if (Abs(sinX + sinY) <= _eps && Abs(cosX + cosY) <= _eps)
                    c = -1;
                else
                    c = 0;
            }
            else if (Abs(sinX) <= _eps && Abs(sinY) <= _eps)
            {
                c = cosX * cosY > 0 ? 1 : -1;
                // Coincident geodesics, place intersection at midpoint
                sX = cosX * z / 2; sY = -cosY * z / 2;
                // alt1: sX =  cosX * z; sY = 0;
                // alt2: sY = -cosY * z; sX = 0;
            }
            else
            {
                // General case.  [SKIP: Divide args by |sinz| to avoid possible
                // underflow in {sinX,sinY}*sinz; this is probably not necessary].
                // Definitely need to treat sinz < 0 (z > pi*R) correctly.  Without
                // this we have some convergence failures in Basic.
                sX = _R * Atan2(sinY * sinz, sinY * cosX * cosz - cosY * sinX);
                sY = _R * Atan2(sinX * sinz, -sinX * cosY * cosz + cosX * sinY);
                c = 0;
            }
            return new XPoint(sX, sY, c);
        }

        /// <summary>
        /// The basic algorithm
        /// </summary>
        /// <param name="lineX"></param>
        /// <param name="lineY"></param>
        /// <param name="p0"></param>
        /// <returns></returns>
        private XPoint Basic(IGeodesicLine lineX, IGeodesicLine lineY,
                 XPoint p0)
        {
            ++_cnt1;
            XPoint q = p0;
            for (int n = 0; n < numit_ || GEOGRAPHICLIB_PANIC; ++n)
            {
                ++_cnt0;
                XPoint dq = Spherical(lineX, lineY, q);
                q += dq;
                if (q.C != 0 || !(dq.Dist() > _tol)) break; // break if nan
            }
            return q;
        }

        /// <summary>
        /// The closest intersecton
        /// </summary>
        /// <param name="lineX"></param>
        /// <param name="lineY"></param>
        /// <param name="p0"></param>
        /// <returns></returns>
        private XPoint ClosestInt(IGeodesicLine lineX, IGeodesicLine lineY,
                  XPoint p0)
        {
            const int num = 5;
            Span<int> ix = stackalloc[] { 0, 1, -1, 0, 0 },
                      iy = stackalloc[] { 0, 0, 0, 1, -1 };
            Span<bool> skip = stackalloc bool[] { false, false, false, false, false };
            XPoint q = XPoint.Empty;                    // Best intersection so far
            for (int n = 0; n < num; ++n)
            {
                if (skip[n]) continue;
                XPoint qx = Basic(lineX, lineY, p0 + new XPoint(ix[n] * _d1, iy[n] * _d1));
                qx = fixcoincident(p0, qx);
                if (_comp.Equals(q, qx)) continue;
                if (qx.Dist(p0) < _t1) { q = qx; ++_cnt2; break; }
                if (n == 0 || qx.Dist(p0) < q.Dist(p0)) { q = qx; ++_cnt2; }
                for (int m = n + 1; m < num; ++m)
                    skip[m] = skip[m] ||
                      qx.Dist(p0 + new XPoint(ix[m] * _d1, iy[m] * _d1)) < 2 * _t1 - _d1 - _delta;
            }
            return q;
        }

        /// <summary>
        /// The next intersecton
        /// </summary>
        /// <param name="lineX"></param>
        /// <param name="lineY"></param>
        /// <returns></returns>
        private XPoint NextInt(IGeodesicLine lineX, IGeodesicLine lineY)
        {
            const int num = 8;
            Span<int> ix = stackalloc[] { -1, -1, 1, 1, -2, 0, 2, 0 },
                      iy = stackalloc[] { -1, 1, -1, 1, 0, 2, 0, -2 };
            Span<bool> skip = stackalloc bool[] { false, false, false, false, false, false, false, false };
            XPoint z = XPoint.Zero,                              // for excluding the origin
                   q = new XPoint(double.PositiveInfinity, 0);   // Best intersection so far
            for (int n = 0; n < num; ++n)
            {
                if (skip[n]) continue;
                XPoint qx = Basic(lineX, lineY, new XPoint(ix[n] * _d2, iy[n] * _d2));
                qx = fixcoincident(z, qx);
                bool zerop = _comp.Equals(z, qx);
                if (qx.C == 0 && zerop) continue;
                if (qx.C != 0 && zerop)
                {
                    for (int sgn = -1; sgn <= 1; sgn += 2)
                    {
                        double s = ConjugateDist(lineX, sgn * _d, false);
                        XPoint qa = new XPoint(s, qx.C * s, qx.C);
                        if (qa.Dist() < q.Dist()) { q = qa; ++_cnt2; }
                    }
                }
                else
                {
                    if (qx.Dist() < q.Dist()) { q = qx; ++_cnt2; }
                }
                for (int sgn = -1; sgn <= 1; ++sgn)
                {
                    // if qx.c == 0 only process sgn == 0
                    // if zerop skip sgn == 0
                    if ((qx.C == 0 && sgn != 0) || (zerop && sgn == 0)) continue;
                    XPoint qy = qx.C != 0 ? qx + new XPoint(sgn * _d2, qx.C * sgn * _d2) : qx;
                    for (int m = n + 1; m < num; ++m)
                        skip[m] = skip[m] ||
                          qy.Dist(new XPoint(ix[m] * _d2, iy[m] * _d2)) < 2 * _t1 - _d2 - _delta;
                }
            }
            return q;
        }

        /// <summary>
        /// Segment intersecton
        /// </summary>
        /// <param name="lineX"></param>
        /// <param name="lineY"></param>
        /// <param name="segmode"></param>
        /// <returns></returns>
        private XPoint SegmentInt(IGeodesicLine lineX, IGeodesicLine lineY,
                      out int segmode)
        {
            // The conjecture is that whenever two geodesic segments intersect, the
            // intersection is the one that is closest to the midpoints of segments.
            // If this is proven, set conjectureproved to true.
            const bool conjectureproved = false;
            double sx = lineX.Distance, sy = lineY.Distance;
            // p0 is center of [sx,sy] rectangle, q is intersection closest to p0
            XPoint p0 = new XPoint(sx / 2, sy / 2), q = ClosestInt(lineX, lineY, p0);
            q = fixsegment(sx, sy, q);
            segmode = segmentmode(sx, sy, q);
            // Are corners of [sx,sy] rectangle further from p0 than q?
            if (!conjectureproved && segmode != 0 && p0.Dist() >= p0.Dist(q))
            {
                int segmodex = 1;
                XPoint qx = XPoint.Empty;
                // Cycle through 4 corners of [sx,sy] rectangle
                for (int ix = 0; ix < 2 && segmodex != 0; ++ix)
                {
                    for (int iy = 0; iy < 2 && segmodex != 0; ++iy)
                    {
                        XPoint t = new XPoint(ix * sx, iy * sy); // corner point
                                                                 // Is corner outside next intersection exclusion circle?
                        if (q.Dist(t) >= 2 * _t1)
                        {
                            ++_cnt3;
                            qx = Basic(lineX, lineY, t);
                            // fixsegment is not needed because the coincidence line must just
                            // slice off a corner of the sx x sy rectangle.
                            qx = fixcoincident(t, qx);
                            // No need to check if equal to q, because result is only accepted
                            // if segmode != 0 && segmodex == 0.
                            segmodex = segmentmode(sx, sy, qx);
                        }
                    }
                }
                if (segmodex == 0) { ++_cnt4; segmode = 0; q = qx; }
            }
            return q;
        }

        /// <summary>
        /// All intersectons
        /// </summary>
        /// <param name="lineX"></param>
        /// <param name="lineY"></param>
        /// <param name="maxdist"></param>
        /// <param name="p0"></param>
        /// <returns></returns>
        private XPoint[] AllInt0(IGeodesicLine lineX, IGeodesicLine lineY,
           double maxdist, XPoint p0)
        {
            double maxdistx = maxdist + _delta;
            int m = (int)(Ceiling(maxdistx / _d3)), // process m x m set of tiles
              m2 = m * m + (m - 1) % 2,                // add center tile if m is even
              n = m - 1;                             // Range of i, j = [-n:2:n]
            double d3 = maxdistx / m;                    // d3 <= _d3
            var start = new List<XPoint>(m2);
            var skip = new bool[m2];
            int c0 = 0;
            start.Add(p0);
            for (int i = -n; i <= n; i += 2)
                for (int j = -n; j <= n; j += 2)
                {
                    if (!(i == 0 && j == 0))
                        start.Add(p0 + new XPoint(d3 * (i + j) / 2, d3 * (i - j) / 2));
                }

            HashSet<XPoint> r = new HashSet<XPoint>(_comp), // Intersections found
                            c = new HashSet<XPoint>(_comp); // Closest coincident intersections
            List<XPoint> added = new List<XPoint>();

            for (int k = 0; k < m2; ++k)
            {
                if (skip[k]) continue;
                XPoint q = Basic(lineX, lineY, start[k]);
                if (r.Contains(q)   // intersection already found
                                    // or it's on a line of coincident intersections already processed
                    || (c0 != 0 && c.Contains(fixcoincident(p0, q))))
                    continue;
                added.Clear();
                if (q.C != 0)
                {
                    // This value of q.c must be constitent with c0
                    // assert(c0 == 0 || c0 == q.c);
                    c0 = q.C;
                    // Process coincident intersections
                    q = fixcoincident(p0, q);
                    c.Add(q);
                    // Elimate all existing intersections on this line (which
                    // didn't set c0).
                    r.RemoveWhere(qp => _comp.Equals(fixcoincident(p0, qp, c0), q));
                    double s0 = q.X;
                    XPoint qc;
                    double m12, M12, M21;
                    lineX.GenPosition(false, s0,
                                      GeodesicFlags.ReducedLength |
                                      GeodesicFlags.GeodesicScale,
                                      out _, out _, out _, out _, out m12, out M12, out M21, out _);
                    // Compute line of conjugate points
                    for (int sgn = -1; sgn <= 1; sgn += 2)
                    {
                        double sa = 0;
                        do
                        {
                            sa = ConjugateDist(lineX, s0 + sa + sgn * _d, false, m12, M12, M21)
                              - s0;
                            qc = q + new XPoint(sa, c0 * sa);
                            added.Add(qc);
                            r.Add(qc);
                        } while (qc.Dist(p0) <= maxdistx);
                    }
                }
                added.Add(q);
                r.Add(q);
                foreach (var qp in added)
                {
                    for (int l = k + 1; l < m2; ++l)
                        skip[l] = skip[l] || qp.Dist(start[l]) < 2 * _t1 - d3 - _delta;
                }
            }
            // Trim intersections to maxdist
            r.RemoveWhere(qp => !(qp.Dist(p0) <= maxdist));
            var rankPoint = p0.ToPoint() == Point.Zero ? s_defaultRankPoint : new RankPoint(p0);
            return r.OrderBy(x => x, rankPoint).ToArray();
        }

        private Point[] AllInternal(IGeodesicLine lineX, IGeodesicLine lineY,
                double maxdist, Point p0,
                out int[] c, bool cp)
        {
            var v = AllInt0(lineX, lineY, Max(0, maxdist), new XPoint(p0));
            c = cp ? v.Select(x => x.C).ToArray() : Array.Empty<int>();
            return v.Select(x => x.ToPoint()).ToArray();
        }

        private unsafe double distpolar(double lat1, double* lat2 = null)
        {
            IGeodesicLine line = _geod.Line(lat1, 0, 0,
                                   GeodesicFlags.ReducedLength |
                                   GeodesicFlags.GeodesicScale |
                                   GeodesicFlags.DistanceIn);
            double s = ConjugateDist(line, (1 + _f / 2) * _a * PI / 2, true);
            if (lat2 != null)
            {
                line.GenPosition(false, s, GeodesicFlags.Latitude,
                                 out *lat2, out _, out _, out _, out _, out _, out _, out _);
            }
            return s;
        }

        private unsafe double polarb(double* lata = null, double* latb = null)
        {
            if (_f == 0)
            {
                if (lata != null) *lata = 64;
                if (latb != null) *latb = 90 - 64;
                return _d;
            }
            double
              lat0 = 63, s0 = distpolar(lat0),
              lat1 = 65, s1 = distpolar(lat1),
              lat2 = 64, s2 = distpolar(lat2),
              latx = lat2, sx = s2;
            // Solve for ds(lat)/dlat = 0 with a quadratic fit
            for (int i = 0; i < 10; ++i)
            {
                double den = (lat1 - lat0) * s2 + (lat0 - lat2) * s1 + (lat2 - lat1) * s0;
                if (!(den < 0 || den > 0)) break; // Break if nan
                double latn = ((lat1 - lat0) * (lat1 + lat0) * s2 + (lat0 - lat2) * (lat0 + lat2) * s1 +
                             (lat2 - lat1) * (lat2 + lat1) * s0) / (2 * den);
                lat0 = lat1; s0 = s1;
                lat1 = lat2; s1 = s2;
                lat2 = latn; s2 = distpolar(lat2);
                if (_f < 0 ? (s2 < sx) : (s2 > sx))
                {
                    sx = s2;
                    latx = lat2;
                }
            }
            if (lata != null) *lata = latx;
            if (latb != null) distpolar(latx, latb);
            return 2 * sx;
        }

        /// <summary>
        ///  Find {semi-,}conjugate point which is close to s3.  Optional m12, M12,
        /// M21 use {semi-,}conjugacy relative to point 2
        /// </summary>
        /// <param name="line"></param>
        /// <param name="s3"></param>
        /// <param name="semi"></param>
        /// <param name="m12"></param>
        /// <param name="M12"></param>
        /// <param name="M21"></param>
        /// <returns></returns>
        private double ConjugateDist(IGeodesicLine line, double s3, bool semi,
                                  double m12 = 0, double M12 = 1,
                                  double M21 = 1)
        {
            // semi = false: solve for m23 = 0 using dm23/ds3 = M32
            // semi = true : solve for M23 = 0 using dM23/ds3 = - (1 - M23*M32)/m23
            // Here 2 is point with given m12, M12, M21 and default values s.t. point 2
            // = point 1.
            double s = s3;
            for (int i = 0; i < 100; ++i)
            {
                double m13, M13, M31;
                line.GenPosition(false, s,
                                 GeodesicFlags.ReducedLength |
                                 GeodesicFlags.GeodesicScale,
                                 out _, out _, out _, out _, out m13, out M13, out M31, out _);
                double
                  // See "Algorithms for geodesics", eqs. 31, 32, 33.
                  m23 = m13 * M12 - m12 * M13,
                  // when m12 -> eps, (1 - M12 * M21) -> eps^2, I suppose.
                  M23 = M13 * M21 + (m12 == 0 ? 0 : (1 - M12 * M21) * m13 / m12),
                  M32 = M31 * M12 + (m13 == 0 ? 0 : (1 - M13 * M31) * m12 / m13);
                double ds = semi ? m23 * M23 / (1 - M23 * M32) : -m23 / M32;
                s = s + ds;
                if (!(Abs(ds) > _tol)) break;
            }
            return s;
        }

        private unsafe double conjdist(double azi, double* ds = null,
                        double* sp = null, double* sm = null)
        {
            IGeodesicLine line = _geod.Line(0, 0, azi, LineCaps);
            double s = ConjugateDist(line, _d, false);
            if (ds != null)
            {
                XPoint p = Basic(line, line, new XPoint(s / 2, -3 * s / 2));
                if (sp != null) *sp = p.X;
                if (sm != null) *sm = p.Y;
                *ds = p.Dist() - 2 * s;
            }
            return s;
        }

        private unsafe double distoblique(double* azi = null, double* sp = null, double* sm = null)
        {
            if (_f == 0)
            {
                if (azi != null) *azi = 45;
                if (sp != null) *sp = 0.5;
                if (sm != null) *sm = -1.5;
                return _d;
            }
            double sa, sb,
              azi0 = 46, ds0, s0 = conjdist(azi0, &ds0, &sa, &sb),
              azi1 = 44, ds1, s1 = conjdist(azi1, &ds1, &sa, &sb),
              azix = azi1, dsx = Abs(ds1), sx = s1, sax = sa, sbx = sb;
            // find ds(azi) = 0 by secant method
            // (void)s0;
            for (int i = 0; i < 10 && ds1 != ds0; ++i)
            {
                double azin = (azi0 * ds1 - azi1 * ds0) / (ds1 - ds0);
                azi0 = azi1; s0 = s1; ds0 = ds1;
                azi1 = azin; s1 = conjdist(azi1, &ds1, &sa, &sb);
                if (Abs(ds1) < dsx)
                {
                    azix = azi1; sx = s1; dsx = Abs(ds1);
                    sax = sa; sbx = sb;
                    if (ds1 == 0) break;
                }
            }
            if (azi != null) *azi = azix;
            if (sp != null) *sp = sax;
            if (sm != null) *sm = sbx;
            return sx;
        }

        /// <summary>
        /// p is intersection point on coincident lines orientation = c; p0 is
        /// origin point.  Change p to center point wrt p0, i.e, abs((p-p0)_x) =
        /// abs((p-p0)_y)
        /// </summary>
        /// <param name="p0"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        static XPoint fixcoincident(XPoint p0, XPoint p)
        {
            return fixcoincident(p0, p, p.C);
        }

        static XPoint fixcoincident(XPoint p0, XPoint p, int c)
        {
            if (c == 0) return p;
            // eqs : [p0x-p1x = -c*(p0y-p1y), p1x = px+s, p1y = py+c*s]$
            // sol : solve(eqs,[s,p1x,p1y]);
            // =>
            // sol:[ s = ((p0x+c*p0y) - (px+c*py))/2,
            //       p1x = px +     ((p0x+c*p0y) - (px+c*py))/2,
            //       p1y = py + c * ((p0x+c*p0y) - (px+c*py))/2
            // ];
            double s = ((p0.X + c * p0.Y) - (p.X + c * p.Y)) / 2;
            return p + new XPoint(s, c * s);
        }

        static XPoint fixsegment(double sx, double sy, XPoint p)
        {
            if (p.C == 0) return p;
            // eq0: [p1x = px+s, p1y = py+f*s]$
            // solx0:linsolve(cons(p1x=0 ,eq0),[s,p1x,p1y]);
            // solx1:linsolve(cons(p1x=sx,eq0),[s,p1x,p1y]);
            // soly0:linsolve(cons(p1y=0 ,eq0),[s,p1x,p1y]);
            // soly1:linsolve(cons(p1y=sy,eq0),[s,p1x,p1y]);
            // solx0:[s = -px      ,p1x = 0 ,p1y = py-f*px     ];
            // solx1:[s = sx-px    ,p1x = sx,p1y = py-f*(px-sx)];
            // soly0:[s = -f*py    ,p1x = px-f*py     ,p1y = 0 ];
            // soly1:[s = f*(sy-py),p1x = px-f*(py-sy),p1y = sy];
            double
              pya = p.Y - p.C * p.X, sa = -p.X,  // pxa = 0
              pyb = p.Y - p.C * (p.X - sx), sb = sx - p.X,  // pxb = sx
              pxc = p.X - p.C * p.Y, sc = p.C * -p.Y,  // pyc = 0
              pxd = p.X - p.C * (p.Y - sy), sd = p.C * (sy - p.Y); // pyd = sy
            bool
              ga = 0 <= pya && pya <= sy,
              gb = 0 <= pyb && pyb <= sy,
              gc = 0 <= pxc && pxc <= sx,
              gd = 0 <= pxd && pxd <= sx;
            double s;
            // Test opposite sides of the rectangle first
            if (ga && gb) s = (sa + sb) / 2;
            else if (gc && gd) s = (sc + sd) / 2;
            else if (ga && gc) s = (sa + sc) / 2;
            else if (ga && gd) s = (sa + sd) / 2;
            else if (gb && gc) s = (sb + sc) / 2;
            else if (gb && gd) s = (sb + sd) / 2;
            else
            {
                // Intersection not within segments; place intersection in smallest gap.
                if (p.C > 0)
                {
                    // distance from p to corner p0 is abs( (px - py) - (p0x - p0y) )
                    // consider corners p0 = [0, sy] and p0 = [sx, 0]
                    if (Abs((p.X - p.Y) + sy) < Abs((p.X - p.Y) - sx))
                        s = (sy - (p.X + p.Y)) / 2;
                    else
                        s = (sx - (p.X + p.Y)) / 2;
                }
                else
                {
                    // distance from p to corner p0 is abs( (px + p.y) - (p0x + p0y) )
                    // consider corners p0 = [0, 0] and p0 = [sx, sy]
                    if (Abs(p.X + p.Y) < Abs((p.X + p.Y) - (sx + sy)))
                        s = (0 - (p.X - p.Y)) / 2;
                    else
                        s = ((sx - sy) - (p.X - p.Y)) / 2;
                }
            }
            return p + new XPoint(s, p.C * s);
        }

        private static int segmentmode(double sx, double sy, XPoint p)
        {
            return (p.X < 0 ? -1 : p.X <= sx ? 0 : 1) * 3
              + (p.Y < 0 ? -1 : p.Y <= sy ? 0 : 1);
        }

        private readonly struct XPoint
        {
            public readonly double X, Y;
            public readonly int C;

            public static XPoint Empty = new XPoint(double.NaN, double.NaN, 0);

            public static XPoint Zero = new XPoint(0, 0, 0);

            public XPoint(double x, double y, int c = 0)
            {
                X = x;
                Y = y;
                C = c;
            }

            public XPoint(Point p) : this(p.X, p.Y, 0)
            {

            }

            public static XPoint operator +(XPoint a, XPoint b)
            {
                var x = a.X + b.X;
                var y = a.Y + b.Y;
                var c = b.C == 0 ? a.C : b.C;
                return new XPoint(x, y, c);
            }

            public double Dist() => d1(X, Y);

            public double Dist(XPoint p) => d1(X - p.X, Y - p.Y);

            public Point ToPoint() => new Point(X, Y);
        }

        private class SetComp : IEqualityComparer<XPoint>
        {
            private readonly double _delta;

            public SetComp(double delta) => _delta = delta;

            public bool Equals(XPoint x, XPoint y)
            {
                return d1(x.X - y.X, x.Y - y.Y) <= _delta;
            }

            public int GetHashCode(XPoint obj)
            {
                long x = ((Bit64)obj.X).Int64,
                     y = ((Bit64)obj.Y).Int64;

                unchecked
                {
                    return ((int)(x & 0xffffffff)) ^ ((int)(x >> 32)) ^
                           ((int)(y & 0xffffffff)) ^ ((int)(y >> 32)) ^
                           obj.C;
                }
            }
        }

        private class RankPoint : IComparer<XPoint>
        {
            private readonly double _x, _y;

            public RankPoint(double x, double y)
            {
                (_x, _y) = (x, y);
            }

            public RankPoint(Point p) : this(p.X, p.Y) { }

            public RankPoint(XPoint p) : this(p.X, p.Y) { }

            public int Compare(XPoint p, XPoint q)
            {
                double dp = d1(p.X - _x, p.Y - _y),
                       dq = d1(q.X - _x, q.Y - _y);

                if (dp != dq)
                {
                    return dp.CompareTo(dq);
                }
                else if (p.X != q.X)
                {
                    return p.X.CompareTo(q.X);
                }
                else
                {
                    return p.Y.CompareTo(q.Y);
                }
            }
        }
    }
}
