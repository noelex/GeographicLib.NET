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
    /// Gnomonic projection centered at an arbitrary position C on the ellipsoid.
    /// This projection is derived in Section 8 of
    /// <list type="bullet">
    /// <item>C. F. F. Karney, <a href="https://doi.org/10.1007/s00190-012-0578-z">
    /// Algorithms for geodesics</a>, J. Geodesy 87, 43–55 (2013);
    /// DOI: <a href="https://doi.org/10.1007/s00190-012-0578-z">10.1007/s00190-012-0578-z</a>;
    /// addenda: <a href="https://geographiclib.sourceforge.io/geod-addenda.html">geod-addenda.html</a>.</item>
    /// </list>
    /// The projection of <i>P</i> is defined as follows:
    /// compute the geodesic line from <i>C</i> to <i>P</i>; compute the reduced length <i>m12</i>,
    /// geodesic scale <i>M12</i>, and ρ = <i>m12</i>/<i>M12</i>; finally <i>x</i> = ρ sin <i>azi1</i>;
    /// <i>y</i> = ρ cos <i>azi1</i>, where <i>azi1</i> is the azimuth of the geodesic at <i>C</i>.
    /// The <see cref="Forward(double, double, double, double)"/> and 
    /// <see cref="Reverse(double, double, double, double)"/> methods also return the azimuth <i>azi</i>
    /// of the geodesic at <i>P</i> and reciprocal scale <i>rk</i> in the azimuthal direction.
    /// The scale in the radial direction if 1/<i>rk</i>^2.
    /// <para>
    /// For a sphere, ρ is reduces to a tan(<i>s12</i>/<i>a</i>),
    /// where <i>s12</i> is the length of the geodesic from <i>C</i> to <i>P</i>,
    /// and the gnomonic projection has the property that all geodesics appear as straight lines.
    /// For an ellipsoid, this property holds only for geodesics interesting the centers.
    /// However geodesic segments close to the center are approximately straight.
    /// </para>
    /// <para>
    /// Consider a geodesic segment of length <i>l</i>.
    /// Let <i>T</i> be the point on the geodesic (extended if necessary) closest to <i>C</i> the
    /// center of the projection and t be the distance <i>CT</i>.
    /// To lowest order, the maximum deviation (as a true distance) of the corresponding gnomonic
    /// line segment (i.e., with the same end points) from the geodesic is
    /// </para>
    /// <list type="table">
    /// <item>(<i>K</i>(<i>T</i>) - <i>K</i>(<i>C</i>)) <i>l</i>^2 <i>t</i> / 32.</item>
    /// </list>
    /// where <i>K</i> is the Gaussian curvature.
    /// <para>
    /// This result applies for any surface. For an ellipsoid of revolution,
    /// consider all geodesics whose end points are within a distance <i>r</i> of <i>C</i>.
    /// For a given <i>r</i>, the deviation is maximum when the latitude of <i>C</i> is 45°,
    /// when endpoints are a distance <i>r</i> away, and when their azimuths from the center
    /// are ± 45° or ± 135°. To lowest order in <i>r</i> and the flattening <i>f</i>,
    /// the deviation is <i>f</i> (<i>r</i>/2<i>a</i>)^3 <i>r</i>.
    /// </para>
    /// <para>
    /// The conversions all take place using a <see cref="Geodesic"/> object
    /// (by default <see cref="Geodesic.WGS84"/>).
    /// For more information on geodesics see
    /// <a href="https://geographiclib.sourceforge.io/C++/doc/geodesic.html">
    /// Geodesics on an ellipsoid of revolution</a>.
    /// </para>
    /// <b>Warning</b>
    /// <para>
    /// The definition of this projection for a sphere is standard.
    /// However, there is no standard for how it should be extended to an ellipsoid.
    /// The choices are:
    /// </para>
    /// <list type="bullet">
    /// <item>Declare that the projection is undefined for an ellipsoid.</item>
    /// <item>
    /// Project to a tangent plane from the center of the ellipsoid.
    /// This causes great ellipses to appear as straight lines in the projection;
    /// i.e., it generalizes the spherical great circle to a great ellipse.
    /// This was proposed by independently by Bowring and Williams in 1997.
    /// </item>
    /// <item>
    /// Project to the conformal sphere with the constant of integration chosen so that
    /// the values of the latitude match for the center point and perform a central
    /// projection onto the plane tangent to the conformal sphere at the center point.
    /// This causes normal sections through the center point to appear as straight lines
    /// in the projection; i.e., it generalizes the spherical great circle to a normal
    /// section. This was proposed by I. G. Letoval'tsev, Generalization of the gnomonic
    /// projection for a spheroid and the principal geodetic problems involved in the
    /// alignment of surface routes, Geodesy and Aerophotography (5), 271–274 (1963).
    /// </item>
    /// <item>
    /// The projection given here. This causes geodesics close to the center point to
    /// appear as straight lines in the projection; i.e., it generalizes the spherical
    /// great circle to a geodesic.
    /// </item>
    /// </list>
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
