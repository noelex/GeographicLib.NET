using System;
using System.Collections.Generic;
using System.Text;
using GeographicLib;

using static System.Math;
using static GeographicLib.MathEx;

namespace GeographicLib.Projections
{
    /// <summary>
    /// Cassini-Soldner projection.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Cassini-Soldner projection centered at an arbitrary position, <i>lat0</i>, <i>lon0</i>, on the ellipsoid.
    /// This projection is a transverse cylindrical equidistant projection.
    /// The projection from (<i>lat</i>, <i>lon</i>) to easting and northing (<i>x</i>, <i>y</i>) is defined by geodesics as follows.
    /// Go north along a geodesic a distance <i>y</i> from the central point; then turn clockwise 90° and go a distance <i>x</i> along a geodesic.
    /// (Although the initial heading is north, this changes to south if the pole is crossed.) This procedure uniquely defines the reverse projection.
    /// The forward projection is constructed as follows. Find the point (<i>lat1</i>, <i>lon1</i>) on the meridian closest to (<i>lat</i>, <i>lon</i>).
    /// Here we consider the full meridian so that lon1 may be either <i>lon0</i> or <i>lon0</i> + 180°.
    /// x is the geodesic distance from (<i>lat1</i>, <i>lon1</i>) to (<i>lat</i>, <i>lon</i>), appropriately signed according to which side of the 
    /// central meridian (<i>lat</i>, <i>lon</i>) lies. y is the shortest distance along the meridian from (<i>lat0</i>, <i>lon0</i>) to 
    /// (<i>lat1</i>, <i>lon1</i>), again, appropriately signed according to the initial heading. 
    /// [Note that, in the case of prolate ellipsoids, the shortest meridional path from (<i>lat0</i>, <i>lon0</i>) to (<i>lat1</i>, <i>lon1</i>) 
    /// may not be the shortest path.]
    /// This procedure uniquely defines the forward projection except for a small class of points for which there may be two equally short routes for
    /// either leg of the path.
    /// </para>
    /// <para>
    /// Because of the properties of geodesics, the (<i>x</i>, <i>y</i>) grid is orthogonal.
    /// The scale in the easting direction is unity. The scale, <i>k</i>, in the northing direction is unity on the central meridian and increases away 
    /// from the central meridian.
    /// The projection routines return azi, the true bearing of the easting direction, and <i>rk</i> = 1/<i>k</i>, the reciprocal of the scale in the
    /// northing direction.
    /// </para>
    /// <para>
    /// The conversions all take place using a <see cref="Geodesic"/> object (by default <see cref="Geodesic.WGS84"/>).
    /// For more information on geodesics see <a href="https://geographiclib.sourceforge.io/html/geodesic.html">Geodesics on an ellipsoid of revolution</a>. 
    /// The determination of (<i>lat1</i>, <i>lon1</i>) in the forward projection is by solving the inverse geodesic problem for (<i>lat</i>, <i>lon</i>)
    /// and its twin obtained by reflection in the meridional plane. 
    /// The scale is found by determining where two neighboring geodesics intersecting the central meridian at <i>lat1</i> and <i>lat1</i> + <i>dlat1</i>
    /// intersect and taking the ratio of the reduced lengths for the two geodesics between that point and, respectively, (<i>lat1</i>, <i>lon1</i>) and (<i>lat</i>, <i>lon</i>).
    /// </para>
    /// </remarks>
    public class CassiniSoldner : IEllipsoid
    {
        private const uint maxit_ = 10;

        private readonly IGeodesic _earth;

        private IGeodesicLine _meridian;
        private double _sbet0, _cbet0;

        /// <summary>
        /// Initialize a new <see cref="CassiniSoldner"/> instance with <see cref="Geodesic.WGS84"/>.
        /// </summary>
        public CassiniSoldner() : this(Geodesic.WGS84) { }

        /// <summary>
        /// Initialize a new <see cref="CassiniSoldner"/> instance with specified <see cref="IGeodesic"/> instance.
        /// </summary>
        /// <param name="earth">the <see cref="IGeodesic"/> object to use for geodesic calculations.</param>
        public CassiniSoldner(IGeodesic earth) : this(0, 0, earth) { }

        /// <summary>
        /// Initialize a new <see cref="CassiniSoldner"/> instance with specified center point and <see cref="Geodesic"/> instance.
        /// </summary>
        /// <param name="earth">the <see cref="Geodesic"/> object to use for geodesic calculations.</param>
        /// <param name="lat0">latitude of center point of projection (degrees).</param>
        /// <param name="lon0">longitude of center point of projection (degrees).</param>
        public CassiniSoldner(double lat0, double lon0, IGeodesic earth)
        {
            _earth = earth;
            Reset(lat0, lon0);
        }

        /// <summary>
        /// Gets a value representing the equatorial radius (<i>a</i>) of the ellipsoid.
        /// </summary>
        public double EquatorialRadius => _earth.EquatorialRadius;

        /// <summary>
        /// Gets a value representing the flatterning (<i>f</i>) of the ellipsoid.
        /// </summary>
        public double Flattening => _earth.Flattening;

        /// <summary>
        /// Gets a value representing the latitude of origin (degrees).
        /// </summary>
        public double LatitudeOrigin => (_meridian?.Latitude).GetValueOrDefault();

        /// <summary>
        /// Gets a value representing the longitude of origin (degrees).
        /// </summary>
        public double LongitudeOrigin => (_meridian?.Latitude).GetValueOrDefault();

        /// <summary>
        /// Set the central point of the projection
        /// </summary>
        /// <param name="lat0">latitude of center point of projection (degrees).</param>
        /// <param name="lon0">longitude of center point of projection (degrees).</param>
        /// <remarks>
        /// <paramref name="lat0"/> should be in the range [−90°, 90°].
        /// </remarks>
        public void Reset(double lat0, double lon0)
        {
            _meridian = _earth.Line(lat0, lon0, 0,
                        GeodesicFlags.Latitude | GeodesicFlags.Longitude |
                        GeodesicFlags.Distance | GeodesicFlags.DistanceIn |
                        GeodesicFlags.Azimuth);
            var f = _earth.Flattening;
            SinCosd(LatitudeOrigin, out _sbet0, out _cbet0);
            _sbet0 *= (1 - f);
            Norm(ref _sbet0, ref _cbet0);
        }

        /// <summary>
        /// Forward projection, from geographic to Cassini-Soldner.
        /// </summary>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <param name="azi">azimuth of easting direction at point (degrees).</param>
        /// <param name="rk">reciprocal of azimuthal northing scale at point.</param>
        /// <returns>
        /// <i>x</i>, easting of point and <i>y</i>, northing of point, in meters.
        /// </returns>
        /// <remarks>
        /// <paramref name="lat"/> should be in the range [−90°, 90°].
        /// A call to <see cref="Forward(double, double, out double, out double)"/> 
        /// followed by a call to <see cref="Reverse(double, double, out double, out double)"/> 
        /// will return the original (<paramref name="lat"/>, <paramref name="lon"/>) (to within roundoff).
        /// The routine does nothing if the origin has not been set.
        /// </remarks>
        public (double x, double y) Forward(double lat, double lon, out double azi, out double rk)
        {
            var dlon = AngDiff(LongitudeOrigin, lon);
            var sig12 = 
                _earth.Inverse(lat, -Abs(dlon), lat, Abs(dlon), out var s12, out var azi1, out var azi2);
            sig12 *= 0.5;
            s12 *= 0.5;
            if (s12 == 0)
            {
                var da = AngDiff(azi1, azi2) / 2;
                if (Abs(dlon) <= 90)
                {
                    azi1 = 90 - da;
                    azi2 = 90 + da;
                }
                else
                {
                    azi1 = -90 - da;
                    azi2 = -90 + da;
                }
            }
            if (SignBit(dlon))
            {
                azi2 = azi1;
                s12 = -s12;
                sig12 = -sig12;
            }
            var x = s12;
            azi = AngNormalize(azi2);
            var perp = _earth.Line(lat, dlon, azi, GeodesicFlags.GeodesicScale);

            perp.GenPosition(true, -sig12,
                             GeodesicFlags.GeodesicScale,
                             out _, out _, out _, out _, out _, out _, out rk, out _);

            SinCosd(perp.EquatorialAzimuth, out var salp0, out var calp0);
            double
              sbet1 = lat >= 0 ? calp0 : -calp0,
              cbet1 = Abs(dlon) <= 90 ? Abs(salp0) : -Abs(salp0),
              sbet01 = sbet1 * _cbet0 - cbet1 * _sbet0,
              cbet01 = cbet1 * _cbet0 + sbet1 * _sbet0,
              sig01 = Atan2(sbet01, cbet01) / Degree;

            _meridian.GenPosition(true, sig01,
                                  GeodesicFlags.Distance,
                                  out _, out _, out _, out var y, out _, out _, out _, out _);

            return (x, y);
        }

        /// <summary>
        /// Reverse projection, from Cassini-Soldner to geographic.
        /// </summary>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <param name="azi">azimuth of easting direction at point (degrees).</param>
        /// <param name="rk">reciprocal of azimuthal northing scale at point.</param>
        /// <returns>
        /// <i>lat</i>, latitude of point and <i>lon</i>, longitude of point, in degress.
        /// </returns>
        /// <remarks>
        /// A call to <see cref="Reverse(double, double, out double, out double)"/> followed
        /// by a call to <see cref="Forward(double, double, out double, out double)"/> will 
        /// return the original (<paramref name="x"/>, <paramref name="y"/>) (to within roundoff), 
        /// provided that <paramref name="x"/> and <paramref name="y"/> are sufficiently small not to "wrap around" the earth. 
        /// The routine does nothing if the origin has not been set.
        /// </remarks>
        public (double lat, double lon) Reverse(double x, double y, out double azi, out double rk)
        {
            _meridian.Position(y, out var lat1, out var lon1, out var azi0);
            _earth.Direct(lat1, lon1, azi0 + 90, x, out var lat, out var lon, out azi, out rk, out _);

            return (lat, lon);
        }

        /// <summary>
        /// Forward without returning the azimuth and scale.
        /// </summary>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <returns>
        /// <i>x</i>, easting of point and <i>y</i>, northing of point, in meters.
        /// </returns>
        /// <remarks>
        /// <paramref name="lat"/> should be in the range [−90°, 90°].
        /// A call to <see cref="Forward(double, double)"/> 
        /// followed by a call to <see cref="Reverse(double, double)"/> 
        /// will return the original (<paramref name="lat"/>, <paramref name="lon"/>) (to within roundoff).
        /// The routine does nothing if the origin has not been set.
        /// </remarks>
        public (double x, double y) Forward(double lat, double lon) => Forward(lat, lon, out _, out _);

        /// <summary>
        /// Reverse without returning the azimuth and scale.
        /// </summary>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <returns>
        /// <i>lat</i>, latitude of point and <i>lon</i>, longitude of point, in degress.
        /// </returns>
        /// <remarks>
        /// A call to <see cref="Reverse(double, double)"/> followed
        /// by a call to <see cref="Forward(double, double)"/> will 
        /// return the original (<paramref name="x"/>, <paramref name="y"/>) (to within roundoff), 
        /// provided that <paramref name="x"/> and <paramref name="y"/> are sufficiently small not to "wrap around" the earth. 
        /// The routine does nothing if the origin has not been set.
        /// </remarks>
        public (double lat, double lon) Reverse(double x, double y) => Reverse(x, y, out _, out _);
    }
}
