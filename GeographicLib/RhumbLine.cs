﻿namespace GeographicLib
{
    /// <summary>
    /// Find a sequence of points on a single rhumb line.
    /// </summary>
    /// <remarks>
    /// <see cref="RhumbLine"/> facilitates the determination of a series of points on a single rhumb line.
    /// The starting point (<i>lat1</i>, <i>lon1</i>) and the azimuth <i>azi12</i> are specified in the call to <see cref="Rhumb.Line(double, double, double)"/>
    /// which returns a <see cref="RhumbLine"/> object. <see cref="Position(double, out double, out double)"/> returns the location of point 2
    /// (and, optionally, the corresponding area, <i>S12</i>) a distance <i>s12</i> along the rhumb line.
    /// <para>
    /// There is no public constructor for this class.
    /// (Use <see cref="Rhumb.Line(double, double, double)"/> to create an instance.)
    /// The <see cref="Rhumb"/> object used to create a <see cref="RhumbLine"/> must stay in scope as long as the <see cref="RhumbLine"/>.
    /// </para>
    /// </remarks>
    public partial class RhumbLine : IEllipsoid
    {
        private Priv _priv;

        internal RhumbLine(Rhumb rh, double lat1, double lon1, double azi12)
        {
            _priv = new Priv(rh, lat1, lon1, azi12);
        }

        /// <summary>
        /// The general position routine. <see cref="Position(double, out double, out double)"/> is defined in term so this function.
        /// </summary>
        /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
        /// <param name="outmask">a bitor'ed combination of <see cref="GeodesicFlags"/> values specifying which of the following parameters should be set.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="S12">area under the rhumb line (meters^2).</param>
        /// <remarks>
        /// The <see cref="GeodesicFlags"/> values possible for <paramref name="outmask"/> are
        /// <list type="bullet">
        /// <item>outmask |= <see cref="GeodesicFlags.Latitude"/> for the latitude <paramref name="lat2"/>;</item>
        /// <item>outmask |= <see cref="GeodesicFlags.Longitude"/> for the longitude <paramref name="lon2"/>;</item>
        /// <item>outmask |= <see cref="GeodesicFlags.Area"/> for the area <paramref name="S12"/>;</item>
        /// <item>outmask |= <see cref="GeodesicFlags.All"/> for all of the above;</item>
        /// <item>outmask |= <see cref="GeodesicFlags.LongUnroll"/> to unroll <paramref name="lon2"/> instead of wrapping it into the range [−180°, 180°].</item>
        /// </list>
        /// With the <see cref="GeodesicFlags.LongUnroll"/> bit set, the quantity <paramref name="lon2"/> − <i>lon1</i> indicates
        /// how many times and in what sense the rhumb line encircles the ellipsoid.
        /// <para>
        /// If s12 is large enough that the rhumb line crosses a pole,
        /// the longitude of point 2 is indeterminate (a <see cref="double.NaN"/> is returned for <paramref name="lon2"/> and <paramref name="S12"/>).
        /// </para>
        /// </remarks>
        public void GenPosition(double s12, GeodesicFlags outmask, out double lat2, out double lon2, out double S12)
        {
            _priv.GenPosition(s12, outmask, out lat2, out lon2, out S12);
        }

        /// <summary>
        /// Compute the position of point 2 which is a distance <paramref name="s12"/> (meters) from point 1.
        /// </summary>
        /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
        /// <param name="outmask">
        /// a bitor'ed combination of <see cref="GeodesicFlags"/> values specifying
        /// which of the properties in returned <see cref="DirectRhumbResult"/> instance should be set.
        /// </param>
        /// <returns>A <see cref="DirectRhumbResult"/> instance containing the result of the calcutation.</returns>
        /// <remarks>
        /// The <see cref="GeodesicFlags"/> values possible for <paramref name="outmask"/> are
        /// <list type="bullet">
        /// <item><i>outmask</i> |= <see cref="GeodesicFlags.Latitude"/> for the latitude returned in <see cref="DirectRhumbResult.Latitude"/>;</item>
        /// <item><i>outmask</i> |= <see cref="GeodesicFlags.Longitude"/> for the longitude returned in  <see cref="DirectRhumbResult.Longitude"/>;</item>
        /// <item><i>outmask</i> |= <see cref="GeodesicFlags.Area"/> for the area returned in  <see cref="RhumbResult.Area"/>;</item>
        /// <item><i>outmask</i> |= <see cref="GeodesicFlags.All"/> for all of the above;</item>
        /// <item><i>outmask</i> |= <see cref="GeodesicFlags.LongUnroll"/> to unroll <see cref="DirectRhumbResult.Longitude"/> instead of wrapping it into the range [−180°, 180°].</item>
        /// </list>
        /// With the <see cref="GeodesicFlags.LongUnroll"/> bit set, the quantity <see cref="DirectRhumbResult.Longitude"/> − <i>lon1</i> indicates
        /// how many times and in what sense the rhumb line encircles the ellipsoid.
        /// <para>
        /// If <paramref name="s12"/> is large enough that the rhumb line crosses a pole,
        /// the longitude of point 2 is indeterminate (a <see cref="double.NaN"/> is returned for <see cref="DirectRhumbResult.Longitude"/> and <see cref="RhumbResult.Area"/>).
        /// </para>
        /// </remarks>
        public DirectRhumbResult Position(double s12, GeodesicFlags outmask = GeodesicFlags.All)
        {
            GenPosition(s12, outmask, out var lat2, out var lon2, out var S12);
            return new DirectRhumbResult
            {
                Latitude = lat2,
                Longitude = lon2,
                Area = S12
            };
        }

        /// <summary>
        /// Compute the position of point 2 which is a distance <paramref name="s12"/> (meters) from point 1. The area is not computed.
        /// </summary>
        /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <remarks>
        /// The value of <paramref name="lon2"/> returned is in the range [−180°, 180°].
        /// <para>
        /// If s12 is large enough that the rhumb line crosses a pole,
        /// the longitude of point 2 is indeterminate (a <see cref="double.NaN"/> is returned for <paramref name="lon2"/>).
        /// </para>
        /// </remarks>
        public void Position(double s12, out double lat2, out double lon2)
            => GenPosition(s12, GeodesicFlags.Latitude | GeodesicFlags.Longitude, out lat2, out lon2, out _);

        /// <summary>
        /// Compute the position of point 2 which is a distance <paramref name="s12"/> (meters) from point 1. The area is also computed.
        /// </summary>
        /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="S12">area under the rhumb line (meters^2).</param>
        /// <remarks>
        /// The value of <paramref name="lon2"/> returned is in the range [−180°, 180°].
        /// <para>
        /// If s12 is large enough that the rhumb line crosses a pole,
        /// the longitude of point 2 is indeterminate (a <see cref="double.NaN"/> is returned for <paramref name="lon2"/> and <paramref name="S12"/>).
        /// </para>
        /// </remarks>
        public void Position(double s12, out double lat2, out double lon2, out double S12)
            => GenPosition(s12, GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Area, out lat2, out lon2, out S12);

        /// <summary>
        /// Gets a value representing the latitude of point 1 in degrees (<i>lat1</i>).
        /// </summary>
        public double Latitude => _priv._lat1;

        /// <summary>
        /// Gets a value representing the longitude of point 1 in degrees (<i>lon1</i>).
        /// </summary>
        public double Longitude => _priv._lon1;

        /// <summary>
        /// Gets a value representing the azimuth of the rhumb line in degrees (<i>azi12</i>).
        /// </summary>
        public double Azimuth => _priv._azi12;

        /// <inheritdoc/>
        public double EquatorialRadius => _priv._rh.EquatorialRadius;

        /// <inheritdoc/>
        public double Flattening => _priv._rh.Flattening;
    }
}
