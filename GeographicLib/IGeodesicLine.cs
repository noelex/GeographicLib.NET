namespace GeographicLib
{
    /// <summary>
    /// Represents a geodesic line.
    /// </summary>
    public interface IGeodesicLine : IEllipsoid
    {
        #region Properties

        /// <summary>
        /// Gets or sets a value representing the arc length to point 3 (degrees).
        /// </summary>
        double Arc { get; set; }

        /// <summary>
        /// Gets a value representing the azimuth (degrees) of the geodesic line at point 1.
        /// </summary>
        double Azimuth { get; }

        /// <summary>
        /// Gets a avalue representing the computational capabilities that this object was constructed with.
        /// <see cref="GeodesicFlags.Latitude"/> and <see cref="GeodesicFlags.Azimuth"/> are always included.
        /// </summary>
        GeodesicFlags Capabilities { get; }

        /// <summary>
        /// Gets a value representing cosine of <see cref="Azimuth"/>.
        /// </summary>
        double CosineAzimuth { get; }

        /// <summary>
        /// Gets a value representing cosine of <see cref="EquatorialAzimuth"/>.
        /// </summary>
        double CosineEquatorialAzimuth { get; }

        /// <summary>
        /// Gets or sets a value representing the distance to point 3 (meters).
        /// </summary>
        double Distance { get; set; }

        /// <summary>
        /// Gets a value representing the arc length (degrees) between the northward equatorial crossing and point 1.
        /// </summary>
        /// <remarks>
        /// The result lies in (−180°, 180°].
        /// </remarks>
        double EquatorialArc { get; }

        /// <summary>
        /// Gets a value representing the azimuth (degrees) of the geodesic line as it crosses
        /// the equator in a northward direction.
        /// </summary>
        /// <remarks>
        /// The result lies in [−90°, 90°].
        /// </remarks>
        double EquatorialAzimuth { get; }

        /// <summary>
        /// Gets a value representing the latitude of point 1 (degrees).
        /// </summary>
        double Latitude { get; }

        /// <summary>
        /// Gets a value representing the longitude of point 1 (degrees).
        /// </summary>
        double Longitude { get; }

        /// <summary>
        /// Gets a value representing sine of <see cref="Azimuth"/>.
        /// </summary>
        double SineAzimuth { get; }

        /// <summary>
        /// Gets a value representing sine of <see cref="EquatorialAzimuth"/>.
        /// </summary>
        double SineEquatorialAzimuth { get; }

        #endregion

        /// <summary>
        /// The general position function.
        /// <see cref="Position(double, out double, out double, out double, out double, out double, out double, out double)"/>
        /// and <see cref="ArcPosition(double, out double, out double, out double, out double, out double, out double, out double, out double)"/> are defined in terms of this function.
        /// </summary>
        /// <param name="arcmode">
        /// boolean flag determining the meaning of the second parameter; 
        /// if <paramref name="arcmode"/> is <see langword="false"/>, 
        /// then the <see cref="IGeodesicLine"/> object must have been constructed with 
        /// <i>caps</i> |= <see cref="GeodesicFlags.DistanceIn"/>.
        /// </param>
        /// <param name="s12_a12">
        /// if <i>arcmode</i> is false, this is the distance between point 1 and point 2 (meters);
        /// otherwise it is the arc length between point 1 and point 2 (degrees); it can be negative.
        /// </param>
        /// <param name="outmask">a bitor'ed combination of <see cref="GeodesicFlags"/> values specifying which of the following
        /// parameters should be set.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">	longitude of point 2 (degrees); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Longitude"/>.</param>
        /// <param name="azi2">	(forward) azimuth at point 2 (degrees).</param>
        /// <param name="s12">	distance from point 1 to point 2 (meters); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Distance"/>.</param>
        /// <param name="m12">	reduced length of geodesic (meters);
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.ReducedLength"/>.</param>
        /// <param name="M12">geodesic scale of point 2 relative to point 1 (dimensionless); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.GeodesicScale"/>.</param>
        /// <param name="M21">geodesic scale of point 1 relative to point 2 (dimensionless); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.GeodesicScale"/>.</param>
        /// <param name="S12">area under the geodesic (meters2);
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Area"/>.</param>
        /// <returns><i>a12</i>, arc length from point 1 to point 2 (degrees).</returns>
        /// <remarks>
        /// <para>
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
        /// </para>
        /// <para>
        /// Requesting a value which the <see cref="IGeodesicLine"/> object is not capable of computing is not an error; 
        /// the corresponding argument will not be altered. Note, however, that the arc length is always
        /// computed and returned as the function value.
        /// </para>
        /// <para>
        /// With the <see cref="GeodesicFlags.LongUnroll"/> bit set, the quantity <i>lon2</i> − <i>lon1</i> indicates
        /// how many times and in what sense the geodesic encircles the ellipsoid.
        /// </para>
        /// </remarks>
        double GenPosition(bool arcmode, double s12_a12, GeodesicFlags outmask, out double lat2, out double lon2, out double azi2, out double s12, out double m12, out double M12, out double M21, out double S12);

        /// <summary>
        /// Gets the distance or arc length to point 3.
        /// </summary>
        /// <param name="arcmode">boolean flag determining the meaning of returned value.</param>
        /// <returns><i>s13</i> if <paramref name="arcmode"/> is <see langword="false"/>; <i>a13</i> if <paramref name="arcmode"/> is <see langword="true"/>.</returns>
        double GetDistance(bool arcmode);

        /// <summary>
        /// Specify position of point 3 in terms of either distance or arc length.
        /// </summary>
        /// <param name="arcmode">
        /// boolean flag determining the meaning of the second parameter; 
        /// if <paramref name="arcmode"/> is <see langword="false"/>, then the <see cref="IGeodesicLine"/> object
        /// must have been constructed with <i>caps</i> |= <see cref="GeodesicFlags.DistanceIn"/>.
        /// </param>
        /// <param name="s13_a13">
        /// if <paramref name="arcmode"/> is <see langword="false"/>, this is the distance from point 1 to point 3 (meters);
        /// otherwise it is the arc length from point 1 to point 3 (degrees); it can be negative.</param>
        void SetDistance(bool arcmode, double s13_a13);

        /// <summary>
        /// Test what capabilities are available.
        /// </summary>
        /// <param name="testcaps">a set of bitor'ed <see cref="GeodesicFlags"/> values.</param>
        /// <returns><see langword="true"/> if the <see cref="IGeodesicLine"/> object has all these capabilities.</returns>
        bool HasCapability(GeodesicFlags testcaps);

        #region Position in terms of arc length

        /// <summary>
        /// Compute the position of point 2 which is an arc length <i>a12</i> (degrees) from point 1.
        /// </summary>
        /// <param name="a12">arc length from point 1 to point 2 (degrees); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">	longitude of point 2 (degrees); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Longitude"/>.</param>
        /// <param name="azi2">	(forward) azimuth at point 2 (degrees).</param>
        /// <param name="s12">	distance from point 1 to point 2 (meters); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Distance"/>.</param>
        /// <param name="m12">	reduced length of geodesic (meters);
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.ReducedLength"/>.</param>
        /// <param name="M12">geodesic scale of point 2 relative to point 1 (dimensionless); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.GeodesicScale"/>.</param>
        /// <param name="M21">geodesic scale of point 1 relative to point 2 (dimensionless); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.GeodesicScale"/>.</param>
        /// <param name="S12">area under the geodesic (meters2);
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Area"/>.</param>
        /// <remarks>
        /// <para>
        /// The values of <i>lon2</i> and <i>azi2</i> returned are in the range [−180°, 180°].
        /// </para>
        /// <para>
        /// Requesting a value which the <see cref="IGeodesicLine"/> object is not capable of computing is not an error;
        /// the corresponding argument will not be altered.
        /// </para>
        /// <para>
        /// The following functions are overloaded versions of <see cref="ArcPosition(double, out double, out double, out double, out double, out double, out double, out double, out double)"/>
        /// which omit some of the output parameters.
        /// </para>
        /// </remarks>
        void ArcPosition(double a12, out double lat2, out double lon2, out double azi2,
                     out double s12, out double m12, out double M12, out double M21,
                     out double S12);

        /// <summary>
        /// Compute the position of point 2 which is an arc length <i>a12</i> (degrees) from point 1.
        /// </summary>
        /// <param name="a12">arc length from point 1 to point 2 (degrees); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">	longitude of point 2 (degrees); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Longitude"/>.</param>
        void ArcPosition(double a12, out double lat2, out double lon2);

        /// <summary>
        /// Compute the position of point 2 which is an arc length <i>a12</i> (degrees) from point 1.
        /// </summary>
        /// <param name="a12">arc length from point 1 to point 2 (degrees); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">	longitude of point 2 (degrees); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Longitude"/>.</param>
        /// <param name="azi2">	(forward) azimuth at point 2 (degrees).</param>
        void ArcPosition(double a12,
                         out double lat2, out double lon2, out double azi2);

        /// <summary>
        /// Compute the position of point 2 which is an arc length <i>a12</i> (degrees) from point 1.
        /// </summary>
        /// <param name="a12">arc length from point 1 to point 2 (degrees); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">	longitude of point 2 (degrees); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Longitude"/>.</param>
        /// <param name="azi2">	(forward) azimuth at point 2 (degrees).</param>
        /// <param name="s12">	distance from point 1 to point 2 (meters); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Distance"/>.</param>
        void ArcPosition(double a12, out double lat2, out double lon2, out double azi2,
                          out double s12);

        /// <summary>
        /// Compute the position of point 2 which is an arc length <i>a12</i> (degrees) from point 1.
        /// </summary>
        /// <param name="a12">arc length from point 1 to point 2 (degrees); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">	longitude of point 2 (degrees); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Longitude"/>.</param>
        /// <param name="azi2">	(forward) azimuth at point 2 (degrees).</param>
        /// <param name="s12">	distance from point 1 to point 2 (meters); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Distance"/>.</param>
        /// <param name="m12">	reduced length of geodesic (meters);
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.ReducedLength"/>.</param>
        void ArcPosition(double a12, out double lat2, out double lon2, out double azi2,
                          out double s12, out double m12);

        /// <summary>
        /// Compute the position of point 2 which is an arc length <i>a12</i> (degrees) from point 1.
        /// </summary>
        /// <param name="a12">arc length from point 1 to point 2 (degrees); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">	longitude of point 2 (degrees); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Longitude"/>.</param>
        /// <param name="azi2">	(forward) azimuth at point 2 (degrees).</param>
        /// <param name="s12">	distance from point 1 to point 2 (meters); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Distance"/>.</param>
        /// <param name="M12">geodesic scale of point 2 relative to point 1 (dimensionless); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.GeodesicScale"/>.</param>
        /// <param name="M21">geodesic scale of point 1 relative to point 2 (dimensionless); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.GeodesicScale"/>.</param>
        void ArcPosition(double a12, out double lat2, out double lon2, out double azi2,
                         out double s12, out double M12, out double M21);

        /// <summary>
        /// Compute the position of point 2 which is an arc length <i>a12</i> (degrees) from point 1.
        /// </summary>
        /// <param name="a12">arc length from point 1 to point 2 (degrees); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">	longitude of point 2 (degrees); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Longitude"/>.</param>
        /// <param name="azi2">	(forward) azimuth at point 2 (degrees).</param>
        /// <param name="s12">	distance from point 1 to point 2 (meters); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Distance"/>.</param>
        /// <param name="m12">	reduced length of geodesic (meters);
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.ReducedLength"/>.</param>
        /// <param name="M12">geodesic scale of point 2 relative to point 1 (dimensionless); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.GeodesicScale"/>.</param>
        /// <param name="M21">geodesic scale of point 1 relative to point 2 (dimensionless); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.GeodesicScale"/>.</param>
        void ArcPosition(double a12, out double lat2, out double lon2, out double azi2,
                         out double s12, out double m12, out double M12, out double M21);

        #endregion

        #region Position in terms of distance

        /// <summary>
        /// Compute the position of point 2 which is a distance <i>s12</i> (meters) from point 1.
        /// </summary>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">	longitude of point 2 (degrees); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Longitude"/>.</param>
        /// <param name="azi2">	(forward) azimuth at point 2 (degrees).</param>
        /// <param name="s12">	distance from point 1 to point 2 (meters); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Distance"/>.</param>
        /// <param name="m12">	reduced length of geodesic (meters);
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.ReducedLength"/>.</param>
        /// <param name="M12">geodesic scale of point 2 relative to point 1 (dimensionless); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.GeodesicScale"/>.</param>
        /// <param name="M21">geodesic scale of point 1 relative to point 2 (dimensionless); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.GeodesicScale"/>.</param>
        /// <param name="S12">area under the geodesic (meters2);
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Area"/>.</param>
        /// <returns><i>a12</i>, arc length from point 1 to point 2 (degrees).</returns>
        /// <remarks>
        /// The values of <i>lon2</i> and <i>azi2</i> returned are in the range [−180°, 180°].
        /// <para>
        /// The <see cref="IGeodesicLine"/> object must have been constructed with
        /// <i>caps</i> |= <see cref="GeodesicFlags.DistanceIn"/>; otherwise <see cref="double.NaN"/> is returned 
        /// and no parameters are set. Requesting a value which the <see cref="IGeodesicLine"/> object is not capable
        /// of computing is not an error; the corresponding argument will not be altered.
        /// </para>
        /// <para>
        /// The following functions are overloaded versions of <see cref="Position(double, out double, out double, out double, out double, out double, out double, out double)"/>
        /// which omit some of the output parameters.
        /// Note, however, that the arc length is always computed and returned as the function value.
        /// </para>
        /// </remarks>
        double Position(double s12,
                        out double lat2, out double lon2, out double azi2,
                        out double m12, out double M12, out double M21,
                        out double S12);

        /// <summary>
        /// Compute the position of point 2 which is a distance <i>s12</i> (meters) from point 1.
        /// </summary>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">	longitude of point 2 (degrees); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Longitude"/>.</param>
        /// <param name="s12">	distance from point 1 to point 2 (meters); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Distance"/>.</param>
        /// <returns><i>a12</i>, arc length from point 1 to point 2 (degrees).</returns>
        double Position(double s12, out double lat2, out double lon2);

        /// <summary>
        /// Compute the position of point 2 which is a distance <i>s12</i> (meters) from point 1.
        /// </summary>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">	longitude of point 2 (degrees); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Longitude"/>.</param>
        /// <param name="azi2">	(forward) azimuth at point 2 (degrees).</param>
        /// <param name="s12">	distance from point 1 to point 2 (meters); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Distance"/>.</param> 
        /// <returns><i>a12</i>, arc length from point 1 to point 2 (degrees).</returns>
        double Position(double s12, out double lat2, out double lon2,
                            out double azi2);

        /// <summary>
        /// Compute the position of point 2 which is a distance <i>s12</i> (meters) from point 1.
        /// </summary>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">	longitude of point 2 (degrees); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Longitude"/>.</param>
        /// <param name="azi2">	(forward) azimuth at point 2 (degrees).</param>
        /// <param name="s12">	distance from point 1 to point 2 (meters); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Distance"/>.</param>
        /// <param name="m12">	reduced length of geodesic (meters);
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.ReducedLength"/>.</param>        /// <returns><i>a12</i>, arc length from point 1 to point 2 (degrees).</returns>
        double Position(double s12, out double lat2, out double lon2,
                            out double azi2, out double m12);

        /// <summary>
        /// Compute the position of point 2 which is a distance <i>s12</i> (meters) from point 1.
        /// </summary>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">	longitude of point 2 (degrees); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Longitude"/>.</param>
        /// <param name="azi2">	(forward) azimuth at point 2 (degrees).</param>
        /// <param name="s12">	distance from point 1 to point 2 (meters); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Distance"/>.</param>
        /// <param name="M12">geodesic scale of point 2 relative to point 1 (dimensionless); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.GeodesicScale"/>.</param>
        /// <param name="M21">geodesic scale of point 1 relative to point 2 (dimensionless); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.GeodesicScale"/>.</param>
        /// <returns><i>a12</i>, arc length from point 1 to point 2 (degrees).</returns>
        double Position(double s12, out double lat2, out double lon2,
                            out double azi2, out double M12, out double M21);

        /// <summary>
        /// Compute the position of point 2 which is a distance <i>s12</i> (meters) from point 1.
        /// </summary>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">	longitude of point 2 (degrees); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Longitude"/>.</param>
        /// <param name="azi2">	(forward) azimuth at point 2 (degrees).</param>
        /// <param name="s12">	distance from point 1 to point 2 (meters); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.Distance"/>.</param>
        /// <param name="m12">	reduced length of geodesic (meters);
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.ReducedLength"/>.</param>
        /// <param name="M12">geodesic scale of point 2 relative to point 1 (dimensionless); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.GeodesicScale"/>.</param>
        /// <param name="M21">geodesic scale of point 1 relative to point 2 (dimensionless); 
        /// requires that the <see cref="IGeodesicLine"/> object was constructed with <i>caps</i> |= <see cref="GeodesicFlags.GeodesicScale"/>.</param>
        /// <returns><i>a12</i>, arc length from point 1 to point 2 (degrees).</returns>
        double Position(double s12,
                            out double lat2, out double lon2, out double azi2,
                            out double m12, out double M12, out double M21);

        #endregion

    }
}