namespace GeographicLib
{
    /// <summary>
    /// Exposes geodesic calculations.
    /// </summary>
    public interface IGeodesic : IGeodesicLike
    {
        #region Direct geodesic problem specified in terms of arc length

        /// <summary>
        /// Solve the direct geodesic problem where the length of the geodesic is specified in terms of arc length.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="a12">arc length of between point 1 and point 2 (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="azi2">(forward) azimuth at point 2 (degrees).</param>
        /// <param name="m12">reduced length of geodesic (meters).</param>
        /// <param name="M12">geodesic scale of point 2 relative to point 1 (dimensionless).</param>
        /// <param name="M21">geodesic scale of point 1 relative to point 2 (dimensionless).</param>
        /// <param name="S12">area under the geodesic (meters^2).</param>
        /// <remarks>
        /// <para>
        /// <i>lat1</i> should be in the range [−90°, 90°]. 
        /// The values of <i>lon2</i> and <i>azi2</i> returned are in the range [−180°, 180°].
        /// </para>
        /// <para>
        /// If either point is at a pole, the azimuth is defined by keeping the longitude fixed,
        /// writing lat = ±(90° − ε), and taking the limit ε → 0+.
        /// An arc length greater that 180° signifies a geodesic which is not a shortest path.
        /// (For a prolate ellipsoid, an additional condition is necessary for a shortest path:
        /// the longitudinal extent must not exceed of 180°.)
        /// </para>
        /// <para>
        /// The following functions are overloaded versions of 
        /// <see cref="ArcDirect(double, double, double, double, out double, out double, out double, out double, out double, out double, out double, out double)"/>
        /// which omit some of the output parameters.
        /// </para>
        /// </remarks>
        void ArcDirect(double lat1, double lon1, double azi1, double a12,
                out double lat2, out double lon2, out double azi2, out double s12,
                out double m12, out double M12, out double M21, out double S12);

        /// <summary>
        /// Solve the direct geodesic problem where the length of the geodesic is specified in terms of arc length.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="a12">arc length of between point 1 and point 2 (degrees).</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="azi2">(forward) azimuth at point 2 (degrees).</param>
        void ArcDirect(double lat1, double lon1, double azi1, double a12,
                   out double lat2, out double lon2, out double azi2);

        /// <summary>
        /// Solve the direct geodesic problem where the length of the geodesic is specified in terms of arc length.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="a12">arc length of between point 1 and point 2 (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="azi2">(forward) azimuth at point 2 (degrees).</param>
        void ArcDirect(double lat1, double lon1, double azi1, double a12,
                   out double lat2, out double lon2, out double azi2, out double s12);

        /// <summary>
        /// Solve the direct geodesic problem where the length of the geodesic is specified in terms of arc length.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="a12">arc length of between point 1 and point 2 (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="azi2">(forward) azimuth at point 2 (degrees).</param>
        /// <param name="m12">reduced length of geodesic (meters).</param>
        void ArcDirect(double lat1, double lon1, double azi1, double a12,
                       out double lat2, out double lon2, out double azi2,
                       out double s12, out double m12);

        /// <summary>
        /// Solve the direct geodesic problem where the length of the geodesic is specified in terms of arc length.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="a12">arc length of between point 1 and point 2 (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="azi2">(forward) azimuth at point 2 (degrees).</param>
        /// <param name="M12">geodesic scale of point 2 relative to point 1 (dimensionless).</param>
        /// <param name="M21">geodesic scale of point 1 relative to point 2 (dimensionless).</param>
        void ArcDirect(double lat1, double lon1, double azi1, double a12,
                       out double lat2, out double lon2, out double azi2, out double s12,
                       out double M12, out double M21);

        /// <summary>
        /// Solve the direct geodesic problem where the length of the geodesic is specified in terms of arc length.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="a12">arc length of between point 1 and point 2 (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="azi2">(forward) azimuth at point 2 (degrees).</param>
        /// <param name="m12">reduced length of geodesic (meters).</param>
        /// <param name="M12">geodesic scale of point 2 relative to point 1 (dimensionless).</param>
        /// <param name="M21">geodesic scale of point 1 relative to point 2 (dimensionless).</param>
        void ArcDirect(double lat1, double lon1, double azi1, double a12,
                       out double lat2, out double lon2, out double azi2, out double s12,
                       out double m12, out double M12, out double M21);

        #endregion

        #region Direct geodesic problem specified in terms of distance

        /// <summary>
        /// Solve the direct geodesic problem where the length of the geodesic is specified in terms of distance.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <returns><i>a12</i>, arc length of between point 1 and point 2 (degrees).</returns>
        double Direct(double lat1, double lon1, double azi1, double s12, out double lat2, out double lon2);

        /// <summary>
        /// Solve the direct geodesic problem where the length of the geodesic is specified in terms of distance.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="azi2">(forward) azimuth at point 2 (degrees).</param>
        /// <returns><i>a12</i>, arc length of between point 1 and point 2 (degrees).</returns>
        double Direct(double lat1, double lon1, double azi1, double s12, out double lat2, out double lon2, out double azi2);

        /// <summary>
        /// Solve the direct geodesic problem where the length of the geodesic is specified in terms of distance.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="azi2">(forward) azimuth at point 2 (degrees).</param>
        /// <param name="m12">reduced length of geodesic (meters).</param>
        /// <returns><i>a12</i>, arc length of between point 1 and point 2 (degrees).</returns>
        double Direct(double lat1, double lon1, double azi1, double s12, out double lat2, out double lon2, out double azi2, out double m12);

        /// <summary>
        /// Solve the direct geodesic problem where the length of the geodesic is specified in terms of distance.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="azi2">(forward) azimuth at point 2 (degrees).</param>
        /// <param name="M12">geodesic scale of point 2 relative to point 1 (dimensionless).</param>
        /// <param name="M21">geodesic scale of point 1 relative to point 2 (dimensionless).</param>
        /// <returns><i>a12</i>, arc length of between point 1 and point 2 (degrees).</returns>
        double Direct(double lat1, double lon1, double azi1, double s12, out double lat2, out double lon2, out double azi2, out double M12, out double M21);

        /// <summary>
        /// Solve the direct geodesic problem where the length of the geodesic is specified in terms of distance.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="azi2">(forward) azimuth at point 2 (degrees).</param>
        /// <param name="m12">reduced length of geodesic (meters).</param>
        /// <param name="M12">geodesic scale of point 2 relative to point 1 (dimensionless).</param>
        /// <param name="M21">geodesic scale of point 1 relative to point 2 (dimensionless).</param>
        /// <returns><i>a12</i>, arc length of between point 1 and point 2 (degrees).</returns>
        double Direct(double lat1, double lon1, double azi1, double s12,
            out double lat2, out double lon2, out double azi2, out double m12, out double M12, out double M21);

        /// <summary>
        /// Solve the direct geodesic problem where the length of the geodesic is specified in terms of distance.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="azi2">(forward) azimuth at point 2 (degrees).</param>
        /// <param name="m12">reduced length of geodesic (meters).</param>
        /// <param name="M12">geodesic scale of point 2 relative to point 1 (dimensionless).</param>
        /// <param name="M21">geodesic scale of point 1 relative to point 2 (dimensionless).</param>
        /// <param name="S12">area under the geodesic (meters^2).</param>
        /// <returns><i>a12</i>, arc length of between point 1 and point 2 (degrees).</returns>
        /// <remarks>
        /// <i>lat1</i> should be in the range [−90°, 90°]. The values of <i>lon2</i> and <i>azi2</i> returned are in the range [−180°, 180°].
        /// <para>
        /// If either point is at a pole, the azimuth is defined by keeping the longitude fixed, writing <i>lat</i> = ±(90° − ε), 
        /// and taking the limit ε → 0+. An arc length greater that 180° signifies a geodesic which is not a shortest path.
        /// (For a prolate ellipsoid, an additional condition is necessary for a shortest path: the longitudinal extent must not exceed of 180°.)
        /// </para>
        /// <para>
        /// The following functions are overloaded versions of 
        /// <see cref="Direct(double, double, double, double, out double, out double, out double, out double, out double, out double, out double)"/>
        /// which omit some of the output parameters. Note, however, that the arc length is always computed and returned as the function value.
        /// </para>
        /// </remarks>
        double Direct(double lat1, double lon1, double azi1, double s12,
            out double lat2, out double lon2, out double azi2, out double m12, out double M12, out double M21, out double S12);
        #endregion

        #region Inverse geodesic problem

        /// <summary>
        /// Solve the inverse geodesic problem.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="azi2">(forward) azimuth at point 2 (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters).</param>
        /// <param name="m12">reduced length of geodesic (meters).</param>
        /// <param name="M12">geodesic scale of point 2 relative to point 1 (dimensionless).</param>
        /// <param name="M21">geodesic scale of point 1 relative to point 2 (dimensionless).</param>
        /// <param name="S12">area under the geodesic (meters^2).</param>
        /// <returns><i>a12</i>, arc length of between point 1 and point 2 (degrees).</returns>
        /// <remarks>
        /// <para>
        /// <i>lat1</i> and <i>lat2</i> should be in the range [−90°, 90°]. 
        /// The values of <i>azi1</i> and <i>azi2</i> returned are in the range [−180°, 180°].
        /// </para>
        /// <para>
        /// If either point is at a pole, the azimuth is defined by keeping the longitude fixed,
        /// writing lat = ±(90° − ε), and taking the limit ε → 0+.
        /// </para>
        /// <para>
        /// The solution to the inverse problem is found using Newton's method. 
        /// If this fails to converge (this is very unlikely in geodetic applications but does occur for very eccentric ellipsoids), 
        /// then the bisection method is used to refine the solution.
        /// </para>
        /// <para>
        /// The following functions are overloaded versions of <see cref="Inverse(double, double, double, double, out double)"/>
        /// which omit some of the output parameters.
        /// Note, however, that the arc length is always computed and returned as the function value.
        /// </para>
        /// </remarks>
        double Inverse(double lat1, double lon1, double lat2, double lon2,
                       out double s12, out double azi1, out double azi2, out double m12,
                       out double M12, out double M21, out double S12);

        /// <summary>
        /// Solve the inverse geodesic problem.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters).</param>
        /// <returns><i>a12</i>, arc length of between point 1 and point 2 (degrees).</returns>
        double Inverse(double lat1, double lon1, double lat2, double lon2,
                           out double s12);

        /// <summary>
        /// Solve the inverse geodesic problem.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="azi2">(forward) azimuth at point 2 (degrees).</param>
        /// <returns><i>a12</i>, arc length of between point 1 and point 2 (degrees).</returns>
        double Inverse(double lat1, double lon1, double lat2, double lon2,
                       out double azi1, out double azi2);

        /// <summary>
        /// Solve the inverse geodesic problem.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="azi2">(forward) azimuth at point 2 (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters).</param>
        /// <returns><i>a12</i>, arc length of between point 1 and point 2 (degrees).</returns>
        double Inverse(double lat1, double lon1, double lat2, double lon2,
                           out double s12, out double azi1, out double azi2);

        /// <summary>
        /// Solve the inverse geodesic problem.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="azi2">(forward) azimuth at point 2 (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters).</param>
        /// <param name="m12">reduced length of geodesic (meters).</param>
        /// <returns><i>a12</i>, arc length of between point 1 and point 2 (degrees).</returns>
        double Inverse(double lat1, double lon1, double lat2, double lon2,
                           out double s12, out double azi1, out double azi2, out double m12);

        /// <summary>
        /// Solve the inverse geodesic problem.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="azi2">(forward) azimuth at point 2 (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters).</param>
        /// <param name="M12">geodesic scale of point 2 relative to point 1 (dimensionless).</param>
        /// <param name="M21">geodesic scale of point 1 relative to point 2 (dimensionless).</param>
        /// <returns><i>a12</i>, arc length of between point 1 and point 2 (degrees).</returns>
        double Inverse(double lat1, double lon1, double lat2, double lon2,
                           out double s12, out double azi1, out double azi2,
                           out double M12, out double M21);

        /// <summary>
        /// Solve the inverse geodesic problem.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="azi2">(forward) azimuth at point 2 (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters).</param>
        /// <param name="m12">reduced length of geodesic (meters).</param>
        /// <param name="M12">geodesic scale of point 2 relative to point 1 (dimensionless).</param>
        /// <param name="M21">geodesic scale of point 1 relative to point 2 (dimensionless).</param>
        /// <returns><i>a12</i>, arc length of between point 1 and point 2 (degrees).</returns>
        double Inverse(double lat1, double lon1, double lat2, double lon2,
                           out double s12, out double azi1, out double azi2, out double m12,
                           out double M12, out double M21);

        #endregion

        #region Interface to GeodesicLine

        /// <summary>
        /// Define a <see cref="IGeodesicLine"/> in terms of the direct geodesic problem
        /// specified in terms of either distance or arc length.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="arcmode">boolean flag determining the meaning of the <i>s12_a12</i>.</param>
        /// <param name="s12_a12">
        /// if <paramref name="arcmode"/> is false, this is the distance between point 1 and point 2 (meters);
        /// otherwise it is the arc length between point 1 and point 2 (degrees); it can be negative.
        /// </param>
        /// <param name="caps">
        /// bitor'ed combination of <see cref="GeodesicFlags"/> values specifying the capabilities
        /// the <see cref="IGeodesicLine"/> object should possess, i.e., which quantities can be returned
        /// in calls to <see cref="IGeodesicLine.Position(double, out double, out double, out double, out double, out double, out double, out double)"/>.
        /// </param>
        /// <returns>a <see cref="IGeodesicLine"/> object.</returns>
        /// <remarks>
        /// This function sets point 3 of the <see cref="IGeodesicLine"/> to correspond to point 2 of the direct geodesic problem.
        /// <para><paramref name="lat1"/> should be in the range [−90°, 90°].</para>
        /// </remarks>
        IGeodesicLine GenDirectLine(double lat1, double lon1, double azi1,
                           bool arcmode, double s12_a12,
                           GeodesicFlags caps = GeodesicFlags.All);

        /// <summary>
        /// Set up to compute several points on a single geodesic.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="caps">
        /// bitor'ed combination of <see cref="GeodesicFlags"/> values specifying the capabilities
        /// the <see cref="IGeodesicLine"/> object should possess, i.e., which quantities can be returned
        /// in calls to <see cref="IGeodesicLine.Position(double, out double, out double, out double, out double, out double, out double, out double)"/>.
        /// </param>
        /// <returns>a <see cref="IGeodesicLine"/> object.</returns>
        /// <remarks>
        /// <para><paramref name="lat1"/> should be in the range [−90°, 90°].</para>
        /// <para>
        /// The <see cref="GeodesicFlags"/> values are
        /// <list type="bullet">
        /// <item><i>caps</i> |= <see cref="GeodesicFlags.Latitude"/> for the latitude <i>lat2</i>;</item>
        /// <item><i>caps</i> |= <see cref="GeodesicFlags.Longitude"/> for the latitude <i>lon2</i>;</item>
        /// <item><i>caps</i> |= <see cref="GeodesicFlags.Azimuth"/> for the latitude <i>azi2</i>;</item>
        /// <item><i>caps</i> |= <see cref="GeodesicFlags.Distance"/> for the distance <i>s12</i>;</item>
        /// <item><i>caps</i> |= <see cref="GeodesicFlags.ReducedLength"/> for the reduced length <i>m12</i>;</item>
        /// <item><i>caps</i> |= <see cref="GeodesicFlags.GeodesicScale"/> for the geodesic scales <i>M12</i> and <i>M21</i>;</item>
        /// <item><i>caps</i> |= <see cref="GeodesicFlags.Area"/> for the area <i>S12</i>;</item>
        /// <item><i>caps</i> |= <see cref="GeodesicFlags.All"/> for all of the above;</item>
        /// <item><i>caps</i> |= <see cref="GeodesicFlags.LongUnroll"/> to unroll <i>lon2</i> instead of wrapping it into the range [−180°, 180°].</item>
        /// </list>
        /// </para>
        /// <para>The default value of caps is <see cref="GeodesicFlags.All"/>.</para>
        /// <para>
        /// If the point is at a pole, the azimuth is defined by keeping lon1 fixed, 
        /// writing <i>lat1</i> = ±(90 − ε), and taking the limit ε → 0+.
        /// </para>
        /// </remarks>
        IGeodesicLine Line(double lat1, double lon1, double azi1, GeodesicFlags caps = GeodesicFlags.All);

        /// <summary>
        /// Define a <see cref="IGeodesicLine"/> in terms of the inverse geodesic problem.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="lat2">latitude of point 2 (degrees).</param>
        /// <param name="lon2">longitude of point 2 (degrees).</param>
        /// <param name="caps">
        /// bitor'ed combination of <see cref="GeodesicFlags"/> values specifying the capabilities
        /// the <see cref="IGeodesicLine"/> object should possess, i.e., which quantities can be returned
        /// in calls to <see cref="IGeodesicLine.Position(double, out double, out double, out double, out double, out double, out double, out double)"/>.
        /// </param>
        /// <returns>a <see cref="IGeodesicLine"/> object.</returns>
        /// <remarks>
        /// <para>
        /// This function sets point 3 of the <see cref="IGeodesicLine"/> to correspond to point 2 of the inverse geodesic problem.
        /// </para>
        /// <para>
        /// <i>lat1</i> and <i>lat2</i> should be in the range [−90°, 90°].
        /// </para>
        /// </remarks>
        IGeodesicLine InverseLine(double lat1, double lon1, double lat2, double lon2, GeodesicFlags caps = GeodesicFlags.All);

        /// <summary>
        /// Define a <see cref="IGeodesicLine"/> in terms of the direct geodesic problem specified in terms of distance.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
        /// <param name="caps">
        /// bitor'ed combination of <see cref="GeodesicFlags"/> values specifying the capabilities
        /// the <see cref="IGeodesicLine"/> object should possess, i.e., which quantities can be returned
        /// in calls to <see cref="IGeodesicLine.Position(double, out double, out double, out double, out double, out double, out double, out double)"/>.
        /// </param>
        /// <returns>a <see cref="IGeodesicLine"/> object.</returns>
        /// <remarks>
        /// <para>
        /// This function sets point 3 of the <see cref="IGeodesicLine"/> to correspond to point 2 of the direct geodesic problem.
        /// </para>
        /// <para><i>lat1</i> should be in the range [−90°, 90°].</para>
        /// </remarks>
        IGeodesicLine DirectLine(double lat1, double lon1, double azi1, double s12,
            GeodesicFlags caps = GeodesicFlags.All);

        /// <summary>
        /// Define a <see cref="IGeodesicLine"/> in terms of the direct geodesic problem specified in terms of arc length.
        /// </summary>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="a12">arc length between point 1 and point 2 (degrees); it can be negative.</param>
        /// <param name="caps">
        /// bitor'ed combination of <see cref="GeodesicFlags"/> values specifying the capabilities
        /// the <see cref="IGeodesicLine"/> object should possess, i.e., which quantities can be returned
        /// in calls to <see cref="IGeodesicLine.Position(double, out double, out double, out double, out double, out double, out double, out double)"/>.
        /// </param>
        /// <returns>a <see cref="IGeodesicLine"/> object.</returns>
        /// <remarks>
        /// <para>
        /// This function sets point 3 of the <see cref="IGeodesicLine"/> to correspond to point 2 of the direct geodesic problem.
        /// </para>
        /// <para>
        /// <i>lat1</i> should be in the range [−90°, 90°].
        /// </para>
        /// </remarks>
        IGeodesicLine ArcDirectLine(double lat1, double lon1, double azi1, double a12,
                           GeodesicFlags caps = GeodesicFlags.All);

        #endregion
    }
}