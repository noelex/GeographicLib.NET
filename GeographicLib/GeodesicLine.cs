using static GeographicLib.MathEx;

namespace GeographicLib
{
    /// <summary>
    /// Represents a geodesic line.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <see cref="GeodesicLine"/> facilitates the determination of a series of points on a single geodesic.
    /// The starting point (<i>lat1</i>, <i>lon1</i>) and the azimuth azi1 are specified in the constructor;
    /// alternatively, the <see cref="Geodesic.Line"/> method can be used to create a <see cref="GeodesicLine"/>.
    /// <see cref="IGeodesicLine.Position(double, out double, out double, out double, out double, out double, out double, out double)"/> returns the location of point 2 a distance s12 along the geodesic.
    /// In addition, <see cref="IGeodesicLine.ArcPosition(double, out double, out double, out double, out double, out double, out double, out double, out double)"/> gives the position of point 2 an arc length <i>a12</i> along the geodesic.
    /// </para>
    /// <para>
    /// You can register the position of a reference point 3 a distance (arc length), <i>s13</i> (<i>a13</i>)
    /// along the geodesic with the <see cref="Distance"/> (<see cref="Arc"/>) functions.
    /// Points a fractional distance along the line can be found by providing, for example, 0.5 * <see cref="Distance"/>
    /// as an argument to GeodesicLine.Position. The <see cref="Geodesic.InverseLine"/> or <see cref="GeodesicBase.DirectLine"/> methods return
    /// <see cref="GeodesicLine"/> objects with point 3 set to the point 2 of the corresponding geodesic problem.
    /// <see cref="GeodesicLine"/> objects created with the public constructor or with <see cref="Geodesic.Line"/>
    /// have <i>s13</i> and <i>a13</i> set to <see cref="double.NaN"/>.
    /// </para>
    /// <para>
    /// The default copy constructor and assignment operators work with this class. 
    /// Similarly, a vector can be used to hold <see cref="GeodesicLine"/> objects.
    /// </para>
    /// <para>
    /// The calculations are accurate to better than 15 nm (15 nanometers).
    /// See Sec. 9 of <a href="https://arxiv.org/abs/1102.1215v1">arXiv:1102.1215v1</a> for details.
    /// With <i>exact</i> = <see langword="false"/> (the default) in the constructor for the
    /// <see cref="Geodesic"/> object, the algorithms used by this class are based on series
    /// expansions using the flattening <i>f</i> as a small parameter.These are only
    /// accurate for |<i>f</i>| &lt; 0.02; however reasonably accurate results
    /// will be obtained for |<i>f</i>| &lt; 0.2.  For very eccentric ellipsoids,
    /// set <i>exact</i> = <see langword="true"/> in the constructor for the <see cref="Geodesic"/> object; this will
    /// delegate the calculations to <see cref="GeodesicLineExact"/>.
    /// </para>
    /// <para>
    /// The algorithms are described in
    /// <list type="bullet">
    /// <item>
    /// C. F. F. Karney, <a href="https://doi.org/10.1007/s00190-012-0578-z">Algorithms for geodesics</a>,
    /// J. Geodesy 87, 43–55 (2013); DOI: <a href="https://doi.org/10.1007/s00190-012-0578-z">10.1007/s00190-012-0578-z</a>; 
    /// addenda: <a href="https://geographiclib.sourceforge.io/geod-addenda.html">geod-addenda.html</a>.
    /// </item>
    /// </list>
    /// </para>
    /// <para>
    /// For more information on geodesics see <a href="https://geographiclib.sourceforge.io/html/geodesic.html">Geodesics on an ellipsoid of revolution</a>.
    /// </para>
    /// </remarks>
    public partial class GeodesicLine : GeodesicLineBase
    {
        private Priv _priv;
        private readonly GeodesicLineExact _lineexact;

        internal GeodesicLine(Geodesic g,
                 double lat1, double lon1,
                 double azi1, double salp1, double calp1,
                 GeodesicFlags caps, bool arcmode, double s13_a13)
        {
            if (g.IsExact)
            {
                _lineexact = new GeodesicLineExact(g.GeodesicExact, lat1, lon1, azi1, salp1, calp1, caps, arcmode, s13_a13);
            }
            else
            {
                _priv.Init(g, lat1, lon1, azi1, salp1, calp1, caps);
                SetDistance(arcmode, s13_a13);
            }
        }

        /// <summary>
        /// Constructor for a geodesic line staring at latitude <i>lat1</i>, longitude <i>lon1</i>, 
        /// and azimuth <i>azi1</i> (all in degrees).
        /// </summary>
        /// <param name="g">
        /// A <see cref="Geodesic"/> object used to compute the necessary information about the <see cref="GeodesicLine"/>.
        /// </param>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="caps">bitor'ed combination of <see cref="GeodesicFlags"/> values specifying the capabilities
        /// the <see cref="GeodesicLine"/> object should possess, i.e., which quantities can be returned in calls to
        /// <see cref="IGeodesicLine.Position(double, out double, out double, out double, out double, out double, out double, out double)"/>.</param>
        /// <remarks>
        /// <i>lat1</i> should be in the range [−90°, 90°].
        /// <para>
        /// The <see cref="GeodesicFlags"/> values possible for <i>caps</i> are
        /// <list type="bullet">
        /// <item><i>caps</i> |= <see cref="GeodesicFlags.Latitude"/> for the latitude <i>lat2</i>; this is added automatically;</item>
        /// <item><i>caps</i> |= <see cref="GeodesicFlags.Longitude"/> for the latitude <i>lon2</i>;</item>
        /// <item><i>caps</i> |= <see cref="GeodesicFlags.Azimuth"/> for the latitude <i>azi2</i>;this is added automatically;</item>
        /// <item><i>caps</i> |= <see cref="GeodesicFlags.Distance"/> for the distance <i>s12</i>;</item>
        /// <item><i>caps</i> |= <see cref="GeodesicFlags.ReducedLength"/> for the reduced length <i>m12</i>;</item>
        /// <item><i>caps</i> |= <see cref="GeodesicFlags.GeodesicScale"/> for the geodesic scales <i>M12</i> and <i>M21</i>;</item>
        /// <item><i>caps</i> |= <see cref="GeodesicFlags.Area"/> for the area <i>S12</i>;</item>
        /// <item><i>caps</i> |= <see cref="GeodesicFlags.DistanceIn"/> 
        /// permits the length of the geodesic to be given in terms of <i>s12</i>;
        /// without this capability the length can only be specified in terms of arc length;</item>
        /// <item><i>caps</i> |= <see cref="GeodesicFlags.All"/> for all of the above.</item>
        /// </list>
        /// The default value of <i>caps</i> is <see cref="GeodesicFlags.All"/>.
        /// </para>
        /// <para>
        /// If the point is at a pole, the azimuth is defined by keeping lon1 fixed,
        /// writing <i>lat1</i> = ±(90° − ε), and taking the limit ε → 0+.
        /// </para>
        /// </remarks>
        public GeodesicLine(Geodesic g, double lat1, double lon1, double azi1,
                 GeodesicFlags caps = GeodesicFlags.All)
        {
            if (g.IsExact)
            {
                _lineexact = new GeodesicLineExact(g.GeodesicExact, lat1, lon1, azi1, caps);
            }
            else
            {
                _priv.Init(g, lat1, lon1, azi1, caps);
            }
        }

        /// <inheritdoc/>
        public override void SetDistance(bool arcmode, double s13_a13)
        {
            if (_lineexact != null)
            {
                _lineexact.SetDistance(arcmode, s13_a13);
            }
            else
            {
                if (arcmode)
                {
                    _priv._a13 = s13_a13;
                    // In case the GeodesicLine doesn't have the DISTANCE capability.
                    _priv._s13 = double.NaN;

                    GenPosition(true, _priv._a13, GeodesicFlags.Distance, out _, out _, out _, out _priv._s13, out _, out _, out _, out _);
                }
                else
                {
                    _priv._s13 = s13_a13;

                    // This will set _a13 to NaN if the GeodesicLine doesn't have the
                    // DISTANCE_IN capability.
                    _priv._a13 = GenPosition(false, _priv._s13, 0u, out _, out _, out _, out _, out _, out _, out _, out _);
                }
            }
        }

        /// <inheritdoc/>
        public override double GetDistance(bool arcmode) => _lineexact != null ? _lineexact.GetDistance(arcmode) : (arcmode ? _priv._a13 : _priv._s13);

        /// <inheritdoc/>
        public override double GenPosition(bool arcmode, double s12_a12, GeodesicFlags outmask, out double lat2, out double lon2, out double azi2, out double s12, out double m12, out double M12, out double M21, out double S12)
        {
            if (_lineexact != null)
                return _lineexact.GenPosition(arcmode, s12_a12, outmask,
                                              out lat2, out lon2, out azi2,
                                              out s12, out m12, out M12, out M21, out S12);

            return _priv.GenPosition(arcmode, s12_a12, outmask, out lat2, out lon2, out azi2, out s12, out m12, out M12, out M21, out S12);
        }

        #region Public properties

        /// <inheritdoc/>
        public override double Latitude => _lineexact != null ? _lineexact.Latitude : _priv._lat1;

        /// <inheritdoc/>
        public override double Longitude => _lineexact != null ? _lineexact.Longitude : _priv._lon1;

        /// <inheritdoc/>
        public override double Azimuth => _lineexact != null ? _lineexact.Azimuth : _priv._azi1;

        /// <inheritdoc/>
        public override double SineAzimuth => _lineexact != null ? _lineexact.SineAzimuth : _priv._salp1;

        /// <inheritdoc/>
        public override double CosineAzimuth => _lineexact != null ? _lineexact.CosineAzimuth : _priv._calp1;

        /// <inheritdoc/>
        public override double EquatorialAzimuth => _lineexact != null ? _lineexact.EquatorialAzimuth : Atan2d(_priv._salp0, _priv._calp0);

        /// <inheritdoc/>
        public override double SineEquatorialAzimuth => _lineexact != null ? _lineexact.SineEquatorialAzimuth : _priv._salp0;

        /// <inheritdoc/>
        public override double CosineEquatorialAzimuth => _lineexact != null ? _lineexact.CosineEquatorialAzimuth : _priv._calp0;

        /// <inheritdoc/>
        public override double EquatorialArc => _lineexact != null ? _lineexact.EquatorialArc : Atan2d(_priv._ssig1, _priv._csig1);

        /// <inheritdoc/>
        public override double EquatorialRadius => _lineexact != null ? _lineexact.EquatorialRadius : _priv._a;

        /// <inheritdoc/>
        public override double Flattening => _lineexact != null ? _lineexact.Flattening : _priv._f;

        /// <inheritdoc/>
        public override GeodesicFlags Capabilities => _lineexact != null ? _lineexact.Capabilities : _priv._caps;

        /// <inheritdoc/>
        public override double Distance
        {
            get => GetDistance(false);
            set => SetDistance(false, value);
        }

        /// <inheritdoc/>
        public override double Arc
        {
            get => GetDistance(true);
            set => SetDistance(true, value);
        }

        #endregion
    }
}
