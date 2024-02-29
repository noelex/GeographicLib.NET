using static GeographicLib.MathEx;
using static System.Math;

namespace GeographicLib
{
    /// <summary>
    /// Represent an exact geodesic line.
    /// </summary>
    /// <remarks>
    /// <see cref="GeodesicLineExact"/> facilitates the determination of a series of points on a
    /// single geodesic.This is a companion to the <see cref="GeodesicExact"/> class.  For
    /// additional information on this class see the documentation on the
    /// <see cref="GeodesicLine"/> class.
    /// </remarks>
    public partial class GeodesicLineExact : GeodesicLineBase
    {
        private Priv _priv;

        internal GeodesicLineExact(GeodesicExact g,
                double lat1, double lon1,
                double azi1, double salp1, double calp1,
                GeodesicFlags caps, bool arcmode, double s13_a13)
            : this(g, lat1, lon1, azi1, salp1, calp1, caps)
        {
            SetDistance(arcmode, s13_a13);
        }

        internal GeodesicLineExact(GeodesicExact g,
                double lat1, double lon1,
                double azi1, double salp1, double calp1,
                GeodesicFlags caps)
        {
            _priv.Init(g, lat1, lon1, azi1, salp1, calp1, caps);
        }

        /// <summary>
        /// Constructor for a exact geodesic line staring at latitude <i>lat1</i>, longitude <i>lon1</i>, 
        /// and azimuth <i>azi1</i> (all in degrees).
        /// </summary>
        /// <param name="g">
        /// A <see cref="GeodesicExact"/> object used to compute the necessary information about the <see cref="GeodesicLineExact"/>.
        /// </param>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="caps">bitor'ed combination of <see cref="GeodesicFlags"/> values specifying the capabilities
        /// the <see cref="GeodesicLineExact"/> object should possess, i.e., which quantities can be returned in calls to
        /// <see cref="GeodesicLineBase.Position(double, out double, out double, out double, out double, out double, out double, out double)"/>.</param>
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
        public GeodesicLineExact(GeodesicExact g,
                                 double lat1, double lon1, double azi1,
                                 GeodesicFlags caps)
        {
            _priv.Init(g, lat1, lon1, azi1, caps);
        }

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

        /// <inheritdoc/>
        public override double Azimuth => _priv._azi1;

        /// <inheritdoc/>
        public override GeodesicFlags Capabilities => _priv._caps;

        /// <inheritdoc/>
        public override double CosineAzimuth => _priv._calp1;

        /// <inheritdoc/>
        public override double CosineEquatorialAzimuth => _priv._calp0;

        /// <inheritdoc/>
        public override double EquatorialArc => Atan2(_priv._ssig1, _priv._csig1) / Degree;

        /// <inheritdoc/>
        public override double EquatorialAzimuth => Atan2d(_priv._salp0, _priv._calp0);

        /// <inheritdoc/>
        public override double Latitude => _priv._lat1;

        /// <inheritdoc/>
        public override double Longitude => _priv._lon1;

        /// <inheritdoc/>
        public override double SineAzimuth => _priv._salp1;

        /// <inheritdoc/>
        public override double SineEquatorialAzimuth => _priv._salp0;

        /// <inheritdoc/>
        public override double EquatorialRadius => _priv._a;

        /// <inheritdoc/>
        public override double Flattening => _priv._f;

        /// <inheritdoc/>
        public override double GenPosition(bool arcmode, double s12_a12, GeodesicFlags outmask, out double lat2,
            out double lon2, out double azi2, out double s12, out double m12, out double M12, out double M21, out double S12)
        {
            return _priv.GenPosition(arcmode, s12_a12, outmask, out lat2, out lon2, out azi2, out s12, out m12, out M12, out M21, out S12);
        }

        /// <inheritdoc/>
        public override double GetDistance(bool arcmode) => (arcmode ? _priv._a13 : _priv._s13);

        /// <inheritdoc/>
        public override void SetDistance(bool arcmode, double s13_a13)
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
}
