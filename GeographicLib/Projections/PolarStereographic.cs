using static GeographicLib.MathEx;
using static System.Math;

namespace GeographicLib.Projections
{
    /// <summary>
    /// Polar stereographic projection.
    /// </summary>
    /// <remarks>
    /// Implementation taken from the report,
    /// <list type="bullet">
    /// <item>
    /// J. P. Snyder, <a href="http://pubs.er.usgs.gov/usgspubs/pp/pp1395">Map Projections: A Working Manual</a>,
    /// USGS Professional Paper 1395 (1987), pp. 160–163.
    /// </item>
    /// </list>
    /// <para>
    /// This is a straightforward implementation of the equations in Snyder except that Newton's method is used to invert the projection.
    /// </para>
    /// <para>
    /// This class also returns the meridian convergence <i>gamma</i> and scale <i>k</i>. 
    /// The meridian convergence is the bearing of grid north (the <i>y</i> axis) measured clockwise from true north.
    /// </para>
    /// </remarks>
    public class PolarStereographic : IEllipsoid
    {
        private readonly double _a, _f, _e2, _es, _e2m, _c;
        private double _k0;

        /// <summary>
        /// Initialize a new <see cref="PolarStereographic"/> instance with specified ellipsoid and scale.
        /// </summary>
        /// <param name="ellipsoid">the ellipsoid.</param>
        /// <param name="k0">central scale factor.</param>
        public PolarStereographic(IEllipsoid ellipsoid, double k0)
            : this(ellipsoid.EquatorialRadius, ellipsoid.Flattening, k0)
        {
        }

        /// <summary>
        /// Initialize a new <see cref="PolarStereographic"/> instance with specified equatorial radius, flattening of ellipsoid and scale.
        /// </summary>
        /// <param name="a">equatorial radius (meters).</param>
        /// <param name="f">flattening of ellipsoid. Setting <i>f</i> = 0 gives a sphere. Negative <i>f</i> gives a prolate ellipsoid.</param>
        /// <param name="k0">central scale factor.</param>
        public PolarStereographic(double a, double f, double k0)
        {
            _a = a;
            _f = f;
            _e2 = _f * (2 - _f);
            _es = (_f < 0 ? -1 : 1) * Sqrt(Abs(_e2));
            _e2m = 1 - _e2;
            _c = (1 - _f) * Exp(EAtanhE(1, _es));
            _k0 = k0;

            if (!(IsFinite(_a) && _a > 0))
                throw new GeographicException("Equatorial radius is not positive");
            if (!(IsFinite(_f) && _f < 1))
                throw new GeographicException("Polar semi-axis is not positive");
            if (!(IsFinite(_k0) && _k0 > 0))
                throw new GeographicException("Scale is not positive");
        }

        /// <summary>
        /// A global instantiation of <see cref="PolarStereographic"/> with
        /// the WGS84 ellipsoid and the UPS scale factor. However, unlike UPS, no false easting or northing is added.
        /// </summary>
        public static PolarStereographic UPS { get; } = new PolarStereographic(Ellipsoid.WGS84, Constants.UPS_k0);

        /// <summary>
        /// Gets a value representing central scale for the projection.  This is the azimuthal scale on the latitude of origin.
        /// </summary>
        public double CentralScale => _k0;

        /// <summary>
        /// Gets a value representing the equatorial radius (<i>a</i>) of the ellipsoid.
        /// </summary>
        public double EquatorialRadius => _a;

        /// <summary>
        /// Gets a value representing the flatterning (<i>f</i>) of the ellipsoid.
        /// </summary>
        public double Flattening => _f;

        /// <summary>
        /// Gets or sets an value representing that whether current <see cref="PolarStereographic"/> is frozen.
        /// </summary>
        public bool IsFrozen { get; private set; }

        /// <summary>
        /// Freeze current <see cref="PolarStereographic"/> instance to prevent its scale being modified.
        /// </summary>
        public PolarStereographic Freeze() { IsFrozen = true; return this; }

        /// <summary>
        /// Forward projection, from geographic to projected coordinate system.
        /// </summary>
        /// <param name="northp">the pole which is the center of projection (<see langword="true"/> means north, <see langword="false"/> means south).</param>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <param name="gamma">meridian convergence at point (degrees).</param>
        /// <param name="k">scale of projection at point.</param>
        /// <returns>
        /// <i>x</i>, easting of point (meters) and <i>y</i>, northing of point (meters).
        /// </returns>
        public (double x, double y) Forward(bool northp, double lat, double lon, out double gamma, out double k)
        {
            double x, y;
            lat = LatFix(lat);
            lat *= northp ? 1 : -1;
            double
              tau = Tand(lat),
              secphi = Hypot(1, tau),
              taup = Taupf(tau, _es),
              rho = Hypot(1, taup) + Abs(taup);
            rho = taup >= 0 ? (lat != qd ? 1 / rho : 0) : rho;
            rho *= 2 * _k0 * _a / _c;
            k = lat != qd ?
                (rho / _a) * secphi * Sqrt(_e2m + _e2 / Sq(secphi)) :
              _k0;
            SinCosd(lon, out x, out y);
            x *= rho;
            y *= (northp ? -rho : rho);
            gamma = AngNormalize(northp ? lon : -lon);

            return (x, y);
        }

        /// <summary>
        /// Reverse projection, from projected coordinate system to geographic.
        /// </summary>
        /// <param name="northp">the pole which is the center of projection (<see langword="true"/> means north, <see langword="false"/> means south).</param>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <param name="gamma">meridian convergence at point (degrees).</param>
        /// <param name="k">azimuthal scale of projection at point; the radial scale is the 1/<paramref name="k"/>.</param>
        /// <returns>
        /// <i>lat</i>, latitude of point (degrees) and <i>lon</i>, longitude of point (degrees).
        /// </returns>
        public (double lat, double lon) Reverse(bool northp, double x, double y, out double gamma, out double k)
        {
            double
              lat, lon,
              rho = Hypot(x, y),
              t = rho != 0 ? rho / (2 * _k0 * _a / _c) :
              Sq(DBL_EPSILON),
              taup = (1 / t - t) / 2,
              tau = Tauf(taup, _es),
              secphi = Hypot(1, tau);
            k = rho != 0 ? (rho / _a) * secphi * Sqrt(_e2m + _e2 / Sq(secphi)) :
              _k0;
            lat = (northp ? 1 : -1) * Atand(tau);
            lon = Atan2d(x, northp ? -y : y);
            gamma = AngNormalize(northp ? lon : -lon);

            return (lat, lon);
        }

        /// <summary>
        /// Forward without returning convergence and scale.
        /// </summary>
        /// <param name="northp">the pole which is the center of projection (<see langword="true"/> means north, <see langword="false"/> means south).</param>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <returns>
        /// <i>x</i>, easting of point (meters) and <i>y</i>, northing of point (meters).
        /// </returns>
        public (double x, double y) Forward(bool northp, double lat, double lon) => Forward(northp, lat, lon, out _, out _);

        /// <summary>
        /// Reverse without returning convergence and scale.
        /// </summary>
        /// <param name="northp">the pole which is the center of projection (<see langword="true"/> means north, <see langword="false"/> means south).</param>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <returns>
        /// <i>lat</i>, latitude of point (degrees) and <i>lon</i>, longitude of point (degrees).
        /// </returns>
        public (double lat, double lon) Reverse(bool northp, double x, double y) => Reverse(northp, x, y, out _, out _);

        /// <summary>
        /// Set the scale for the projection.
        /// </summary>
        /// <param name="lat">(degrees).</param>
        /// <param name="k">scale at latitude <paramref name="lat"/> (default 1).</param>
        public void SetScale(double lat, double k = 1)
        {
            if (IsFrozen)
                throw new GeographicException("Projection is frozen");
            if (!(IsFinite(k) && k > 0))
                throw new GeographicException("Scale is not positive");
            if (!(-qd < lat && lat <= qd))
                throw new GeographicException($"Latitude must be in (-{qd}d, {qd}d]");

            _k0 = 1;
            Forward(true, lat, 0, out _, out var kold);
            _k0 *= k / kold;
        }
    }
}
