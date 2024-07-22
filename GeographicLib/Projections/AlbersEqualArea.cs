using static GeographicLib.Macros;
using static GeographicLib.MathEx;
using static System.Math;

namespace GeographicLib.Projections
{
    /// <summary>
    /// Albers equal area conic projection.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Implementation taken from the report,
    /// - J.P.Snyder,
    ///   <a href = "https://pubs.usgs.gov/publication/pp1395" > Map Projections: A
    /// Working Manual</a>, USGS Professional Paper 1395 (1987),
    /// pp. 101--102.
    /// </para>
    /// <para>
    /// This is a implementation of the equations in Snyder except that divided
    /// differences will be[have been] used to transform the expressions into
    /// ones which may be evaluated accurately.  [In this implementation, the
    /// projection correctly becomes the cylindrical equal area or the azimuthal
    /// equal area projection when the standard latitude is the equator or a
    /// pole.]
    /// </para>
    /// </remarks>
    public class AlbersEqualArea : IEllipsoid
    {
        private const double eps_ = DBL_EPSILON;
        private static readonly double
            lg2eps_ = -Log2(DBL_EPSILON / 2),
            epsx_ = Sq(eps_),
            epsx2_ = Sq(epsx_),
            tol_ = Sqrt(eps_),
            tol0_ = tol_ * Sqrt(Sqrt(eps_));

        private readonly double _a, _f, _fm, _e2, _e, _e2m, _qZ, _qx;
        private readonly double _sign, _lat0;
        private readonly double _n0, _m02, _nrho0, _txi0, _scxi0, _sxi0;

        private double _k0, _k2;

        /// <summary>
        /// Newton iterations in Reverse
        /// </summary>
        private const int _numit = 5;

        /// <summary>
        /// Newton iterations in Init
        /// </summary>
        private const int _numit0 = 20;

        private AlbersEqualArea(double a, double f)
        {
            _a = a;
            _f = f;
            _fm = 1 - f;
            _e2 = _f * (2 - _f);
            _e = Sqrt(Abs(_e2));
            _e2m = 1 - _e2;
            _qZ = 1 + _e2m * AtanhEE(1d);
            _qx = _qZ / (2 * _e2m);
        }

        /// <summary>
        /// Constructor with a single standard parallel.
        /// </summary>
        /// <param name="ellipsoid"><see cref="IEllipsoid"/> instance to be used in projection.</param>
        /// <param name="stdlat">standard parallel (degrees), the circle of tangency.</param>
        /// <param name="k0">azimuthal scale on the standard parallel.</param>
        public AlbersEqualArea(IEllipsoid ellipsoid, double stdlat, double k0)
            : this(ellipsoid.EquatorialRadius, ellipsoid.Flattening) { }

        /// <summary>
        /// Constructor with two standard parallels.
        /// </summary>
        /// <param name="ellipsoid"><see cref="IEllipsoid"/> instance to be used in projection.</param>
        /// <param name="stdlat1">first standard parallel (degrees).</param>
        /// <param name="stdlat2">second standard parallel (degrees).</param>
        /// <param name="k1">azimuthal scale on the standard parallels.</param>
        public AlbersEqualArea(IEllipsoid ellipsoid, double stdlat1, double stdlat2, double k1)
            : this(ellipsoid.EquatorialRadius, ellipsoid.Flattening, stdlat1, stdlat2, k1) { }

        /// <summary>
        /// Constructor with a single standard parallel.
        /// </summary>
        /// <param name="a">equatorial radius of ellipsoid (meters).</param>
        /// <param name="f">flattening of ellipsoid.  Setting <i>f</i> = 0 gives a sphere.</param>
        /// <param name="stdlat">standard parallel (degrees), the circle of tangency.</param>
        /// <param name="k0">azimuthal scale on the standard parallel.</param>
        public AlbersEqualArea(double a, double f, double stdlat, double k0)
            : this(a, f)
        {
            if (!(IsFinite(_a) && _a > 0))
                throw new GeographicException("Equatorial radius is not positive");
            if (!(IsFinite(_f) && _f < 1))
                throw new GeographicException("Polar semi-axis is not positive");
            if (!(IsFinite(k0) && k0 > 0))
                throw new GeographicException("Scale is not positive");
            if (!(Abs(stdlat) <= QD))
                throw new GeographicException($"Standard latitude not in [-{QD}d, {QD}d]");

            SinCosd(stdlat, out var sphi, out var cphi);
            Init(sphi, cphi, sphi, cphi, k0,
                ref _sign, ref _txi0, ref _scxi0, ref _sxi0, ref _n0, ref _m02, ref _nrho0, ref _lat0);
        }

        /// <summary>
        /// Constructor with two standard parallels.
        /// </summary>
        /// <param name="a">equatorial radius of ellipsoid (meters).</param>
        /// <param name="f">flattening of ellipsoid.  Setting <i>f</i> = 0 gives a sphere.</param>
        /// <param name="stdlat1">first standard parallel (degrees).</param>
        /// <param name="stdlat2">second standard parallel (degrees).</param>
        /// <param name="k1">azimuthal scale on the standard parallels.</param>
        public AlbersEqualArea(double a, double f, double stdlat1, double stdlat2, double k1)
            : this(a, f)
        {
            if (!(IsFinite(_a) && _a > 0))
                throw new GeographicException("Equatorial radius is not positive");
            if (!(IsFinite(_f) && _f < 1))
                throw new GeographicException("Polar semi-axis is not positive");
            if (!(IsFinite(k1) && k1 > 0))
                throw new GeographicException("Scale is not positive");
            if (!(Abs(stdlat1) <= QD))
                throw new GeographicException($"Standard latitude 1 not in [-{QD}d, {QD}d]");
            if (!(Abs(stdlat2) <= QD))
                throw new GeographicException($"Standard latitude 2 not in [-{QD}d, {QD}d]");

            SinCosd(stdlat1, out var sphi1, out var cphi1);
            SinCosd(stdlat2, out var sphi2, out var cphi2);

            Init(sphi1, cphi1, sphi2, cphi2, k1,
                ref _sign, ref _txi0, ref _scxi0, ref _sxi0, ref _n0, ref _m02, ref _nrho0, ref _lat0);
        }

        /// <summary>
        /// Constructor with two standard parallels specified by sines and cosines.
        /// </summary>
        /// <param name="a">equatorial radius of ellipsoid (meters).</param>
        /// <param name="f">flattening of ellipsoid.  Setting <i>f</i> = 0 gives a sphere.</param>
        /// <param name="sinlat1">sine of first standard parallel.</param>
        /// <param name="coslat1">cosine of first standard parallel.</param>
        /// <param name="sinlat2">sine of second standard parallel.</param>
        /// <param name="coslat2">cosine of second standard parallel.</param>
        /// <param name="k1">azimuthal scale on the standard parallels.</param>
        /// <remarks>
        /// This allows parallels close to the poles to be specified accurately.
        /// This routine computes the latitude of origin and the azimuthal scale at
        /// this latitude.If <i>dlat</i> = abs(<i>lat2</i> - <i>lat1</i>) &lt;= 160°,
        /// then the error in the latitude of origin is less than 4.5 * 10^-14 d;.
        /// </remarks>
        public AlbersEqualArea(double a, double f,
            double sinlat1, double coslat1,
            double sinlat2, double coslat2,
            double k1) : this(a, f)
        {
            if (!IsFinite(_a) && _a > 0)
                throw new GeographicException("Equatorial radius is not positive");
            if (!IsFinite(_f) && _f < 1)
                throw new GeographicException("Polar semi-axis is not positive");
            if (!IsFinite(k1) && k1 > 0)
                throw new GeographicException("Scale is not positive");
            if (SignBit(coslat1))
                throw new GeographicException($"Standard latitude 1 not in [-{QD}d, {QD}d]");
            if (SignBit(coslat2))
                throw new GeographicException($"Standard latitude 2 not in [-{QD}d, {QD}d]");
            if (!(Abs(sinlat1) <= 1 && coslat1 <= 1) || (coslat1 == 0 && sinlat1 == 0))
                throw new GeographicException("Bad sine/cosine of standard latitude 1");
            if (!(Abs(sinlat2) <= 1 && coslat2 <= 1) || (coslat2 == 0 && sinlat2 == 0))
                throw new GeographicException("Bad sine/cosine of standard latitude 2");
            if (coslat1 == 0 && coslat2 == 0 && sinlat1 * sinlat2 <= 0)
                throw new GeographicException("Standard latitudes cannot be opposite poles");

            Init(sinlat1, coslat1, sinlat2, coslat2, k1,
                ref _sign, ref _txi0, ref _scxi0, ref _sxi0, ref _n0, ref _m02, ref _nrho0, ref _lat0);
        }

        /// <summary>
        /// Gets a value representing the equatorial radius of the ellipsoid (meters).  This is the value used in the constructor.
        /// </summary>
        public double EquatorialRadius => _a;

        /// <summary>
        /// Gets a value representing the flattening of the ellipsoid.  This is the value used in the constructor.
        /// </summary>
        public double Flattening => _f;

        /// <summary>
        /// Gets a value representing the latitude of the origin for the projection (degrees).
        /// </summary>
        /// <remarks>
        /// This is the latitude of minimum azimuthal scale and equals the <i>stdlat</i> 
        /// in the 1-parallel constructor and lies between <i>stdlat1</i>  and <i>stdlat2</i> 
        /// in the 2-parallel constructors.
        /// </remarks>
        public double OriginLatitude => _lat0;

        /// <summary>
        /// Gets a value representing central scale for the projection.  This is the azimuthal scale on the latitude of origin.
        /// </summary>
        public double CentralScale => _k0;

        /// <summary>
        /// Gets or sets an value representing that whether current <see cref="AlbersEqualArea"/> is frozen.
        /// </summary>
        public bool IsFrozen { get; protected set; }

        /// <summary>
        /// A global instantiation of <see cref="AlbersEqualArea"/> with the WGS84 ellipsoid, 
        /// <i>stdlat</i> = 0, and <i>k0</i> = 1.This degenerates to the cylindrical equalarea projection.
        /// </summary>
        /// <returns></returns>
        public static AlbersEqualArea Cylindrical { get; } = new AlbersEqualArea(Constants.WGS84_a, Constants.WGS84_f, 0, 1, 0, 1, 1).Freeze();

        /// <summary>
        /// A global instantiation of <see cref="AlbersEqualArea"/> with the WGS84 ellipsoid, 
        /// <i>stdlat</i> = 90°, and <i>k0</i> = 1.This degenerates to the Lambert azimuthal equal area projection.
        /// </summary>
        /// <returns></returns>
        public static AlbersEqualArea North { get; } = new AlbersEqualArea(Constants.WGS84_a, Constants.WGS84_f, 1, 0, 1, 0, 1).Freeze();

        /// <summary>
        /// A global instantiation of <see cref="AlbersEqualArea"/> with the WGS84 ellipsoid, 
        /// <i>stdlat</i> = -90°, and <i>k0</i> = 1.This degenerates to the Lambert azimuthal equal area projection.
        /// </summary>
        /// <returns></returns>
        public static AlbersEqualArea South { get; } = new AlbersEqualArea(Constants.WGS84_a, Constants.WGS84_f, -1, 0, -1, 0, 1).Freeze();

        /// <summary>
        /// Freeze current <see cref="AlbersEqualArea"/> instance to prevent its scale being modified.
        /// </summary>
        public AlbersEqualArea Freeze() { IsFrozen = true; return this; }

        /// <summary>
        /// Set the azimuthal scale for the projection.
        /// </summary>
        /// <param name="lat">lat (degrees).</param>
        /// <param name="k">azimuthal scale at latitude <paramref name="lat"/> (default 1).</param>
        /// <remarks>
        /// This allows a "latitude of conformality" to be specified.
        /// </remarks>
        public void SetScale(double lat, double k = 1d)
        {
            if (IsFrozen)
                throw new GeographicException("Projection is frozen");
            if (!(IsFinite(k) && k > 0))
                throw new GeographicException("Scale is not positive");
            if (!(Abs(lat) < QD))
                throw new GeographicException($"Latitude for SetScale not in (-{QD}d, {QD}d)");

            Forward(0, lat, 0, out _, out var kold);
            k /= kold;
            _k0 *= k;
            _k2 = Sq(_k0);
        }

        /// <summary>
        /// Forward projection, from geographic to Lambert conformal conic.
        /// </summary>
        /// <param name="lon0">central meridian longitude (degrees).</param>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <param name="gamma">meridian convergence at point (degrees).</param>
        /// <param name="k">azimuthal scale of projection at point; the radial scale is the 1/<paramref name="k"/>.</param>
        /// <returns>
        /// <i>x</i>, easting of point and <i>y</i>, northing of point, in meters.
        /// </returns>
        /// <remarks>
        /// The latitude origin is given by <see cref="OriginLatitude"/>.  No
        /// false easting or northing is added and <paramref name="lat"/> should be in the range
        /// [-90°, 90°].  The values of <i>x</i> and <i>y</i> returned for
        /// points which project to infinity (i.e., one or both of the poles) will
        /// be large but finite.
        /// </remarks>
        public (double x, double y) Forward(double lon0, double lat, double lon, out double gamma, out double k)
        {
            lon = AngDiff(lon0, lon);
            lat *= _sign;
            double sphi, cphi;
            SinCosd(LatFix(lat) * _sign, out sphi, out cphi);
            cphi = Max(epsx_, cphi);
            double
              lam = lon * Degree,
              tphi = sphi / cphi, txi = Txif(tphi), sxi = txi / Hypot(txi),
              dq = _qZ * Dsn(txi, _txi0, sxi, _sxi0) * (txi - _txi0),
              drho = -_a * dq / (Sqrt(_m02 - _n0 * dq) + _nrho0 / _a),
              theta = _k2 * _n0 * lam, stheta = Sin(theta), ctheta = Cos(theta),
              t = _nrho0 + _n0 * drho;
            var x = t * (_n0 != 0 ? stheta / _n0 : _k2 * lam) / _k0;
            var y = (_nrho0 *
                 (_n0 != 0 ?
                  (ctheta < 0 ? 1 - ctheta : Sq(stheta) / (1 + ctheta)) / _n0 :
                  0)
                 - drho * ctheta) / _k0;
            k = _k0 * (t != 0 ? t * Hypot(_fm * tphi) / _a : 1);
            y *= _sign;
            gamma = _sign * theta / Degree;

            return (x, y);
        }

        /// <summary>
        /// Forward without returning convergence and scale.
        /// </summary>
        /// <param name="lon0">central meridian longitude (degrees).</param>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <returns>
        /// <i>x</i>, easting of point and <i>y</i>, northing of point, in meters.
        /// </returns>
        public (double x, double y) Forward(double lon0, double lat, double lon)
            => Forward(lon0, lat, lon, out _, out _);

        /// <summary>
        /// Reverse projection, from Lambert conformal conic to geographic.
        /// </summary>
        /// <param name="lon0">central meridian longitude (degrees).</param>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <param name="gamma">meridian convergence at point (degrees).</param>
        /// <param name="k">azimuthal scale of projection at point; the radial scale is the 1/<paramref name="k"/>.</param>
        /// <returns>
        /// <i>lat</i>, latitude of point and <i>lon</i>, longitude of point, in degress.
        /// </returns>
        /// <remarks>
        /// The latitude origin is given by <see cref="OriginLatitude"/>.  No
        /// false easting or northing is added. The value of <i>lon</i> returned is in
        /// the range[-180°, 180°].  The value of <i>lat</i> returned is
        /// in the range[-90°, 90°].  If the input point is outside
        /// the legal projected space the nearest pole is returned.
        /// </remarks>
        public (double lat, double lon) Reverse(double lon0, double x, double y, out double gamma, out double k)
        {
            y *= _sign;
            double
              nx = _k0 * _n0 * x, ny = _k0 * _n0 * y, y1 = _nrho0 - ny,
              den = Hypot(nx, y1) + _nrho0, // 0 implies origin with polar aspect
              drho = den != 0 ? (_k0 * x * nx - 2 * _k0 * y * _nrho0 + _k0 * y * ny) / den : 0,
              // dsxia = scxi0 * dsxi
              dsxia = -_scxi0 * (2 * _nrho0 + _n0 * drho) * drho /
                      (Sq(_a) * _qZ),
              txi = (_txi0 + dsxia) / Sqrt(Max(1 - dsxia * (2 * _txi0 + dsxia), epsx2_)),
              tphi = Tphif(txi),
              theta = Atan2(nx, y1),
              lam = _n0 != 0 ? theta / (_k2 * _n0) : x / (y1 * _k0);

            gamma = _sign * theta / Degree;
            var lat = Atand(_sign * tphi);
            var lon = lam / Degree;
            lon = AngNormalize(lon + AngNormalize(lon0));
            k = _k0 * (den != 0 ? (_nrho0 + _n0 * drho) * Hypot(_fm * tphi) / _a : 1);

            return (lat, lon);
        }

        /// <summary>
        /// Reverse without returning convergence and scale.
        /// </summary>
        /// <param name="lon0">central meridian longitude (degrees).</param>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <returns>
        /// <i>lat</i>, latitude of point and <i>lon</i>, longitude of point, in degress.
        /// </returns>
        public (double lat, double lon) Reverse(double lon0, double x, double y)
            => Reverse(lon0, x, y, out _, out _);

        /// <summary>
        /// <para>Returns atanh(      e   * x)/      e   when f > 0,</para>
        /// <para>Returns atan (sqrt(-e2) * x)/sqrt(-e2) when f &lt; 0,</para>
        /// <para>Returns x                              when f = 0.</para>
        /// </summary>
        /// <param name="x"></param>
        private double AtanhEE(double x) =>
            _f > 0 ? Atanh(_e * x) / _e : (_f < 0 ? (Atan(_e * x) / _e) : x);

        /// <summary>
        ///   return atanh(sqrt(x))/sqrt(x) - 1 = y/3 + y^2/5 + y^3/7 + ...
        ///   typical x &lt; e^2 = 2*f
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private double AtanhXm1(double x)
        {
            var s = 0d;
            if (Abs(x) < 0.5)
            {
                Frexp(x, out var e);
                e = -e;
                // x = [0.5,1) * 2^(-e)
                // estimate n s.t. x^n/(2*n+1) < x/3 * epsilon/2
                // a stronger condition is x^(n-1) < epsilon/2
                // taking log2 of both sides, a stronger condition is
                // (n-1)*(-e) < -lg2eps or (n-1)*e > lg2eps or n > ceiling(lg2eps/e)+1
                int n = x == 0 ? 1 : (int)Ceiling(lg2eps_ / e) + 1;
                while (n-- != 0)               // iterating from n-1 down to 0
                    s = x * s + (n != 0 ? 1 : 0) / (double)(2 * n + 1);
            }
            else
            {
                var xs = Sqrt(Abs(x));
                s = (x > 0 ? Atanh(xs) : Atan(xs)) / xs - 1;
            }
            return s;
        }

        /// <summary>
        /// Datanhee(x,y) = atanhee((x-y)/(1-e^2*x*y))/(x-y)
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        private double DAtanhEE(double x, double y)
        {
            double t = x - y, d = 1 - _e2 * x * y;
            return t == 0 ? 1 / d :
              (x * y < 0 ? AtanhEE(x) - AtanhEE(y) : AtanhEE(t / d)) / t;
        }

        /// <summary>
        /// Returns (Datanhee(1,y) - Datanhee(1,x))/(y-x)
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        private double DDAtanhEE(double x, double y)
        {
            // This function is called with x = sphi1, y = sphi2, phi1 <= phi2, sphi2
            // >= 0, abs(sphi1) <= phi2.  However for safety's sake we enforce x <= y.
            if (y < x) Swap(ref x, ref y);      // ensure that x <= y
            double q1 = Abs(_e2),
              q2 = Abs(2 * _e / _e2m * (1 - x));
            return
              x <= 0 || !(Min(q1, q2) < 0.75) ? DDAtanhEE0(x, y) :
              (q1 < q2 ? DDAtanhEE1(x, y) : DDAtanhEE2(x, y));
        }

        // Rearrange difference so that 1 - x is in the denominator, then do a
        // straight divided difference.
        private double DDAtanhEE0(double x, double y)
            => (DAtanhEE(1, y) - DAtanhEE(x, y)) / (1 - x);

        // The expansion for e2 small
        private double DDAtanhEE1(double x, double y)
        {
            // The series in e2 is
            //   sum( c[l] * e2^l, l, 1, N)
            // where
            //   c[l] = sum( x^i * y^j; i >= 0, j >= 0, i+j < 2*l) / (2*l + 1)
            //        = ( (x-y) - (1-y) * x^(2*l+1) + (1-x) * y^(2*l+1) ) /
            //          ( (2*l+1) * (x-y) * (1-y) * (1-x) )
            // For x = y = 1, c[l] = l
            //
            // In the limit x,y -> 1,
            //
            //   DDatanhee -> e2/(1-e2)^2 = sum(l * e2^l, l, 1, inf)
            //
            // Use if e2 is sufficiently small.
            var s = 0d;
            double z = 1, k = 1, t = 0, c = 0, en = 1;
            while (true)
            {
                t = y * t + z; c += t; z *= x;
                t = y * t + z; c += t; z *= x;
                k += 2; en *= _e2;
                // Here en[l] = e2^l, k[l] = 2*l + 1,
                // c[l] = sum( x^i * y^j; i >= 0, j >= 0, i+j < 2*l) / (2*l + 1)
                // Taylor expansion is
                // s = sum( c[l] * e2^l, l, 1, N)
                var ds = en * c / k;
                s += ds;
                if (!(Abs(ds) > Abs(s) * eps_ / 2))
                    break;            // Iterate until the added term is sufficiently small
            }
            return s;
        }

        // The expansion for x (and y) close to 1
        private double DDAtanhEE2(double x, double y)
        {
            // If x and y are both close to 1, expand in Taylor series in dx = 1-x and
            // dy = 1-y:
            //
            // DDatanhee = sum(C_m * (dx^(m+1) - dy^(m+1)) / (dx - dy), m, 0, inf)
            //
            // where
            //
            // C_m = sum( (m+2)!! / (m+2-2*k)!! *
            //            ((m+1)/2)! / ((m+1)/2-k)! /
            //            (k! * (2*k-1)!!) *
            //            e2^((m+1)/2+k),
            //           k, 0, (m+1)/2) * (-1)^m / ((m+2) * (1-e2)^(m+2))
            // for m odd, and
            //
            // C_m = sum( 2 * (m+1)!! / (m+1-2*k)!! *
            //            (m/2+1)! / (m/2-k)! /
            //            (k! * (2*k+1)!!) *
            //            e2^(m/2+1+k),
            //           k, 0, m/2)) * (-1)^m / ((m+2) * (1-e2)^(m+2))
            // for m even.
            //
            // Here i!! is the double factorial extended to negative i with
            // i!! = (i+2)!!/(i+2).
            //
            // Note that
            //   (dx^(m+1) - dy^(m+1)) / (dx - dy) =
            //     dx^m + dx^(m-1)*dy ... + dx*dy^(m-1) + dy^m
            //
            // Leading (m = 0) term is e2 / (1 - e2)^2
            //
            // Magnitude of mth term relative to the leading term scales as
            //
            //   2*(2*e/(1-e2)*dx)^m
            //
            // So use series if (2*e/(1-e2)*dx) is sufficiently small
            double s, dx = 1 - x, dy = 1 - y, xy = 1, yy = 1, ee = _e2 / Sq(_e2m);
            s = ee;
            for (int m = 1; ; ++m)
            {
                double c = m + 2, t = c;
                yy *= dy;               // yy = dy^m
                xy = dx * xy + yy;
                // Now xy = dx^m + dx^(m-1)*dy ... + dx*dy^(m-1) + dy^m
                //        = (dx^(m+1) - dy^(m+1)) / (dx - dy)
                // max value = (m+1) * max(dx,dy)^m
                ee /= -_e2m;
                if (m % 2 == 0) ee *= _e2;
                // Now ee = (-1)^m * e2^(floor(m/2)+1) / (1-e2)^(m+2)
                int kmax = (m + 1) / 2;
                for (int k = kmax - 1; k >= 0; --k)
                {
                    // max coeff is less than 2^(m+1)
                    c *= (k + 1) * (2 * (k + m - 2 * kmax) + 3);
                    c /= (kmax - k) * (2 * (kmax - k) + 1);
                    // Horner sum for inner _e2 series
                    t = _e2 * t + c;
                }
                // Straight sum for outer m series
                var ds = t * ee * xy / (m + 2);
                s = s + ds;
                if (!(Abs(ds) > Abs(s) * eps_ / 2))
                    break;            // Iterate until the added term is sufficiently small
            }
            return s;
        }

        /// <summary>
        /// Divided differences
        /// <para>Definition: Df(x,y) = (f(x)-f(y))/(x-y)</para>
        /// <para>
        /// See:
        ///   W. M. Kahan and R. J. Fateman,
        ///   Symbolic computation of divided differences,
        ///   SIGSAM Bull. 33(2), 7-28 (1999)
        ///   https://doi.org/10.1145/334714.334716
        ///   http://www.cs.berkeley.edu/~fateman/papers/divdiff.pdf
        /// </para>
        /// <para>General rules</para>
        /// <para>h(x) = f(g(x)): Dh(x,y) = Df(g(x),g(y))*Dg(x,y)</para>
        /// <para>h(x) = f(x)*g(x):</para>
        /// <para>Dh(x,y) = Df(x,y)*g(x) + Dg(x,y)*f(y)</para>
        /// <para>        = Df(x,y)*g(y) + Dg(x,y)*f(x)</para>
        /// <para>        = Df(x,y)*(g(x)+g(y))/2 + Dg(x,y)*(f(x)+f(y))/2</para>
        ///
        /// <para>sn(x) = x/sqrt(1+x^2): Dsn(x,y) = (x+y)/((sn(x)+sn(y))*(1+x^2)*(1+y^2))</para>
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="sx"></param>
        /// <param name="sy"></param>
        /// <returns></returns>
        private double Dsn(double x, double y, double sx, double sy)
        {
            // sx = x/hyp(x)
            var t = x * y;
            return t > 0 ? (x + y) * Sq((sx * sy) / t) / (sx + sy) :
              (x - y != 0 ? (sx - sy) / (x - y) : 1);
        }

        internal double Txif(double tphi)
        {
            // sxi = ( sphi/(1-e2*sphi^2) + atanhee(sphi) ) /
            //       ( 1/(1-e2) + atanhee(1) )
            //
            // txi = ( sphi/(1-e2*sphi^2) + atanhee(sphi) ) /
            //       sqrt( ( (1+e2*sphi)*(1-sphi)/( (1-e2*sphi^2) * (1-e2) ) +
            //               atanhee((1-sphi)/(1-e2*sphi)) ) *
            //             ( (1-e2*sphi)*(1+sphi)/( (1-e2*sphi^2) * (1-e2) ) +
            //               atanhee((1+sphi)/(1+e2*sphi)) ) )
            //     = ( tphi/(1-e2*sphi^2) + atanhee(sphi, e2)/cphi ) /
            //       sqrt(
            //       ( (1+e2*sphi)/( (1-e2*sphi^2) * (1-e2) ) + Datanhee(1,  sphi)  ) *
            //       ( (1-e2*sphi)/( (1-e2*sphi^2) * (1-e2) ) + Datanhee(1, -sphi)  ) )
            //
            // This function maintains odd parity
            double
              cphi = 1 / Sqrt(1 + Sq(tphi)),
              sphi = tphi * cphi,
              es1 = _e2 * sphi,
              es2m1 = 1 - es1 * sphi,   // 1 - e2 * sphi^2
              es2m1a = _e2m * es2m1;    // (1 - e2 * sphi^2) * (1 - e2)
            return (tphi / es2m1 + AtanhEE(sphi) / cphi) /
              Sqrt(((1 + es1) / es2m1a + DAtanhEE(1, sphi)) *
                    ((1 - es1) / es2m1a + DAtanhEE(1, -sphi)));
        }

        internal double Tphif(double txi)
        {
            double
              tphi = txi,
              stol = tol_ * Max(1d, Abs(txi));

            // CHECK: min iterations = 1, max iterations = 2; mean = 1.99
            for (int i = 0; i < _numit || GEOGRAPHICLIB_PANIC; ++i)
            {
                // dtxi/dtphi = (scxi/scphi)^3 * 2*(1-e^2)/(qZ*(1-e^2*sphi^2)^2)
                double
                  txia = Txif(tphi),
                  tphi2 = Sq(tphi),
                  scphi2 = 1 + tphi2,
                  scterm = scphi2 / (1 + Sq(txia)),
                  dtphi = (txi - txia) * scterm * Sqrt(scterm) *
                  _qx * Sq(1 - _e2 * tphi2 / scphi2);

                tphi += dtphi;
                if (!(Abs(dtphi) >= stol))
                    break;
            }
            return tphi;
        }

        private void Init(double sphi1, double cphi1,
                 double sphi2, double cphi2, double k1,
                 ref double _sign, ref double _txi0, ref double _scxi0, ref double _sxi0, ref double _n0, ref double _m02,
                 ref double _nrho0, ref double _lat0)
        {
            {
                double r;
                r = Hypot(sphi1, cphi1);
                sphi1 /= r; cphi1 /= r;
                r = Hypot(sphi2, cphi2);
                sphi2 /= r; cphi2 /= r;
            }

            bool polar = (cphi1 == 0);
            cphi1 = Max(epsx_, cphi1);   // Avoid singularities at poles
            cphi2 = Max(epsx_, cphi2);
            // Determine hemisphere of tangent latitude
            _sign = sphi1 + sphi2 >= 0 ? 1 : -1;
            // Internally work with tangent latitude positive
            sphi1 *= _sign; sphi2 *= _sign;
            if (sphi1 > sphi2)
            {
                Swap(ref sphi1, ref sphi2); Swap(ref cphi1, ref cphi2); // Make phi1 < phi2
            }

            double
              tphi1 = sphi1 / cphi1, tphi2 = sphi2 / cphi2;

            // q = (1-e^2)*(sphi/(1-e^2*sphi^2) - atanhee(sphi))
            // qZ = q(pi/2) = (1 + (1-e^2)*atanhee(1))
            // atanhee(x) = atanh(e*x)/e
            // q = sxi * qZ
            // dq/dphi = 2*(1-e^2)*cphi/(1-e^2*sphi^2)^2
            //
            // n = (m1^2-m2^2)/(q2-q1) -> sin(phi0) for phi1, phi2 -> phi0
            // C = m1^2 + n*q1 = (m1^2*q2-m2^2*q1)/(q2-q1)
            // let
            //   rho(pi/2)/rho(-pi/2) = (1-s)/(1+s)
            //   s = n*qZ/C
            //     = qZ * (m1^2-m2^2)/(m1^2*q2-m2^2*q1)
            //     = qZ * (scbet2^2 - scbet1^2)/(scbet2^2*q2 - scbet1^2*q1)
            //     = (scbet2^2 - scbet1^2)/(scbet2^2*sxi2 - scbet1^2*sxi1)
            //     = (tbet2^2 - tbet1^2)/(scbet2^2*sxi2 - scbet1^2*sxi1)
            // 1-s = -((1-sxi2)*scbet2^2 - (1-sxi1)*scbet1^2)/
            //         (scbet2^2*sxi2 - scbet1^2*sxi1)
            //
            // Define phi0 to give same value of s, i.e.,
            //  s = sphi0 * qZ / (m0^2 + sphi0*q0)
            //    = sphi0 * scbet0^2 / (1/qZ + sphi0 * scbet0^2 * sxi0)

            double tphi0, C;
            if (polar || tphi1 == tphi2)
            {
                tphi0 = tphi2;
                C = 1;                    // ignored
            }
            else
            {
                double
                  tbet1 = _fm * tphi1, scbet12 = 1 + Sq(tbet1),
                  tbet2 = _fm * tphi2, scbet22 = 1 + Sq(tbet2),
                  txi1 = Txif(tphi1), cxi1 = 1 / Hypot(txi1), sxi1 = txi1 * cxi1,
                  txi2 = Txif(tphi2), cxi2 = 1 / Hypot(txi2), sxi2 = txi2 * cxi2,
                  dtbet2 = _fm * (tbet1 + tbet2),
                  es1 = 1 - _e2 * Sq(sphi1), es2 = 1 - _e2 * Sq(sphi2),
                  /*
                  dsxi = ( (_e2 * sq(sphi2 + sphi1) + es2 + es1) / (2 * es2 * es1) +
                           Datanhee(sphi2, sphi1) ) * Dsn(tphi2, tphi1, sphi2, sphi1) /
                  ( 2 * _qx ),
                  */
                  dsxi = ((1 + _e2 * sphi1 * sphi2) / (es2 * es1) +
                           DAtanhEE(sphi2, sphi1)) * Dsn(tphi2, tphi1, sphi2, sphi1) /
                  (2 * _qx),
                  den = (sxi2 + sxi1) * dtbet2 + (scbet22 + scbet12) * dsxi,
                  // s = (sq(tbet2) - sq(tbet1)) / (scbet22*sxi2 - scbet12*sxi1)
                  s = 2 * dtbet2 / den,
                  // 1-s = -(sq(scbet2)*(1-sxi2) - sq(scbet1)*(1-sxi1)) /
                  //        (scbet22*sxi2 - scbet12*sxi1)
                  // Write
                  //   sq(scbet)*(1-sxi) = sq(scbet)*(1-sphi) * (1-sxi)/(1-sphi)
                  sm1 = -Dsn(tphi2, tphi1, sphi2, sphi1) *
                  (-(((sphi2 <= 0 ? (1 - sxi2) / (1 - sphi2) :
                         Sq(cxi2 / cphi2) * (1 + sphi2) / (1 + sxi2)) +
                        (sphi1 <= 0 ? (1 - sxi1) / (1 - sphi1) :
                         Sq(cxi1 / cphi1) * (1 + sphi1) / (1 + sxi1)))) *
                    (1 + _e2 * (sphi1 + sphi2 + sphi1 * sphi2)) /
                    (1 + (sphi1 + sphi2 + sphi1 * sphi2)) +
                    (scbet22 * (sphi2 <= 0 ? 1 - sphi2 :
                                Sq(cphi2) / (1 + sphi2)) +
                     scbet12 * (sphi1 <= 0 ? 1 - sphi1 : Sq(cphi1) / (1 + sphi1)))
                    * (_e2 * (1 + sphi1 + sphi2 + _e2 * sphi1 * sphi2) / (es1 * es2)
                    + _e2m * DDAtanhEE(sphi1, sphi2)) / _qZ) / den;
                // C = (scbet22*sxi2 - scbet12*sxi1) / (scbet22 * scbet12 * (sx2 - sx1))
                C = den / (2 * scbet12 * scbet22 * dsxi);
                tphi0 = (tphi2 + tphi1) / 2;
                double stol = tol0_ * Max(1d, Abs(tphi0));
                for (int i = 0; i < 2 * _numit0 || GEOGRAPHICLIB_PANIC; ++i)
                {
                    // Solve (scbet0^2 * sphi0) / (1/qZ + scbet0^2 * sphi0 * sxi0) = s
                    // for tphi0 by Newton's method on
                    // v(tphi0) = (scbet0^2 * sphi0) - s * (1/qZ + scbet0^2 * sphi0 * sxi0)
                    //          = 0
                    // Alt:
                    // (scbet0^2 * sphi0) / (1/qZ - scbet0^2 * sphi0 * (1-sxi0))
                    //          = s / (1-s)
                    // w(tphi0) = (1-s) * (scbet0^2 * sphi0)
                    //             - s  * (1/qZ - scbet0^2 * sphi0 * (1-sxi0))
                    //          = (1-s) * (scbet0^2 * sphi0)
                    //             - S/qZ  * (1 - scbet0^2 * sphi0 * (qZ-q0))
                    // Now
                    // qZ-q0 = (1+e2*sphi0)*(1-sphi0)/(1-e2*sphi0^2) +
                    //         (1-e2)*atanhee((1-sphi0)/(1-e2*sphi0))
                    // In limit sphi0 -> 1, qZ-q0 -> 2*(1-sphi0)/(1-e2), so wrte
                    // qZ-q0 = 2*(1-sphi0)/(1-e2) + A + B
                    // A = (1-sphi0)*( (1+e2*sphi0)/(1-e2*sphi0^2) - (1+e2)/(1-e2) )
                    //   = -e2 *(1-sphi0)^2 * (2+(1+e2)*sphi0) / ((1-e2)*(1-e2*sphi0^2))
                    // B = (1-e2)*atanhee((1-sphi0)/(1-e2*sphi0)) - (1-sphi0)
                    //   = (1-sphi0)*(1-e2)/(1-e2*sphi0)*
                    //     ((atanhee(x)/x-1) - e2*(1-sphi0)/(1-e2))
                    // x = (1-sphi0)/(1-e2*sphi0), atanhee(x)/x = atanh(e*x)/(e*x)
                    //
                    // 1 - scbet0^2 * sphi0 * (qZ-q0)
                    //   = 1 - scbet0^2 * sphi0 * (2*(1-sphi0)/(1-e2) + A + B)
                    //   = D - scbet0^2 * sphi0 * (A + B)
                    // D = 1 - scbet0^2 * sphi0 * 2*(1-sphi0)/(1-e2)
                    //   = (1-sphi0)*(1-e2*(1+2*sphi0*(1+sphi0)))/((1-e2)*(1+sphi0))
                    // dD/dsphi0 = -2*(1-e2*sphi0^2*(2*sphi0+3))/((1-e2)*(1+sphi0)^2)
                    // d(A+B)/dsphi0 = 2*(1-sphi0^2)*e2*(2-e2*(1+sphi0^2))/
                    //                 ((1-e2)*(1-e2*sphi0^2)^2)

                    double
                      scphi02 = 1 + Sq(tphi0), scphi0 = Sqrt(scphi02),
                      // sphi0m = 1-sin(phi0) = 1/( sec(phi0) * (tan(phi0) + sec(phi0)) )
                      sphi0 = tphi0 / scphi0, sphi0m = 1 / (scphi0 * (tphi0 + scphi0)),
                      // scbet0^2 * sphi0
                      g = (1 + Sq(_fm * tphi0)) * sphi0,
                      // dg/dsphi0 = dg/dtphi0 * scphi0^3
                      dg = _e2m * scphi02 * (1 + 2 * Sq(tphi0)) + _e2,
                      D = sphi0m * (1 - _e2 * (1 + 2 * sphi0 * (1 + sphi0))) / (_e2m * (1 + sphi0)),
                      // dD/dsphi0
                      dD = -2 * (1 - _e2 * Sq(sphi0) * (2 * sphi0 + 3)) /
                           (_e2m * Sq(1 + sphi0)),
                      A = -_e2 * Sq(sphi0m) * (2 + (1 + _e2) * sphi0) /
                          (_e2m * (1 - _e2 * Sq(sphi0))),
                      B = (sphi0m * _e2m / (1 - _e2 * sphi0) *
                           (AtanhXm1(_e2 *
                                     Sq(sphi0m / (1 - _e2 * sphi0))) - _e2 * sphi0m / _e2m)),
                      // d(A+B)/dsphi0
                      dAB = (2 * _e2 * (2 - _e2 * (1 + Sq(sphi0))) /
                             (_e2m * Sq(1 - _e2 * Sq(sphi0)) * scphi02)),
                      u = sm1 * g - s / _qZ * (D - g * (A + B)),
                      // du/dsphi0
                      du = sm1 * dg - s / _qZ * (dD - dg * (A + B) - g * dAB),
                      dtu = -u / du * (scphi0 * scphi02);
                    tphi0 += dtu;
                    if (!(Abs(dtu) >= stol))
                        break;
                }
            }
            _txi0 = Txif(tphi0); _scxi0 = Hypot(_txi0); _sxi0 = _txi0 / _scxi0;
            _n0 = tphi0 / Hypot(tphi0);
            _m02 = 1 / (1 + Sq(_fm * tphi0));
            _nrho0 = polar ? 0 : _a * Sqrt(_m02);
            _k0 = Sqrt(tphi1 == tphi2 ? 1 : C / (_m02 + _n0 * _qZ * _sxi0)) * k1;
            _k2 = Sq(_k0);
            _lat0 = _sign * Atan(tphi0) / Degree;
        }
    }
}
