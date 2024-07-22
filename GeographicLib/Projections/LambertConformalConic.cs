using static GeographicLib.MathEx;
using static System.Math;

namespace GeographicLib.Projections
{
    /// <summary>
    /// Lambert conformal conic projection.
    /// </summary>
    /// <remarks>
    /// Implementation taken from the report,
    /// <list type="bullet">
    /// <item>
    /// J. P. Snyder, <a href="https://pubs.usgs.gov/publication/pp1395">Map Projections: A Working Manual</a>,
    /// USGS Professional Paper 1395 (1987), pp. 107–109.
    /// </item>
    /// </list>
    /// <para>
    /// This is a implementation of the equations in Snyder except that divided differences have been used to 
    /// transform the expressions into ones which may be evaluated accurately and that Newton's method is used
    /// to invert the projection. In this implementation, the projection correctly becomes the Mercator projection
    /// or the polar stereographic projection when the standard latitude is the equator or a pole. The accuracy of
    /// the projections is about 10 nm (10 nanometers).
    /// </para>
    /// <para>
    /// The ellipsoid parameters, the standard parallels, and the scale on the standard parallels are set in the constructor.
    /// Internally, the case with two standard parallels is converted into a single standard parallel, the latitude of tangency
    /// (also the latitude of minimum scale), with a scale specified on this parallel. This latitude is also used as the latitude
    /// of origin which is returned by <see cref="OriginLatitude"/>. The scale on the latitude of origin is given by 
    /// <see cref="CentralScale"/>. The case with two distinct standard parallels where one is a pole is singular and is disallowed.
    /// The central meridian (which is a trivial shift of the longitude) is specified as the <i>lon0</i> argument of the 
    /// <see cref="Forward(double, double, double, out double, out double)"/>
    /// and <see cref="Reverse(double, double, double, out double, out double)"/> functions.
    /// </para>
    /// <para>
    /// This class also returns the meridian convergence <i>gamma</i> and scale <i>k</i>. 
    /// The meridian convergence is the bearing of grid north (the <i>y</i> axis) measured clockwise from true north.
    /// </para>
    /// </remarks>
    public class LambertConformalConic : IEllipsoid
    {
        private const int numit_ = 5;

        private readonly double eps_, epsx_, ahypover_;
        private readonly double _a, _f, _fm, _e2, _es;
        private readonly double _sign, _n, _nc, _t0nm1, _lat0;
        private readonly double _scbet0, _tchi0, _scchi0, _psi0, _nrho0, _drhomax;

        private double _scale, _k0;

        #region Constructors

        private LambertConformalConic(double a, double f)
        {
            eps_ = DBL_EPSILON;
            epsx_ = Sq(eps_);
            ahypover_ = DBL_MANT_DIG * Log(FLT_RADIX) + 2;
            _a = a;
            _f = f;
            _fm = 1 - _f;
            _e2 = _f * (2 - _f);
            _es = (_f < 0 ? -1 : 1) * Sqrt(Abs(_e2));
        }

        /// <summary>
        /// Constructor with a single standard parallel.
        /// </summary>
        /// <param name="ellipsoid"><see cref="IEllipsoid"/> instance to be used in projection.</param>
        /// <param name="stdlat">standard parallel (degrees), the circle of tangency.</param>
        /// <param name="k0">scale on the standard parallel.</param>
        public LambertConformalConic(IEllipsoid ellipsoid, double stdlat, double k0)
            : this(ellipsoid.EquatorialRadius, ellipsoid.Flattening, stdlat, k0) { }

        /// <summary>
        /// Constructor with two standard parallels.
        /// </summary>
        /// <param name="ellipsoid"><see cref="IEllipsoid"/> instance to be used in projection.</param>
        /// <param name="stdlat1">first standard parallel (degrees).</param>
        /// <param name="stdlat2">second standard parallel (degrees).</param>
        /// <param name="k1">scale on the standard parallels.</param>
        /// <exception cref="GeographicException">
        /// When <i>a</i>, (1 − <i>f</i>) <i>a</i>, or <paramref name="k1"/> is not positive. 
        /// When <paramref name="stdlat1"/> or <paramref name="stdlat2"/> is not in [−90°, 90°], 
        /// or if either <paramref name="stdlat1"/> or <paramref name="stdlat2"/> is a pole and <paramref name="stdlat1"/> is not equal <paramref name="stdlat2"/>.
        /// </exception>
        public LambertConformalConic(IEllipsoid ellipsoid, double stdlat1, double stdlat2, double k1)
            : this(ellipsoid.EquatorialRadius, ellipsoid.Flattening, stdlat1, stdlat2, k1) { }

        /// <summary>
        /// Constructor with a single standard parallel.
        /// </summary>
        /// <param name="a">equatorial radius of ellipsoid (meters).</param>
        /// <param name="f">flattening of ellipsoid. Setting <i>f</i> = 0 gives a sphere. Negative <i>f</i> gives a prolate ellipsoid.</param>
        /// <param name="stdlat">standard parallel (degrees), the circle of tangency.</param>
        /// <param name="k0">scale on the standard parallel.</param>
        public LambertConformalConic(double a, double f, double stdlat, double k0)
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
                out _sign, out _n, out _nc, out _scbet0, out _tchi0, out _scchi0, out _psi0, out _lat0,
                out _t0nm1, out _scale, out _k0, out _nrho0, out _drhomax);
        }

        /// <summary>
        /// Constructor with two standard parallels.
        /// </summary>
        /// <param name="a">equatorial radius of ellipsoid (meters).</param>
        /// <param name="f">flattening of ellipsoid. Setting <i>f</i> = 0 gives a sphere. Negative <i>f</i> gives a prolate ellipsoid.</param>
        /// <param name="stdlat1">first standard parallel (degrees).</param>
        /// <param name="stdlat2">second standard parallel (degrees).</param>
        /// <param name="k1">scale on the standard parallels.</param>
        /// <exception cref="GeographicException">
        /// When <paramref name="a"/>, (1 − <paramref name="f"/>) <paramref name="a"/>, or <paramref name="k1"/> is not positive. 
        /// When <paramref name="stdlat1"/> or <paramref name="stdlat2"/> is not in [−90°, 90°], 
        /// or if either <paramref name="stdlat1"/> or <paramref name="stdlat2"/> is a pole and <paramref name="stdlat1"/> is not equal <paramref name="stdlat2"/>.
        /// </exception>
        public LambertConformalConic(double a, double f, double stdlat1, double stdlat2, double k1)
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
                out _sign, out _n, out _nc, out _scbet0, out _tchi0, out _scchi0, out _psi0, out _lat0,
                out _t0nm1, out _scale, out _k0, out _nrho0, out _drhomax);
        }

        /// <summary>
        /// Constructor with two standard parallels specified by sines and cosines.
        /// </summary>
        /// <param name="a">equatorial radius of ellipsoid (meters).</param>
        /// <param name="f">flattening of ellipsoid. Setting <i>f</i> = 0 gives a sphere. Negative <i>f</i> gives a prolate ellipsoid.</param>
        /// <param name="sinlat1">sine of first standard parallel.</param>
        /// <param name="coslat1">cosine of first standard parallel.</param>
        /// <param name="sinlat2">sine of second standard parallel.</param>
        /// <param name="coslat2">cosine of second standard parallel.</param>
        /// <param name="k1">scale on the standard parallels.</param>
        /// <remarks>
        /// This allows parallels close to the poles to be specified accurately.
        /// This routine computes the latitude of origin and the scale at this latitude.
        /// In the case where <i>lat1</i> and <i>lat2</i> are different, the errors in this routines are as follows:
        /// if <i>dlat</i> = abs(<i>lat2</i> − <i>lat1</i>) ≤ 160° and max(abs(<i>lat1</i>), 
        /// abs(<i>lat2</i>)) ≤ 90 − min(0.0002, 2.2 × 10^−6(180 − <i>dlat</i>), 6 * 10^−8 <i>dlat</i>^2) (in degrees), 
        /// then the error in the latitude of origin is less than 4.5 × 10^−14d and the relative error in the scale is less than 7 × 10^−15.
        /// </remarks>
        public LambertConformalConic(double a, double f,
                          double sinlat1, double coslat1,
                          double sinlat2, double coslat2,
                          double k1)
            : this(a, f)
        {
            if (!(IsFinite(_a) && _a > 0))
                throw new GeographicException("Equatorial radius is not positive");
            if (!(IsFinite(_f) && _f < 1))
                throw new GeographicException("Polar semi-axis is not positive");
            if (!(IsFinite(k1) && k1 > 0))
                throw new GeographicException("Scale is not positive");
            if (SignBit(coslat1))
                throw new GeographicException($"Standard latitude 1 not in [-{QD}d, {QD}d]");
            if (SignBit(coslat2))
                throw new GeographicException($"Standard latitude 2 not in [-{QD}d, {QD}d]");
            if (!(Abs(sinlat1) <= 1 && coslat1 <= 1) || (coslat1 == 0 && sinlat1 == 0))
                throw new GeographicException("Bad sine/cosine of standard latitude 1");
            if (!(Abs(sinlat2) <= 1 && coslat2 <= 1) || (coslat2 == 0 && sinlat2 == 0))
                throw new GeographicException("Bad sine/cosine of standard latitude 2");
            if (coslat1 == 0 || coslat2 == 0)
                if (!(coslat1 == coslat2 && sinlat1 == sinlat2))
                    throw new GeographicException
                      ("Standard latitudes must be equal is either is a pole");

            Init(sinlat1, coslat1, sinlat2, coslat2, k1,
                out _sign, out _n, out _nc, out _scbet0, out _tchi0, out _scchi0, out _psi0, out _lat0,
                out _t0nm1, out _scale, out _k0, out _nrho0, out _drhomax);
        }

        #endregion

        #region Private methods

        /// <summary>
        /// Divided differences
        /// Definition: Df(x,y) = (f(x)-f(y))/(x-y)
        /// <para>
        ///  See:
        ///   W. M. Kahan and R. J. Fateman,
        ///  Symbolic computation of divided differences,
        ///   SIGSAM Bull. 33(2), 7-28 (1999)
        ///   https://doi.org/10.1145/334714.334716
        ///   http://www.cs.berkeley.edu/~fateman/papers/divdiff.pdf
        /// </para>
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="hx"></param>
        /// <param name="hy"></param>
        /// <returns></returns>
        // General rules
        // h(x) = f(g(x)): Dh(x,y) = Df(g(x),g(y))*Dg(x,y)
        // h(x) = f(x)*g(x):
        //        Dh(x,y) = Df(x,y)*g(x) + Dg(x,y)*f(y)
        //                = Df(x,y)*g(y) + Dg(x,y)*f(x)
        //                = Df(x,y)*(g(x)+g(y))/2 + Dg(x,y)*(f(x)+f(y))/2
        //
        // hyp(x) = sqrt(1+x^2): Dhyp(x,y) = (x+y)/(hyp(x)+hyp(y))
        private static double Dhyp(double x, double y, double hx, double hy) => (x + y) / (hx + hy);

        // sn(x) = x/sqrt(1+x^2): Dsn(x,y) = (x+y)/((sn(x)+sn(y))*(1+x^2)*(1+y^2))
        private static double Dsn(double x, double y, double sx, double sy)
        {
            // sx = x/hyp(x)
            var t = x * y;
            return t > 0 ? (x + y) * Sq((sx * sy) / t) / (sx + sy) :
              (x - y != 0 ? (sx - sy) / (x - y) : 1);
        }

        // Dlog1p(x,y) = log1p((x-y)/(1+y))/(x-y)
        private static double Dlog1p(double x, double y)
        {
            var t = x - y; if (t < 0) { t = -t; y = x; }
            return t != 0 ? Log1p(t / (1 + y)) / t : 1 / (1 + x);
        }

        // Dexp(x,y) = exp((x+y)/2) * 2*sinh((x-y)/2)/(x-y)
        private static double Dexp(double x, double y)
        {
            var t = (x - y) / 2;
            return (t != 0 ? Sinh(t) / t : 1) * Exp((x + y) / 2);
        }

        // Dsinh(x,y) = 2*sinh((x-y)/2)/(x-y) * cosh((x+y)/2)
        //   cosh((x+y)/2) = (c+sinh(x)*sinh(y)/c)/2
        //   c=sqrt((1+cosh(x))*(1+cosh(y)))
        //   cosh((x+y)/2) = sqrt( (sinh(x)*sinh(y) + cosh(x)*cosh(y) + 1)/2 )
        private static double Dsinh(double x, double y, double sx, double sy, double cx, double cy)
        // sx = sinh(x), cx = cosh(x)
        {
            // double t = (x - y)/2, c = sqrt((1 + cx) * (1 + cy));
            // return (t ? sinh(t)/t : double(1)) * (c + sx * sy / c) /2;

            var t = (x - y) / 2;
            return (t != 0 ? Sinh(t) / t : 1) * Sqrt((sx * sy + cx * cy + 1) / 2);
        }

        // Dasinh(x,y) = asinh((x-y)*(x+y)/(x*sqrt(1+y^2)+y*sqrt(1+x^2)))/(x-y)
        //             = asinh((x*sqrt(1+y^2)-y*sqrt(1+x^2)))/(x-y)
        private static double Dasinh(double x, double y, double hx, double hy)
        {
            // hx = hyp(x)
            var t = x - y;
            return t != 0 ?
              Asinh(x * y > 0 ? t * (x + y) / (x * hy + y * hx) : x * hy - y * hx) / t :
              1 / hx;
        }

        // Deatanhe(x,y) = eatanhe((x-y)/(1-e^2*x*y))/(x-y)
        private double Deatanhe(double x, double y)
        {
            double t = x - y, d = 1 - _e2 * x * y;
            return t != 0 ? EAtanhE(t / d, _es) / t : _e2 / d;
        }

        private void Init(double sphi1, double cphi1, double sphi2, double cphi2, double k1,
            out double _sign, out double _n, out double _nc, out double _scbet0, out double _tchi0,
            out double _scchi0, out double _psi0, out double _lat0, out double _t0nm1, out double _scale,
            out double _k0, out double _nrho0, out double _drhomax)
        {
            {
                var r = Hypot(sphi1, cphi1);
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
              tphi1 = sphi1 / cphi1, tphi2 = sphi2 / cphi2, tphi0;
            //
            // Snyder: 15-8: n = (log(m1) - log(m2))/(log(t1)-log(t2))
            //
            // m = cos(bet) = 1/sec(bet) = 1/sqrt(1+tan(bet)^2)
            // bet = parametric lat, tan(bet) = (1-f)*tan(phi)
            //
            // t = tan(pi/4-chi/2) = 1/(sec(chi) + tan(chi)) = sec(chi) - tan(chi)
            // log(t) = -asinh(tan(chi)) = -psi
            // chi = conformal lat
            // tan(chi) = tan(phi)*cosh(xi) - sinh(xi)*sec(phi)
            // xi = eatanhe(sin(phi)), eatanhe(x) = e * atanh(e*x)
            //
            // n = (log(sec(bet2))-log(sec(bet1)))/(asinh(tan(chi2))-asinh(tan(chi1)))
            //
            // Let log(sec(bet)) = b(tphi), asinh(tan(chi)) = c(tphi)
            // Then n = Db(tphi2, tphi1)/Dc(tphi2, tphi1)
            // In limit tphi2 -> tphi1, n -> sphi1
            //
            double
              tbet1 = _fm * tphi1, scbet1 = Hypot(tbet1),
              tbet2 = _fm * tphi2, scbet2 = Hypot(tbet2);
            double
              scphi1 = 1 / cphi1,
              xi1 = EAtanhE(sphi1, _es), shxi1 = Sinh(xi1), chxi1 = Hypot(shxi1),
              tchi1 = chxi1 * tphi1 - shxi1 * scphi1, scchi1 = Hypot(tchi1),
              scphi2 = 1 / cphi2,
              xi2 = EAtanhE(sphi2, _es), shxi2 = Sinh(xi2), chxi2 = Hypot(shxi2),
              tchi2 = chxi2 * tphi2 - shxi2 * scphi2, scchi2 = Hypot(tchi2),
              psi1 = Asinh(tchi1);
            if (tphi2 - tphi1 != 0)
            {
                // Db(tphi2, tphi1)
                var num = Dlog1p(Sq(tbet2) / (1 + scbet2),
                                 Sq(tbet1) / (1 + scbet1))
                  * Dhyp(tbet2, tbet1, scbet2, scbet1) * _fm;
                // Dc(tphi2, tphi1)
                var den = Dasinh(tphi2, tphi1, scphi2, scphi1)
                  - Deatanhe(sphi2, sphi1) * Dsn(tphi2, tphi1, sphi2, sphi1);
                _n = num / den;

                if (_n < 1 / 4d)
                    _nc = Sqrt((1 - _n) * (1 + _n));
                else
                {
                    // Compute nc = cos(phi0) = sqrt((1 - n) * (1 + n)), evaluating 1 - n
                    // carefully.  First write
                    //
                    // Dc(tphi2, tphi1) * (tphi2 - tphi1)
                    //   = log(tchi2 + scchi2) - log(tchi1 + scchi1)
                    //
                    // then den * (1 - n) =
                    // (log((tchi2 + scchi2)/(2*scbet2)) -
                    //  log((tchi1 + scchi1)/(2*scbet1))) / (tphi2 - tphi1)
                    // = Dlog1p(a2, a1) * (tchi2+scchi2 + tchi1+scchi1)/(4*scbet1*scbet2)
                    //   * fm * Q
                    //
                    // where
                    // a1 = ( (tchi1 - scbet1) + (scchi1 - scbet1) ) / (2 * scbet1)
                    // Q = ((scbet2 + scbet1)/fm)/((scchi2 + scchi1)/D(tchi2, tchi1))
                    //     - (tbet2 + tbet1)/(scbet2 + scbet1)
                    double t;
                    {
                        double
                          // s1 = (scbet1 - scchi1) * (scbet1 + scchi1)
                          s1 = (tphi1 * (2 * shxi1 * chxi1 * scphi1 - _e2 * tphi1) -
                                Sq(shxi1) * (1 + 2 * Sq(tphi1))),
                          s2 = (tphi2 * (2 * shxi2 * chxi2 * scphi2 - _e2 * tphi2) -
                                Sq(shxi2) * (1 + 2 * Sq(tphi2))),
                          // t1 = scbet1 - tchi1
                          t1 = tchi1 < 0 ? scbet1 - tchi1 : (s1 + 1) / (scbet1 + tchi1),
                          t2 = tchi2 < 0 ? scbet2 - tchi2 : (s2 + 1) / (scbet2 + tchi2),
                          a2 = -(s2 / (scbet2 + scchi2) + t2) / (2 * scbet2),
                          a1 = -(s1 / (scbet1 + scchi1) + t1) / (2 * scbet1);
                        t = Dlog1p(a2, a1) / den;
                    }
                    // multiply by (tchi2 + scchi2 + tchi1 + scchi1)/(4*scbet1*scbet2) * fm
                    t *= (((tchi2 >= 0 ? scchi2 + tchi2 : 1 / (scchi2 - tchi2)) +
                             (tchi1 >= 0 ? scchi1 + tchi1 : 1 / (scchi1 - tchi1))) /
                           (4 * scbet1 * scbet2)) * _fm;

                    // Rewrite
                    // Q = (1 - (tbet2 + tbet1)/(scbet2 + scbet1)) -
                    //     (1 - ((scbet2 + scbet1)/fm)/((scchi2 + scchi1)/D(tchi2, tchi1)))
                    //   = tbm - tam
                    // where
                    var tbm = (((tbet1 > 0 ? 1 / (scbet1 + tbet1) : scbet1 - tbet1) +
                                  (tbet2 > 0 ? 1 / (scbet2 + tbet2) : scbet2 - tbet2)) /
                                 (scbet1 + scbet2));

                    // tam = (1 - ((scbet2+scbet1)/fm)/((scchi2+scchi1)/D(tchi2, tchi1)))
                    //
                    // Let
                    //   (scbet2 + scbet1)/fm = scphi2 + scphi1 + dbet
                    //   (scchi2 + scchi1)/D(tchi2, tchi1) = scphi2 + scphi1 + dchi
                    // then
                    //   tam = D(tchi2, tchi1) * (dchi - dbet) / (scchi1 + scchi2)
                    double
                      // D(tchi2, tchi1)
                      dtchi = den / Dasinh(tchi2, tchi1, scchi2, scchi1),
                      // (scbet2 + scbet1)/fm - (scphi2 + scphi1)
                      dbet = (_e2 / _fm) * (1 / (scbet2 + _fm * scphi2) +
                                           1 / (scbet1 + _fm * scphi1));

                    // dchi = (scchi2 + scchi1)/D(tchi2, tchi1) - (scphi2 + scphi1)
                    // Let
                    //    tzet = chxiZ * tphi - shxiZ * scphi
                    //    tchi = tzet + nu
                    //    scchi = sczet + mu
                    // where
                    //    xiZ = eatanhe(1), shxiZ = sinh(xiZ), chxiZ = cosh(xiZ)
                    //    nu =   scphi * (shxiZ - shxi) - tphi * (chxiZ - chxi)
                    //    mu = - scphi * (chxiZ - chxi) + tphi * (shxiZ - shxi)
                    // then
                    // dchi = ((mu2 + mu1) - D(nu2, nu1) * (scphi2 + scphi1)) /
                    //         D(tchi2, tchi1)
                    double
                      xiZ = EAtanhE(1, _es),
                      shxiZ = Sinh(xiZ), chxiZ = Hypot(shxiZ),
                      // These are differences not divided differences
                      // dxiZ1 = xiZ - xi1; dshxiZ1 = shxiZ - shxi; dchxiZ1 = chxiZ - chxi
                      dxiZ1 = Deatanhe(1, sphi1) / (scphi1 * (tphi1 + scphi1)),
                      dxiZ2 = Deatanhe(1, sphi2) / (scphi2 * (tphi2 + scphi2)),
                      dshxiZ1 = Dsinh(xiZ, xi1, shxiZ, shxi1, chxiZ, chxi1) * dxiZ1,
                      dshxiZ2 = Dsinh(xiZ, xi2, shxiZ, shxi2, chxiZ, chxi2) * dxiZ2,
                      dchxiZ1 = Dhyp(shxiZ, shxi1, chxiZ, chxi1) * dshxiZ1,
                      dchxiZ2 = Dhyp(shxiZ, shxi2, chxiZ, chxi2) * dshxiZ2,
                      // mu1 + mu2
                      amu12 = (-scphi1 * dchxiZ1 + tphi1 * dshxiZ1
                               - scphi2 * dchxiZ2 + tphi2 * dshxiZ2),
                      // D(xi2, xi1)
                      dxi = Deatanhe(sphi1, sphi2) * Dsn(tphi2, tphi1, sphi2, sphi1),
                      // D(nu2, nu1)
                      dnu12 =
                      ((_f * 4 * scphi2 * dshxiZ2 > _f * scphi1 * dshxiZ1 ?
                         // Use divided differences
                         (dshxiZ1 + dshxiZ2) / 2 * Dhyp(tphi1, tphi2, scphi1, scphi2)
                         - ((scphi1 + scphi2) / 2
                             * Dsinh(xi1, xi2, shxi1, shxi2, chxi1, chxi2) * dxi) :
                         // Use ratio of differences
                         (scphi2 * dshxiZ2 - scphi1 * dshxiZ1) / (tphi2 - tphi1))
                        + ((tphi1 + tphi2) / 2 * Dhyp(shxi1, shxi2, chxi1, chxi2)
                            * Dsinh(xi1, xi2, shxi1, shxi2, chxi1, chxi2) * dxi)
                        - (dchxiZ1 + dchxiZ2) / 2),
                      // dtchi * dchi
                      dchia = (amu12 - dnu12 * (scphi2 + scphi1)),
                      tam = (dchia - dtchi * dbet) / (scchi1 + scchi2);
                    t *= tbm - tam;
                    _nc = Sqrt(Max(0, t) * (1 + _n));
                }
                {
                    var r = Hypot(_n, _nc);
                    _n /= r;
                    _nc /= r;
                }
                tphi0 = _n / _nc;
            }
            else
            {
                tphi0 = tphi1;
                _nc = 1 / Hypot(tphi0);
                _n = tphi0 * _nc;
                if (polar)
                    _nc = 0;
            }

            _scbet0 = Hypot(_fm * tphi0);
            var shxi0 = Sinh(EAtanhE(_n, _es));
            _tchi0 = tphi0 * Hypot(shxi0) - shxi0 * Hypot(tphi0); _scchi0 = Hypot(_tchi0);
            _psi0 = Asinh(_tchi0);

            _lat0 = Atan(_sign * tphi0) / Degree;
            _t0nm1 = Expm1(-_n * _psi0); // Snyder's t0^n - 1
                                         // a * k1 * m1/t1^n = a * k1 * m2/t2^n = a * k1 * n * (Snyder's F)
                                         // = a * k1 / (scbet1 * exp(-n * psi1))
            _scale = _a * k1 / scbet1 *
              // exp(n * psi1) = exp(- (1 - n) * psi1) * exp(psi1)
              // with (1-n) = nc^2/(1+n) and exp(-psi1) = scchi1 + tchi1
              Exp(-(Sq(_nc) / (1 + _n)) * psi1)
              * (tchi1 >= 0 ? scchi1 + tchi1 : 1 / (scchi1 - tchi1));
            // Scale at phi0 = k0 = k1 * (scbet0*exp(-n*psi0))/(scbet1*exp(-n*psi1))
            //                    = k1 * scbet0/scbet1 * exp(n * (psi1 - psi0))
            // psi1 - psi0 = Dasinh(tchi1, tchi0) * (tchi1 - tchi0)
            _k0 = k1 * (_scbet0 / scbet1) *
              Exp(-(Sq(_nc) / (1 + _n)) *
                   Dasinh(tchi1, _tchi0, scchi1, _scchi0) * (tchi1 - _tchi0))
              * (tchi1 >= 0 ? scchi1 + tchi1 : 1 / (scchi1 - tchi1)) /
              (_scchi0 + _tchi0);
            _nrho0 = polar ? 0 : _a * _k0 / _scbet0;
            {
                // Figure _drhomax using code at beginning of Forward with lat = -90
                double
                  sphi = -1, cphi = epsx_,
                  tphi = sphi / cphi,
                  scphi = 1 / cphi, shxi = Sinh(EAtanhE(sphi, _es)),
                  tchi = Hypot(shxi) * tphi - shxi * scphi, scchi = Hypot(tchi),
                  psi = Asinh(tchi),
                  dpsi = Dasinh(tchi, _tchi0, scchi, _scchi0) * (tchi - _tchi0);
                _drhomax = -_scale * (2 * _nc < 1 && dpsi != 0 ?
                                       (Exp(Sq(_nc) / (1 + _n) * psi) *
                                        (tchi > 0 ? 1 / (scchi + tchi) : (scchi - tchi))
                                        - (_t0nm1 + 1)) / (-_n) :
                                       Dexp(-_n * psi, -_n * _psi0) * dpsi);
            }
        }

        #endregion

        #region Public properties

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
        /// Gets a value representing the equatorial radius (<i>a</i>) of the ellipsoid.
        /// </summary>
        public double EquatorialRadius => _a;

        /// <summary>
        /// Gets a value representing the flatterning (<i>f</i>) of the ellipsoid.
        /// </summary>
        public double Flattening => _f;

        /// <summary>
        /// Gets or sets an value representing that whether current <see cref="LambertConformalConic"/> is frozen.
        /// </summary>
        public bool IsFrozen { get; private set; }

        #endregion

        #region Public methods

        /// <summary>
        /// Freeze current <see cref="LambertConformalConic"/> instance to prevent its scale being modified.
        /// </summary>
        public LambertConformalConic Freeze() { IsFrozen = true; return this; }

        /// <summary>
        /// Set the scale for the projection.
        /// </summary>
        /// <param name="lat">(degrees).</param>
        /// <param name="k">scale at latitude <paramref name="lat"/> (default 1).</param>
        public void SetScale(double lat, double k = 1d)
        {
            if (IsFrozen)
                throw new GeographicException("Projection is frozen");
            if (!(IsFinite(k) && k > 0))
                throw new GeographicException("Scale is not positive");
            if (!(Abs(lat) <= QD))
                throw new GeographicException($"Latitude for SetScale not in [-{QD}d, {QD}d]");
            if (Abs(lat) == QD && !(_nc == 0 && lat * _n > 0))
                throw new GeographicException("Incompatible polar latitude in SetScale");

            Forward(0, lat, 0, out _, out var kold);
            k /= kold;
            _scale *= k;
            _k0 *= k;
        }

        /// <summary>
        /// Forward projection, from geographic to Lambert conformal conic.
        /// </summary>
        /// <param name="lon0">central meridian longitude (degrees).</param>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <param name="gamma">meridian convergence at point (degrees).</param>
        /// <param name="k">scale of projection at point.</param>
        /// <returns>
        /// <i>x</i>, easting of point and <i>y</i>, northing of point, in meters.
        /// </returns>
        /// <remarks>
        /// The latitude origin is given by <see cref="OriginLatitude"/>. 
        /// No false easting or northing is added and <paramref name="lat"/> should be in the range [−90°, 90°].
        /// The error in the projection is less than about 10 nm (10 nanometers), true distance, 
        /// and the errors in the meridian convergence and scale are consistent with this. 
        /// The values of <i>x</i> and <i>y</i> returned for points which project to infinity
        /// (i.e., one or both of the poles) will be large but finite.
        /// </remarks>
        public (double x, double y) Forward(double lon0, double lat, double lon, out double gamma, out double k)
        {
            lon = AngDiff(lon0, lon);
            // From Snyder, we have
            //
            // theta = n * lambda
            // x = rho * sin(theta)
            //   = (nrho0 + n * drho) * sin(theta)/n
            // y = rho0 - rho * cos(theta)
            //   = nrho0 * (1-cos(theta))/n - drho * cos(theta)
            //
            // where nrho0 = n * rho0, drho = rho - rho0
            // and drho is evaluated with divided differences
            SinCosd(LatFix(lat) * _sign, out var sphi, out var cphi);
            cphi = Max(epsx_, cphi);
            double
              lam = lon * Degree,
              tphi = sphi / cphi, scbet = Hypot(_fm * tphi),
              scphi = 1 / cphi, shxi = Sinh(EAtanhE(sphi, _es)),
              tchi = Hypot(shxi) * tphi - shxi * scphi, scchi = Hypot(tchi),
              psi = Asinh(tchi),
              theta = _n * lam, stheta = Sin(theta), ctheta = Cos(theta),
              dpsi = Dasinh(tchi, _tchi0, scchi, _scchi0) * (tchi - _tchi0),
              drho = -_scale * (2 * _nc < 1 && dpsi != 0 ?
                                 (Exp(Sq(_nc) / (1 + _n) * psi) *
                                  (tchi > 0 ? 1 / (scchi + tchi) : (scchi - tchi))
                                  - (_t0nm1 + 1)) / (-_n) :
                                 Dexp(-_n * psi, -_n * _psi0) * dpsi);
            var x = (_nrho0 + _n * drho) * (_n != 0 ? stheta / _n : lam);
            var y = _nrho0 *
              (_n != 0 ?
               (ctheta < 0 ? 1 - ctheta : Sq(stheta) / (1 + ctheta)) / _n : 0)
              - drho * ctheta;
            k = _k0 * (scbet / _scbet0) /
              (Exp(-(Sq(_nc) / (1 + _n)) * dpsi)
               * (tchi >= 0 ? scchi + tchi : 1 / (scchi - tchi)) / (_scchi0 + _tchi0));
            y *= _sign;
            gamma = _sign * theta / Degree;

            return (x, y);
        }

        /// <summary>
        /// Reverse projection, from Lambert conformal conic to geographic.
        /// </summary>
        /// <param name="lon0">central meridian longitude (degrees).</param>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <param name="gamma">meridian convergence at point (degrees).</param>
        /// <param name="k">scale of projection at point.</param>
        /// <returns>
        /// <i>lat</i>, latitude of point and <i>lon</i>, longitude of point, in degress.
        /// </returns>
        /// <remarks>
        /// The latitude origin is given by <see cref="OriginLatitude"/>.
        /// No false easting or northing is added. 
        /// The value of <i>lon</i> returned is in the range [−180°, 180°].
        /// The error in the projection is less than about 10 nm (10 nanometers), true distance, 
        /// and the errors in the meridian convergence and scale are consistent with this.
        /// </remarks>
        public (double lat, double lon) Reverse(double lon0, double x, double y, out double gamma, out double k)
        {
            // From Snyder, we have
            //
            //        x = rho * sin(theta)
            // rho0 - y = rho * cos(theta)
            //
            // rho = hypot(x, rho0 - y)
            // drho = (n*x^2 - 2*y*nrho0 + n*y^2)/(hypot(n*x, nrho0-n*y) + nrho0)
            // theta = atan2(n*x, nrho0-n*y)
            //
            // From drho, obtain t^n-1
            // psi = -log(t), so
            // dpsi = - Dlog1p(t^n-1, t0^n-1) * drho / scale
            y *= _sign;
            double
              // Guard against 0 * inf in computation of ny
              nx = _n * x, ny = _n != 0 ? _n * y : 0, y1 = _nrho0 - ny,
              den = Hypot(nx, y1) + _nrho0, // 0 implies origin with polar aspect
                                            // isfinite test is to avoid inf/inf
              drho = ((den != 0 && IsFinite(den))
                      ? (x * nx + y * (ny - 2 * _nrho0)) / den
                      : den);
            drho = Min(drho, _drhomax);
            if (_n == 0)
                drho = Max(drho, -_drhomax);
            double
              tnm1 = _t0nm1 + _n * drho / _scale,
              dpsi = (den == 0 ? 0 :
                      (tnm1 + 1 != 0 ? -Dlog1p(tnm1, _t0nm1) * drho / _scale :
                       ahypover_));
            double tchi;
            if (2 * _n <= 1)
            {
                // tchi = sinh(psi)
                double
                  psi = _psi0 + dpsi, tchia = Sinh(psi), scchip = Hypot(tchia),
                  dtchi = Dsinh(psi, _psi0, tchia, _tchi0, scchip, _scchi0) * dpsi;
                tchi = _tchi0 + dtchi;    // Update tchi using divided difference
            }
            else
            {
                // tchi = sinh(-1/n * log(tn))
                //      = sinh((1-1/n) * log(tn) - log(tn))
                //      = + sinh((1-1/n) * log(tn)) * cosh(log(tn))
                //        - cosh((1-1/n) * log(tn)) * sinh(log(tn))
                // (1-1/n) = - nc^2/(n*(1+n))
                // cosh(log(tn)) = (tn + 1/tn)/2; sinh(log(tn)) = (tn - 1/tn)/2
                double
                  tn = tnm1 + 1 == 0 ? epsx_ : tnm1 + 1,
                  sh = Sinh(-Sq(_nc) / (_n * (1 + _n)) *
                             (2 * tn > 1 ? Log1p(tnm1) : Log(tn)));
                tchi = sh * (tn + 1 / tn) / 2 - Hypot(sh) * (tnm1 * (tn + 1) / tn) / 2;
            }

            // log(t) = -asinh(tan(chi)) = -psi
            gamma = Atan2(nx, y1);
            double
              tphi = Tauf(tchi, _es),
              scbet = Hypot(_fm * tphi), scchi = Hypot(tchi),
              lam = _n != 0 ? gamma / _n : x / y1;
            var lat = Atand(_sign * tphi);
            var lon = lam / Degree;
            lon = AngNormalize(lon + AngNormalize(lon0));
            k = _k0 * (scbet / _scbet0) /
              (Exp(_nc != 0 ? -(Sq(_nc) / (1 + _n)) * dpsi : 0)
               * (tchi >= 0 ? scchi + tchi : 1 / (scchi - tchi)) / (_scchi0 + _tchi0));
            gamma /= _sign * Degree;

            return (lat, lon);
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
        public (double x, double y) Forward(double lon0, double lat, double lon) => Forward(lon0, lat, lon, out _, out _);

        /// <summary>
        /// Reverse without returning convergence and scale.
        /// </summary>
        /// <param name="lon0">central meridian longitude (degrees).</param>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <returns>
        /// <i>lat</i>, latitude of point and <i>lon</i>, longitude of point, in degress.
        /// </returns>
        public (double lat, double lon) Reverse(double lon0, double x, double y) => Reverse(lon0, x, y, out _, out _);

        #endregion
    }
}