using System;
using System.Collections.Generic;
using System.Text;
using System.Numerics;

using static System.Math;
using static GeographicLib.MathEx;
using static GeographicLib.Macros;

namespace GeographicLib.Projections
{
    /// <summary>
    /// Transverse Mercator projection.
    /// This uses Krüger's method which evaluates the projection and its inverse in terms of a series. See 
    /// <para>
    ///   - L. Krüger, <a href = "https://doi.org/10.2312/GFZ.b103-krueger28" > 
    ///   Konforme Abbildung des Erdellipsoids in der Ebene</a>
    ///   (Conformal mapping of the ellipsoidal earth to the plane), 
    ///   Royal Prussian Geodetic Institute, New Series 52, 172 pp. (1912).
    /// </para>
    /// <para>
    /// - C. F. F. Karney,
    /// <a href="https://doi.org/10.1007/s00190-011-0445-3">
    /// Transverse Mercator with an accuracy of a few nanometers,</a>
    /// J. Geodesy 85(8), 475--485 (Aug. 2011);
    /// preprint
    /// <a href="https://arxiv.org/abs/1002.1417">arXiv:1002.1417</a>.
    /// </para>
    /// </summary>
    public class TransverseMercator : IEllipsoid
    {
        private const int _maxpow = GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER;
        // private const int _numit = 5;

        private static readonly Memory<double> b1coeff, alpcoeff, betcoeff;

        private readonly double _a, _f, _k0, _e2, _es, _e2m, _c, _n;
        private readonly double _a1, _b1;

        internal readonly Memory<double> _alp = new double[_maxpow + 1], _bet = new double[_maxpow + 1];

        static TransverseMercator()
        {
#pragma warning disable CS0162
            switch (GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER / 2)
            {
                case 2: b1coeff = new double[] { 1, 16, 64, 64 }; break;
                case 3: b1coeff = new double[] { 1, 4, 64, 256, 256 }; break;
                case 4: b1coeff = new double[] { 25, 64, 256, 4096, 16384, 16384 }; break;
                default: throw new GeographicException("Bad value for GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER");
            }

            switch (GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER)
            {
                case 4:
                    alpcoeff = new double[]
                    {
                          // alp[1]/n^1, polynomial in n of order 3
                          164, 225, -480, 360, 720,
                          // alp[2]/n^2, polynomial in n of order 2
                          557, -864, 390, 1440,
                          // alp[3]/n^3, polynomial in n of order 1
                          -1236, 427, 1680,
                          // alp[4]/n^4, polynomial in n of order 0
                          49561, 161280,
                    };
                    betcoeff = new double[]
                    {
                          // bet[1]/n^1, polynomial in n of order 3
                          -4, 555, -960, 720, 1440,
                          // bet[2]/n^2, polynomial in n of order 2
                          -437, 96, 30, 1440,
                          // bet[3]/n^3, polynomial in n of order 1
                          -148, 119, 3360,
                          // bet[4]/n^4, polynomial in n of order 0
                          4397, 161280,
                    };
                    break;
                case 5:
                    alpcoeff = new double[]
                    {
                          // alp[1]/n^1, polynomial in n of order 4
                          -635, 328, 450, -960, 720, 1440,
                          // alp[2]/n^2, polynomial in n of order 3
                          4496, 3899, -6048, 2730, 10080,
                          // alp[3]/n^3, polynomial in n of order 2
                          15061, -19776, 6832, 26880,
                          // alp[4]/n^4, polynomial in n of order 1
                          -171840, 49561, 161280,
                          // alp[5]/n^5, polynomial in n of order 0
                          34729, 80640,
                    };
                    betcoeff = new double[]
                    {
                          // bet[1]/n^1, polynomial in n of order 4
                          -3645, -64, 8880, -15360, 11520, 23040,
                          // bet[2]/n^2, polynomial in n of order 3
                          4416, -3059, 672, 210, 10080,
                          // bet[3]/n^3, polynomial in n of order 2
                          -627, -592, 476, 13440,
                          // bet[4]/n^4, polynomial in n of order 1
                          -3520, 4397, 161280,
                          // bet[5]/n^5, polynomial in n of order 0
                          4583, 161280,
                    };
                    break;
                case 6:
                    alpcoeff = new double[]
                    {
                          // alp[1]/n^1, polynomial in n of order 5
                          31564, -66675, 34440, 47250, -100800, 75600, 151200,
                          // alp[2]/n^2, polynomial in n of order 4
                          -1983433, 863232, 748608, -1161216, 524160, 1935360,
                          // alp[3]/n^3, polynomial in n of order 3
                          670412, 406647, -533952, 184464, 725760,
                          // alp[4]/n^4, polynomial in n of order 2
                          6601661, -7732800, 2230245, 7257600,
                          // alp[5]/n^5, polynomial in n of order 1
                          -13675556, 3438171, 7983360,
                          // alp[6]/n^6, polynomial in n of order 0
                          212378941, 319334400,
                    };
                    betcoeff = new double[]
                    {
                          // bet[1]/n^1, polynomial in n of order 5
                          384796, -382725, -6720, 932400, -1612800, 1209600, 2419200,
                          // bet[2]/n^2, polynomial in n of order 4
                          -1118711, 1695744, -1174656, 258048, 80640, 3870720,
                          // bet[3]/n^3, polynomial in n of order 3
                          22276, -16929, -15984, 12852, 362880,
                          // bet[4]/n^4, polynomial in n of order 2
                          -830251, -158400, 197865, 7257600,
                          // bet[5]/n^5, polynomial in n of order 1
                          -435388, 453717, 15966720,
                          // bet[6]/n^6, polynomial in n of order 0
                          20648693, 638668800,
                    };
                    break;
                case 7:
                    alpcoeff = new double[]
                    {
                          // alp[1]/n^1, polynomial in n of order 6
                          1804025, 2020096, -4267200, 2204160, 3024000, -6451200, 4838400, 9676800,
                          // alp[2]/n^2, polynomial in n of order 5
                          4626384, -9917165, 4316160, 3743040, -5806080, 2620800, 9676800,
                          // alp[3]/n^3, polynomial in n of order 4
                          -67102379, 26816480, 16265880, -21358080, 7378560, 29030400,
                          // alp[4]/n^4, polynomial in n of order 3
                          155912000, 72618271, -85060800, 24532695, 79833600,
                          // alp[5]/n^5, polynomial in n of order 2
                          102508609, -109404448, 27505368, 63866880,
                          // alp[6]/n^6, polynomial in n of order 1
                          -12282192400L, 2760926233L, 4151347200L,
                          // alp[7]/n^7, polynomial in n of order 0
                          1522256789, 1383782400,
                    };
                    betcoeff = new double[]
                    {
                          // bet[1]/n^1, polynomial in n of order 6
                          -5406467, 6156736, -6123600, -107520, 14918400, -25804800, 19353600,
                          38707200,
                          // bet[2]/n^2, polynomial in n of order 5
                          829456, -5593555, 8478720, -5873280, 1290240, 403200, 19353600,
                          // bet[3]/n^3, polynomial in n of order 4
                          9261899, 3564160, -2708640, -2557440, 2056320, 58060800,
                          // bet[4]/n^4, polynomial in n of order 3
                          14928352, -9132761, -1742400, 2176515, 79833600,
                          // bet[5]/n^5, polynomial in n of order 2
                          -8005831, -1741552, 1814868, 63866880,
                          // bet[6]/n^6, polynomial in n of order 1
                          -261810608, 268433009, 8302694400L,
                          // bet[7]/n^7, polynomial in n of order 0
                          219941297, 5535129600L,
                    };
                    break;
                case 8:
                    alpcoeff = new double[]
                    {
                          // alp[1]/n^1, polynomial in n of order 7
                          -75900428, 37884525, 42422016, -89611200, 46287360, 63504000, -135475200,
                          101606400, 203212800,
                          // alp[2]/n^2, polynomial in n of order 6
                          148003883, 83274912, -178508970, 77690880, 67374720, -104509440,
                          47174400, 174182400,
                          // alp[3]/n^3, polynomial in n of order 5
                          318729724, -738126169, 294981280, 178924680, -234938880, 81164160,
                          319334400,
                          // alp[4]/n^4, polynomial in n of order 4
                          -40176129013L, 14967552000L, 6971354016L, -8165836800L, 2355138720L,
                          7664025600L,
                          // alp[5]/n^5, polynomial in n of order 3
                          10421654396L, 3997835751L, -4266773472L, 1072709352, 2490808320L,
                          // alp[6]/n^6, polynomial in n of order 2
                          175214326799L, -171950693600L, 38652967262L, 58118860800L,
                          // alp[7]/n^7, polynomial in n of order 1
                          -67039739596L, 13700311101L, 12454041600L,
                          // alp[8]/n^8, polynomial in n of order 0
                          1424729850961L, 743921418240L,
                    };
                    betcoeff = new double[]
                    {
                          // bet[1]/n^1, polynomial in n of order 7
                          31777436, -37845269, 43097152, -42865200, -752640, 104428800, -180633600,
                          135475200, 270950400,
                          // bet[2]/n^2, polynomial in n of order 6
                          24749483, 14930208, -100683990, 152616960, -105719040, 23224320, 7257600,
                          348364800,
                          // bet[3]/n^3, polynomial in n of order 5
                          -232468668, 101880889, 39205760, -29795040, -28131840, 22619520,
                          638668800,
                          // bet[4]/n^4, polynomial in n of order 4
                          324154477, 1433121792, -876745056, -167270400, 208945440, 7664025600L,
                          // bet[5]/n^5, polynomial in n of order 3
                          457888660, -312227409, -67920528, 70779852, 2490808320L,
                          // bet[6]/n^6, polynomial in n of order 2
                          -19841813847L, -3665348512L, 3758062126L, 116237721600L,
                          // bet[7]/n^7, polynomial in n of order 1
                          -1989295244, 1979471673, 49816166400L,
                          // bet[8]/n^8, polynomial in n of order 0
                          191773887257L, 3719607091200L,
                    };
                default:
                    throw new GeographicException("Bad value for GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER");
            }
#pragma warning restore

            if (b1coeff.Length != _maxpow / 2 + 2)
            {
                throw new GeographicException("Coefficient array size mismatch for b1");
            }

            if (alpcoeff.Length != (_maxpow * (_maxpow + 3)) / 2)
            {
                throw new GeographicException("Coefficient array size mismatch for alp");
            }

            if (betcoeff.Length != (_maxpow * (_maxpow + 3)) / 2)
            {
                throw new GeographicException("Coefficient array size mismatch for bet");
            }

            UTM = new TransverseMercator(Constants.WGS84_a, Constants.WGS84_f, Constants.UTM_k0);
        }

        /// <summary>
        /// Initialize a new <see cref="TransverseMercator"/> instance with specified equatorial radius and flattening.
        /// </summary>
        /// <param name="a">equatorial radius (meters).</param>
        /// <param name="f">flattening of ellipsoid. Setting <i>f</i> = 0 gives a sphere. Negative <i>f</i> gives a prolate ellipsoid.</param>
        /// <param name="k0">central scale factor.</param>
        public TransverseMercator(double a, double f, double k0)
        {
            _a = a;
            _f = f;
            _k0 = k0;
            _e2 = _f * (2 - _f);
            _es = (_f < 0 ? -1 : 1) * Sqrt(Abs(_e2));
            _e2m = 1 - _e2;

            // _c = sqrt( pow(1 + _e, 1 + _e) * pow(1 - _e, 1 - _e) ) )
            // See, for example, Lee (1976), p 100.
            _c = Sqrt(_e2m) * Exp(EAtanhE(1d, _es));

            _n = _f / (2 - _f);

            if (!(IsFinite(_a) && _a > 0))
                throw new GeographicException("Equatorial radius is not positive");
            if (!(IsFinite(_f) && _f < 1))
                throw new GeographicException("Polar semi-axis is not positive");
            if (!(IsFinite(_k0) && _k0 > 0))
                throw new GeographicException("Scale is not positive");

            var m = _maxpow / 2;
            _b1 = PolyVal(m, b1coeff, Sq(_n)) / (b1coeff.Span[m + 1] * (1 + _n));

            // _a1 is the equivalent radius for computing the circumference of
            // ellipse.
            _a1 = _b1 * _a;
            var o = 0;
            var d = _n;
            for (int l = 1; l <= _maxpow; ++l)
            {
                m = _maxpow - l;
                _alp.Span[l] = d * PolyVal(m, alpcoeff.Slice(o), _n) / alpcoeff.Span[o + m + 1];
                _bet.Span[l] = d * PolyVal(m, betcoeff.Slice(o), _n) / betcoeff.Span[o + m + 1];
                o += m + 2;
                d *= _n;
            }

            // Post condition: o == sizeof(alpcoeff) / sizeof(real) &&
            // o == sizeof(betcoeff) / sizeof(real)
        }

        /// <summary>
        /// Initialize a new <see cref="TransverseMercator"/> instance with specified ellipsoid.
        /// </summary>
        /// <param name="ellipsoid"><see cref="IEllipsoid"/> instance to be used in projection.</param>
        /// <param name="k0">central scale factor.</param>
        public TransverseMercator(IEllipsoid ellipsoid, double k0)
            : this(ellipsoid.EquatorialRadius, ellipsoid.Flattening, k0) { }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static TransverseMercator UTM { get; }

        /// <summary>
        /// Gets a value representing equatorial radius of the ellipsoid. 
        /// </summary>
        public double EquatorialRadius => _a;

        /// <summary>
        /// Gets a value representing flatterning of of the ellipsoid. 
        /// </summary>
        public double Flattening => _f;

        /// <summary>
        /// Gets a value representing central scale for the projection.  This is the azimuthal scale on the latitude of origin.
        /// </summary>
        public double CentralScale => _k0;

        /// <summary>
        /// Forward projection, from geographic to projected coordinate system.
        /// </summary>
        /// <param name="lon0">central meridian of the projection (degrees).</param>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <param name="gamma">meridian convergence at point (degrees).</param>
        /// <param name="k">scale of projection at point.</param>
        /// <returns>
        /// <i>x</i>, easting of point (meters) and <i>y</i>, northing of point (meters).
        /// </returns>
        public (double x, double y) Forward(double lon0, double lat, double lon, out double gamma, out double k)
        {
            lat = LatFix(lat);
            lon = AngDiff(lon0, lon);
            // Explicitly enforce the parity
            int
              latsign = (lat < 0) ? -1 : 1,
              lonsign = (lon < 0) ? -1 : 1;
            lon *= lonsign;
            lat *= latsign;
            bool backside = lon > 90;
            if (backside)
            {
                if (lat == 0)
                    latsign = -1;
                lon = 180 - lon;
            }

            SinCosd(lat, out var sphi, out var cphi);
            SinCosd(lon, out var slam, out var clam);
            // phi = latitude
            // phi' = conformal latitude
            // psi = isometric latitude
            // tau = tan(phi)
            // tau' = tan(phi')
            // [xi', eta'] = Gauss-Schreiber TM coordinates
            // [xi, eta] = Gauss-Krueger TM coordinates
            //
            // We use
            //   tan(phi') = sinh(psi)
            //   sin(phi') = tanh(psi)
            //   cos(phi') = sech(psi)
            //   denom^2    = 1-cos(phi')^2*sin(lam)^2 = 1-sech(psi)^2*sin(lam)^2
            //   sin(xip)   = sin(phi')/denom          = tanh(psi)/denom
            //   cos(xip)   = cos(phi')*cos(lam)/denom = sech(psi)*cos(lam)/denom
            //   cosh(etap) = 1/denom                  = 1/denom
            //   sinh(etap) = cos(phi')*sin(lam)/denom = sech(psi)*sin(lam)/denom
            double etap, xip;
            if (lat != 90)
            {
                double
                  tau = sphi / cphi,
                  taup = Taupf(tau, _es);
                xip = Atan2(taup, clam);
                // Used to be
                //   etap = Math::atanh(sin(lam) / cosh(psi));
                etap = Asinh(slam / Hypot(taup, clam));
                // convergence and scale for Gauss-Schreiber TM (xip, etap) -- gamma0 =
                // atan(tan(xip) * tanh(etap)) = atan(tan(lam) * sin(phi'));
                // sin(phi') = tau'/sqrt(1 + tau'^2)
                // Krueger p 22 (44)
                gamma = Atan2d(slam * taup, clam * Hypot(1d, taup));
                // k0 = sqrt(1 - _e2 * sin(phi)^2) * (cos(phi') / cos(phi)) * cosh(etap)
                // Note 1/cos(phi) = cosh(psip);
                // and cos(phi') * cosh(etap) = 1/hypot(sinh(psi), cos(lam))
                //
                // This form has cancelling errors.  This property is lost if cosh(psip)
                // is replaced by 1/cos(phi), even though it's using "primary" data (phi
                // instead of psip).
                k = Sqrt(_e2m + _e2 * Sq(cphi)) * Hypot(1d, tau)
                  / Hypot(taup, clam);
            }
            else
            {
                xip = PI / 2;
                etap = 0;
                gamma = lon;
                k = _c;
            }
            // {xi',eta'} is {northing,easting} for Gauss-Schreiber transverse Mercator
            // (for eta' = 0, xi' = bet). {xi,eta} is {northing,easting} for transverse
            // Mercator with constant scale on the central meridian (for eta = 0, xip =
            // rectifying latitude).  Define
            //
            //   zeta = xi + i*eta
            //   zeta' = xi' + i*eta'
            //
            // The conversion from conformal to rectifying latitude can be expressed as
            // a series in _n:
            //
            //   zeta = zeta' + sum(h[j-1]' * sin(2 * j * zeta'), j = 1..maxpow_)
            //
            // where h[j]' = O(_n^j).  The reversion of this series gives
            //
            //   zeta' = zeta - sum(h[j-1] * sin(2 * j * zeta), j = 1..maxpow_)
            //
            // which is used in Reverse.
            //
            // Evaluate sums via Clenshaw method.  See
            //    https://en.wikipedia.org/wiki/Clenshaw_algorithm
            //
            // Let
            //
            //    S = sum(a[k] * phi[k](x), k = 0..n)
            //    phi[k+1](x) = alpha[k](x) * phi[k](x) + beta[k](x) * phi[k-1](x)
            //
            // Evaluate S with
            //
            //    b[n+2] = b[n+1] = 0
            //    b[k] = alpha[k](x) * b[k+1] + beta[k+1](x) * b[k+2] + a[k]
            //    S = (a[0] + beta[1](x) * b[2]) * phi[0](x) + b[1] * phi[1](x)
            //
            // Here we have
            //
            //    x = 2 * zeta'
            //    phi[k](x) = sin(k * x)
            //    alpha[k](x) = 2 * cos(x)
            //    beta[k](x) = -1
            //    [ sin(A+B) - 2*cos(B)*sin(A) + sin(A-B) = 0, A = k*x, B = x ]
            //    n = maxpow_
            //    a[k] = _alp[k]
            //    S = b[1] * sin(x)
            //
            // For the derivative we have
            //
            //    x = 2 * zeta'
            //    phi[k](x) = cos(k * x)
            //    alpha[k](x) = 2 * cos(x)
            //    beta[k](x) = -1
            //    [ cos(A+B) - 2*cos(B)*cos(A) + cos(A-B) = 0, A = k*x, B = x ]
            //    a[0] = 1; a[k] = 2*k*_alp[k]
            //    S = (a[0] - b[2]) + b[1] * cos(x)
            //
            // Matrix formulation (not used here):
            //    phi[k](x) = [sin(k * x); k * cos(k * x)]
            //    alpha[k](x) = 2 * [cos(x), 0; -sin(x), cos(x)]
            //    beta[k](x) = -1 * [1, 0; 0, 1]
            //    a[k] = _alp[k] * [1, 0; 0, 1]
            //    b[n+2] = b[n+1] = [0, 0; 0, 0]
            //    b[k] = alpha[k](x) * b[k+1] + beta[k+1](x) * b[k+2] + a[k]
            //    N.B., for all k: b[k](1,2) = 0; b[k](1,1) = b[k](2,2)
            //    S = (a[0] + beta[1](x) * b[2]) * phi[0](x) + b[1] * phi[1](x)
            //    phi[0](x) = [0; 0]
            //    phi[1](x) = [sin(x); cos(x)]
            double
              c0 = Cos(2 * xip), ch0 = Cosh(2 * etap),
              s0 = Sin(2 * xip), sh0 = Sinh(2 * etap);
            var a = new Complex(2 * c0 * ch0, -2 * s0 * sh0); // 2 * cos(2*zeta')
            int n = _maxpow;

            Complex
                y0 = new Complex((n & 1) == 1 ? _alp.Span[n] : 0, 0),
                y1 = new Complex(), // default initializer is 0+i0,
                z0 = new Complex((n & 1) == 1 ? 2 * n * _alp.Span[n] : 0, 0),
                z1 = new Complex();

            if ((n & 1) == 1) --n;

            while (n > 0)
            {
                y1 = a * y0 - y1 + _alp.Span[n];
                z1 = a * z0 - z1 + 2 * n * _alp.Span[n];
                --n;
                y0 = a * y1 - y0 + _alp.Span[n];
                z0 = a * z1 - z0 + 2 * n * _alp.Span[n];
                --n;
            }

            a /= 2d;               // cos(2*zeta')
            z1 = 1d - z1 + a * z0;
            a = new Complex(s0 * ch0, c0 * sh0); // sin(2*zeta')
            y1 = new Complex(xip, etap) + a * y0;

            // Fold in change in convergence and scale for Gauss-Schreiber TM to
            // Gauss-Krueger TM.
            gamma -= Atan2d(z1.Imaginary, z1.Real);
            k *= _b1 * Complex.Abs(z1);
            double xi = y1.Real, eta = y1.Imaginary;

            if (backside)
                gamma = 180 - gamma;
            gamma *= latsign * lonsign;
            gamma = AngNormalize(gamma);
            k *= _k0;

            return (_a1 * _k0 * eta * lonsign, 
                _a1 * _k0 * (backside ? PI - xi : xi) * latsign);
        }

        /// <summary>
        /// Reverse projection, from projected coordinate system to geographic.
        /// </summary>
        /// <param name="lon0">central meridian longitude (degrees).</param>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <param name="gamma">meridian convergence at point (degrees).</param>
        /// <param name="k">azimuthal scale of projection at point; the radial scale is the 1/<paramref name="k"/>.</param>
        /// <returns>
        /// <i>lat</i>, latitude of point (degrees) and <i>lon</i>, longitude of point (degrees).
        /// </returns>
        public (double lat, double lon) Reverse(double lon0, double x, double y, out double gamma, out double k)
        {
            // This undoes the steps in Forward.  The wrinkles are: (1) Use of the
            // reverted series to express zeta' in terms of zeta. (2) Newton's method
            // to solve for phi in terms of tan(phi).
            double
                lat, lon,
              xi = y / (_a1 * _k0),
              eta = x / (_a1 * _k0);

            // Explicitly enforce the parity
            int
              xisign = (xi < 0) ? -1 : 1,
              etasign = (eta < 0) ? -1 : 1;

            xi *= xisign;
            eta *= etasign;
            bool backside = xi > PI / 2;
            if (backside)
                xi = PI - xi;

            double
              c0 = Cos(2 * xi), ch0 = Cosh(2 * eta),
              s0 = Sin(2 * xi), sh0 = Sinh(2 * eta);

            var a = new Complex(2 * c0 * ch0, -2 * s0 * sh0); // 2 * cos(2*zeta)
            int n = _maxpow;
            Complex
                y0 = new Complex((n & 1) == 1 ? -_bet.Span[n] : 0, 0),
                y1 = new Complex(), // default initializer is 0+i0
                z0 = new Complex((n & 1) == 1 ? -2 * n * _bet.Span[n] : 0, 0),
                z1 = new Complex();

            if ((n & 1) == 1) --n;

            while (n > 0)
            {
                y1 = a * y0 - y1 - _bet.Span[n];
                z1 = a * z0 - z1 - 2 * n * _bet.Span[n];
                --n;
                y0 = a * y1 - y0 - _bet.Span[n];
                z0 = a * z1 - z0 - 2 * n * _bet.Span[n];
                --n;
            }
            a /= 2d;               // cos(2*zeta)
            z1 = 1d - z1 + a * z0;
            a = new Complex(s0 * ch0, c0 * sh0); // sin(2*zeta)
            y1 = new Complex(xi, eta) + a * y0;
            // Convergence and scale for Gauss-Schreiber TM to Gauss-Krueger TM.
            gamma = Atan2d(z1.Imaginary, z1.Real);
            k = _b1 / Complex.Abs(z1);
            // JHS 154 has
            //
            //   phi' = asin(sin(xi') / cosh(eta')) (Krueger p 17 (25))
            //   lam = asin(tanh(eta') / cos(phi')
            //   psi = asinh(tan(phi'))
            double
              xip = y1.Real, etap = y1.Imaginary,
              s = Sinh(etap),
              c = Max(0d, Cos(xip)), // cos(pi/2) might be negative
              r = Hypot(s, c);

            if (r != 0)
            {
                lon = Atan2d(s, c); // Krueger p 17 (25)
                                    // Use Newton's method to solve for tau
                double
                  sxip = Sin(xip),
                  tau = Tauf(sxip / r, _es);

                gamma += Atan2d(sxip * Tanh(etap), c); // Krueger p 19 (31)
                lat = Atand(tau);
                // Note cos(phi') * cosh(eta') = r
                k *= Sqrt(_e2m + _e2 / (1 + Sq(tau))) *
                  Hypot(1d, tau) * r;
            }
            else
            {
                lat = 90;
                lon = 0;
                k *= _c;
            }
            lat *= xisign;
            if (backside)
                lon = 180 - lon;
            lon *= etasign;
            lon = AngNormalize(lon + lon0);
            if (backside)
                gamma = 180 - gamma;
            gamma *= xisign * etasign;
            gamma = AngNormalize(gamma);
            k *= _k0;

            return (lat, lon);
        }

        /// <summary>
        /// Forward without returning convergence and scale.
        /// </summary>
        /// <param name="lon0">central meridian longitude (degrees).</param>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <returns>
        /// <i>x</i>, easting of point (meters) and <i>y</i>, northing of point (meters).
        /// </returns>
        public (double x, double y) Forward(double lon0, double lat, double lon) => Forward(lon0, lat, lon, out _, out _);

        /// <summary>
        /// Reverse without returning convergence and scale.
        /// </summary>
        /// <param name="lon0">central meridian longitude (degrees).</param>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <returns>
        /// <i>lat</i>, latitude of point (degrees) and <i>lon</i>, longitude of point (degrees).
        /// </returns>
        public (double lat, double lon) Reverse(double lon0, double x, double y) => Reverse(lon0, x, y, out _, out _);
    }
}
