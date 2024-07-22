using System;
using static GeographicLib.Macros;
using static GeographicLib.MathEx;
using static System.Math;

namespace GeographicLib
{
    /// <summary>
    /// Conversions between auxiliary latitudes.
    /// </summary>
    /// <remarks>
    /// This class is an implementation of the methods described in
    /// 
    /// <list type="bullet">
    /// <item>
    /// C. F. F. Karney,
    /// <a href = "https://doi.org/10.1080/00396265.2023.2217604">On auxiliary latitudes</a>,
    /// Survey Review 56(395), 165--180 (2024);
    /// preprint <a href = "https://arxiv.org/abs/2212.05818">arXiv:2212.05818</a>.
    /// </item>
    /// </list>
    /// 
    /// <para>
    /// The class provides accurate conversions between geographic(<i>phi</i>, φ),
    /// parametric(<i>beta</i>, β), geocentric(<i>theta</i>, θ), rectifying
    /// (<i>mu</i>, µ), conformal(<i>chi</i>, χ), and authalic(<i>xi</i>, ξ)
    /// latitudes for an ellipsoid of revolution. A latitude is represented by
    /// <see cref="AuxAngle"/> in order to maintain precision close to the poles.
    /// </para>
    ///
    /// The class implements two methods for the conversion:
    /// <list type="bullet">
    /// <item>
    /// Direct evaluation of the defining equations, the <b>exact</b> method. These
    /// equations are formulated so as to preserve relative accuracy of the
    /// tangent of the latitude, ensuring high accuracy near the equator and the
    /// poles. Newton's method is used for those conversions that can't be
    /// expressed in closed form.
    /// </item>
    /// <item>
    /// Expansions in powers of <i>n</i>, the third flattening, the <b>series</b> method.
    /// This delivers full accuracy for abs(<i>f</i>) ≤ 1/150. Here, <i>f</i> is the
    /// flattening of the ellipsoid.
    /// </item>
    /// </list>
    ///
    /// The series method is the preferred method of conversion for any conversion
    /// involving µ, χ, or ξ, with abs(<i>f</i>) ≤ 1/150.  The equations
    /// for the conversions between φ, β, and θ are sufficiently
    /// simple that the exact method should be used for such conversions and also
    /// for conversions with with abs(<i>f</i>) &gt; 1/150.
    /// </remarks>
    public class AuxLatitude : IEllipsoid
    {
        /// <summary>
        /// The order of the series expansions.  This is set at compile time to
        /// either 4, 6, or 8, by the preprocessor macro
        /// GEOGRAPHICLIB_AUXLATITUDE_ORDER.
        /// </summary>
        internal const int Lmax = GEOGRAPHICLIB_AUXLATITUDE_ORDER;

        /// <summary>
        /// The total number of auxiliary latitudes.
        /// </summary>
        private const int AUXNUMBER = 6;

        // Maximum number of iterations for Newton's method
        private const int numit_ = 1000;

        // Static consts for Newton's method
        // the function atanh(e * sphi)/e + sphi / (1 - (e * sphi)^2);
        private readonly double tol_, bmin_, bmax_;

        // Ellipsoid parameters
        internal readonly double _a, _b, _f, _fm1, _e2, _e2m1, _e12, _e12p1, _n, _e, _e1, _n2, _q;

        // To hold computed Fourier coefficients
        internal readonly double[] _c = new double[Lmax * AUXNUMBER * AUXNUMBER];

        private static readonly double[] s_coeffs = new double[]
        {
          // C[phi,phi] skipped
          // C[phi,beta]; even coeffs only
          0, 0, 1,
          0, 0, 1/(double)(2),
          0, 1/(double)(3),
          0, 1/(double)(4),
          1/(double)(5),
          1/(double)(6),
          // C[phi,theta]; even coeffs only
          2, -2, 2,
          6, -4, 2,
          -8, 8/(double)(3),
          -16, 4,
          32/(double)(5),
          32/(double)(3),
          // C[phi,mu]; even coeffs only
          269/(double)(512), -27/(double)(32), 3/(double)(2),
          6759/(double)(4096), -55/(double)(32), 21/(double)(16),
          -417/(double)(128), 151/(double)(96),
          -15543/(double)(2560), 1097/(double)(512),
          8011/(double)(2560),
          293393/(double)(61440),
          // C[phi,chi]
          -2854/(double)(675), 26/(double)(45), 116/(double)(45), -2, -2/(double)(3), 2,
          2323/(double)(945), 2704/(double)(315), -227/(double)(45), -8/(double)(5), 7/(double)(3),
          73814/(double)(2835), -1262/(double)(105), -136/(double)(35), 56/(double)(15),
          -399572/(double)(14175), -332/(double)(35), 4279/(double)(630),
          -144838/(double)(6237), 4174/(double)(315),
          601676/(double)(22275),
          // C[phi,xi]
          28112932/(double)(212837625), 60136/(double)(467775), -2582/(double)(14175),
          -16/(double)(35), 4/(double)(45), 4/(double)(3),
          251310128/(double)(638512875), -21016/(double)(51975), -11966/(double)(14175),
          152/(double)(945), 46/(double)(45),
          -8797648/(double)(10945935), -94388/(double)(66825), 3802/(double)(14175),
          3044/(double)(2835),
          -1472637812/(double)(638512875), 41072/(double)(93555), 6059/(double)(4725),
          455935736/(double)(638512875), 768272/(double)(467775),
          4210684958L/(double)(1915538625),
          // C[beta,phi]; even coeffs only
          0, 0, -1,
          0, 0, 1/(double)(2),
          0, -1/(double)(3),
          0, 1/(double)(4),
          -1/(double)(5),
          1/(double)(6),
          // C[beta,beta] skipped
          // C[beta,theta]; even coeffs only
          0, 0, 1,
          0, 0, 1/(double)(2),
          0, 1/(double)(3),
          0, 1/(double)(4),
          1/(double)(5),
          1/(double)(6),
          // C[beta,mu]; even coeffs only
          205/(double)(1536), -9/(double)(32), 1/(double)(2),
          1335/(double)(4096), -37/(double)(96), 5/(double)(16),
          -75/(double)(128), 29/(double)(96),
          -2391/(double)(2560), 539/(double)(1536),
          3467/(double)(7680),
          38081/(double)(61440),
          // C[beta,chi]
          -3118/(double)(4725), -1/(double)(3), 38/(double)(45), -1/(double)(3), -2/(double)(3), 1,
          -247/(double)(270), 50/(double)(21), -7/(double)(9), -14/(double)(15), 5/(double)(6),
          17564/(double)(2835), -5/(double)(3), -34/(double)(21), 16/(double)(15),
          -49877/(double)(14175), -28/(double)(9), 2069/(double)(1260),
          -28244/(double)(4455), 883/(double)(315),
          797222/(double)(155925),
          // C[beta,xi]
          7947332/(double)(212837625), 11824/(double)(467775), -1082/(double)(14175),
          -46/(double)(315), 4/(double)(45), 1/(double)(3),
          39946703/(double)(638512875), -16672/(double)(155925), -338/(double)(2025),
          68/(double)(945), 17/(double)(90),
          -255454/(double)(1563705), -101069/(double)(467775), 1102/(double)(14175),
          461/(double)(2835),
          -189032762/(double)(638512875), 1786/(double)(18711), 3161/(double)(18900),
          80274086/(double)(638512875), 88868/(double)(467775),
          880980241/(double)(3831077250L),
          // C[theta,phi]; even coeffs only
          -2, 2, -2,
          6, -4, 2,
          8, -8/(double)(3),
          -16, 4,
          -32/(double)(5),
          32/(double)(3),
          // C[theta,beta]; even coeffs only
          0, 0, -1,
          0, 0, 1/(double)(2),
          0, -1/(double)(3),
          0, 1/(double)(4),
          -1/(double)(5),
          1/(double)(6),
          // C[theta,theta] skipped
          // C[theta,mu]; even coeffs only
          499/(double)(1536), -23/(double)(32), -1/(double)(2),
          6565/(double)(12288), -5/(double)(96), 5/(double)(16),
          -77/(double)(128), 1/(double)(32),
          -4037/(double)(7680), 283/(double)(1536),
          1301/(double)(7680),
          17089/(double)(61440),
          // C[theta,chi]
          -3658/(double)(4725), 2/(double)(9), 4/(double)(9), -2/(double)(3), -2/(double)(3), 0,
          61/(double)(135), 68/(double)(45), -23/(double)(45), -4/(double)(15), 1/(double)(3),
          9446/(double)(2835), -46/(double)(35), -24/(double)(35), 2/(double)(5),
          -34712/(double)(14175), -80/(double)(63), 83/(double)(126),
          -2362/(double)(891), 52/(double)(45),
          335882/(double)(155925),
          // C[theta,xi]
          216932/(double)(2627625), 109042/(double)(467775), -2102/(double)(14175),
          -158/(double)(315), 4/(double)(45), -2/(double)(3),
          117952358/(double)(638512875), -7256/(double)(155925), 934/(double)(14175),
          -16/(double)(945), 16/(double)(45),
          -7391576/(double)(54729675), -25286/(double)(66825), 922/(double)(14175),
          -232/(double)(2835),
          -67048172/(double)(638512875), 268/(double)(18711), 719/(double)(4725),
          46774256/(double)(638512875), 14354/(double)(467775),
          253129538/(double)(1915538625),
          // C[mu,phi]; even coeffs only
          -3/(double)(32), 9/(double)(16), -3/(double)(2),
          135/(double)(2048), -15/(double)(32), 15/(double)(16),
          105/(double)(256), -35/(double)(48),
          -189/(double)(512), 315/(double)(512),
          -693/(double)(1280),
          1001/(double)(2048),
          // C[mu,beta]; even coeffs only
          -1/(double)(32), 3/(double)(16), -1/(double)(2),
          -9/(double)(2048), 1/(double)(32), -1/(double)(16),
          3/(double)(256), -1/(double)(48),
          3/(double)(512), -5/(double)(512),
          -7/(double)(1280),
          -7/(double)(2048),
          // C[mu,theta]; even coeffs only
          -15/(double)(32), 13/(double)(16), 1/(double)(2),
          -1673/(double)(2048), 33/(double)(32), -1/(double)(16),
          349/(double)(256), -5/(double)(16),
          963/(double)(512), -261/(double)(512),
          -921/(double)(1280),
          -6037/(double)(6144),
          // C[mu,mu] skipped
          // C[mu,chi]
          7891/(double)(37800), -127/(double)(288), 41/(double)(180), 5/(double)(16), -2/(double)(3),
          1/(double)(2),
          -1983433/(double)(1935360), 281/(double)(630), 557/(double)(1440), -3/(double)(5),
          13/(double)(48),
          167603/(double)(181440), 15061/(double)(26880), -103/(double)(140), 61/(double)(240),
          6601661/(double)(7257600), -179/(double)(168), 49561/(double)(161280),
          -3418889/(double)(1995840), 34729/(double)(80640),
          212378941/(double)(319334400),
          // C[mu,xi]
          12674323/(double)(851350500), -384229/(double)(14968800), -1609/(double)(28350),
          121/(double)(1680), 4/(double)(45), -1/(double)(6),
          -31621753811L/(double)(1307674368000L), -431/(double)(17325),
          16463/(double)(453600), 26/(double)(945), -29/(double)(720),
          -32844781/(double)(1751349600), 3746047/(double)(119750400), 449/(double)(28350),
          -1003/(double)(45360),
          10650637121L/(double)(326918592000L), 629/(double)(53460),
          -40457/(double)(2419200),
          205072597/(double)(20432412000L), -1800439/(double)(119750400),
          -59109051671L/(double)(3923023104000L),
          // C[chi,phi]
          4642/(double)(4725), 32/(double)(45), -82/(double)(45), 4/(double)(3), 2/(double)(3), -2,
          -1522/(double)(945), 904/(double)(315), -13/(double)(9), -16/(double)(15), 5/(double)(3),
          -12686/(double)(2835), 8/(double)(5), 34/(double)(21), -26/(double)(15),
          -24832/(double)(14175), -12/(double)(5), 1237/(double)(630),
          109598/(double)(31185), -734/(double)(315),
          444337/(double)(155925),
          // C[chi,beta]
          -998/(double)(4725), 2/(double)(5), -16/(double)(45), 0, 2/(double)(3), -1,
          -2/(double)(27), -22/(double)(105), 19/(double)(45), -2/(double)(5), 1/(double)(6),
          116/(double)(567), -22/(double)(105), 16/(double)(105), -1/(double)(15),
          2123/(double)(14175), -8/(double)(105), 17/(double)(1260),
          128/(double)(4455), -1/(double)(105),
          149/(double)(311850),
          // C[chi,theta]
          1042/(double)(4725), -14/(double)(45), -2/(double)(9), 2/(double)(3), 2/(double)(3), 0,
          -712/(double)(945), -4/(double)(45), 43/(double)(45), 4/(double)(15), -1/(double)(3),
          274/(double)(2835), 124/(double)(105), 2/(double)(105), -2/(double)(5),
          21068/(double)(14175), -16/(double)(105), -55/(double)(126),
          -9202/(double)(31185), -22/(double)(45),
          -90263/(double)(155925),
          // C[chi,mu]
          -96199/(double)(604800), 81/(double)(512), 1/(double)(360), -37/(double)(96), 2/(double)(3),
          -1/(double)(2),
          1118711/(double)(3870720), -46/(double)(105), 437/(double)(1440), -1/(double)(15),
          -1/(double)(48),
          -5569/(double)(90720), 209/(double)(4480), 37/(double)(840), -17/(double)(480),
          830251/(double)(7257600), 11/(double)(504), -4397/(double)(161280),
          108847/(double)(3991680), -4583/(double)(161280),
          -20648693/(double)(638668800),
          // C[chi,chi] skipped
          // C[chi,xi]
          -55271278/(double)(212837625), 27128/(double)(93555), -2312/(double)(14175),
          -88/(double)(315), 34/(double)(45), -2/(double)(3),
          106691108/(double)(638512875), -65864/(double)(155925), 6079/(double)(14175),
          -184/(double)(945), 1/(double)(45),
          5921152/(double)(54729675), -14246/(double)(467775), 772/(double)(14175),
          -106/(double)(2835),
          75594328/(double)(638512875), -5312/(double)(467775), -167/(double)(9450),
          2837636/(double)(638512875), -248/(double)(13365),
          -34761247/(double)(1915538625),
          // C[xi,phi]
          -44732/(double)(2837835), 20824/(double)(467775), 538/(double)(4725), 88/(double)(315),
          -4/(double)(45), -4/(double)(3),
          -12467764/(double)(212837625), -37192/(double)(467775), -2482/(double)(14175),
          8/(double)(105), 34/(double)(45),
          100320856/(double)(1915538625), 54968/(double)(467775), -898/(double)(14175),
          -1532/(double)(2835),
          -5884124/(double)(70945875), 24496/(double)(467775), 6007/(double)(14175),
          -839792/(double)(19348875), -23356/(double)(66825),
          570284222/(double)(1915538625),
          // C[xi,beta]
          -70496/(double)(8513505), 2476/(double)(467775), 34/(double)(675), 32/(double)(315),
          -4/(double)(45), -1/(double)(3),
          53836/(double)(212837625), 3992/(double)(467775), 74/(double)(2025), -4/(double)(315),
          -7/(double)(90),
          -661844/(double)(1915538625), 7052/(double)(467775), 2/(double)(14175),
          -83/(double)(2835),
          1425778/(double)(212837625), 934/(double)(467775), -797/(double)(56700),
          390088/(double)(212837625), -3673/(double)(467775),
          -18623681/(double)(3831077250L),
          // C[xi,theta]
          -4286228/(double)(42567525), -193082/(double)(467775), 778/(double)(4725),
          62/(double)(105), -4/(double)(45), 2/(double)(3),
          -61623938/(double)(70945875), 92696/(double)(467775), 12338/(double)(14175),
          -32/(double)(315), 4/(double)(45),
          427003576/(double)(1915538625), 612536/(double)(467775), -1618/(double)(14175),
          -524/(double)(2835),
          427770788/(double)(212837625), -8324/(double)(66825), -5933/(double)(14175),
          -9153184/(double)(70945875), -320044/(double)(467775),
          -1978771378/(double)(1915538625),
          // C[xi,mu]
          -9292991/(double)(302702400), 7764059/(double)(239500800), 1297/(double)(18900),
          -817/(double)(10080), -4/(double)(45), 1/(double)(6),
          36019108271L/(double)(871782912000L), 35474/(double)(467775),
          -29609/(double)(453600), -2/(double)(35), 49/(double)(720),
          3026004511L/(double)(30648618000L), -4306823/(double)(59875200),
          -2917/(double)(56700), 4463/(double)(90720),
          -368661577/(double)(4036032000L), -102293/(double)(1871100),
          331799/(double)(7257600),
          -875457073/(double)(13621608000L), 11744233/(double)(239500800),
          453002260127L/(double)(7846046208000L),
          // C[xi,chi]
          2706758/(double)(42567525), -55222/(double)(93555), 2458/(double)(4725),
          46/(double)(315), -34/(double)(45), 2/(double)(3),
          -340492279/(double)(212837625), 516944/(double)(467775), 3413/(double)(14175),
          -256/(double)(315), 19/(double)(45),
          4430783356L/(double)(1915538625), 206834/(double)(467775), -15958/(double)(14175),
          248/(double)(567),
          62016436/(double)(70945875), -832976/(double)(467775), 16049/(double)(28350),
          -651151712/(double)(212837625), 15602/(double)(18711),
          2561772812L/(double)(1915538625),
          // C[xi,xi] skipped
        };

        private static readonly int[] s_ptrs = new int[]{
          0, 0, 12, 24, 36, 57, 78, 90, 90, 102, 114, 135, 156, 168, 180, 180, 192,
          213, 234, 246, 258, 270, 270, 291, 312, 333, 354, 375, 396, 396, 417,
          438, 459, 480, 501, 522, 522,
        };

        /// <summary>
        /// Initializes a new instance of the <see cref="AuxLatitude"/> class
        /// with the equatorial radius and flattening. 
        /// </summary>
        /// <param name="a">Equatorial radius.</param>
        /// <param name="f">
        /// Flattening of ellipsoid. Setting <paramref name="f"/> = 0 gives a sphere.
        /// Negative <paramref name="f"/> gives a prolate ellipsoid.
        /// </param>
        /// <remarks>
        /// The constructor does not precompute the coefficients for the
        /// Fourier series for the series conversions. These are computed and saved
        /// when first needed.
        /// </remarks>
        /// <exception cref="GeographicException"/>
        public AuxLatitude(double a, double f)
        {
            tol_ = Sqrt(DBL_EPSILON);
            bmin_ = Log2(double.MinValue);
            bmax_ = Log2(double.MaxValue);
            _a = a;
            _b = _a * (1 - f);
            _f = f;
            _fm1 = 1 - _f;
            _e2 = _f * (2 - _f);
            _e2m1 = _fm1 * _fm1;
            _e12 = _e2 / (1 - _e2);
            _e12p1 = 1 / _e2m1;
            _n = _f / (2 - _f);
            _e = Sqrt(Abs(_e2));
            _e1 = Sqrt(Abs(_e12));
            _n2 = _n * _n;
            _q = _e12p1 + (_f == 0 ? 1 : (_f > 0 ? Asinh(_e1) : Atan(_e)) / _e);

            if (!(IsFinite(_a) && _a > 0))
                throw new GeographicException("Equatorial radius is not positive");
            if (!(IsFinite(_b) && _b > 0))
                throw new GeographicException("Polar semi-axis is not positive");

            _c.AsSpan().Fill(double.NaN);
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="AuxLatitude"/> class
        /// with the given <see cref="IEllipsoid"/> object. 
        /// </summary>
        /// <param name="ellipsoid">Source <see cref="IEllipsoid"/> object.</param>
        public AuxLatitude(IEllipsoid ellipsoid) : this(ellipsoid.EquatorialRadius, ellipsoid.Flattening)
        {
        }

        private AuxLatitude((double a, double b) axes)
        {
            tol_ = Sqrt(DBL_EPSILON);
            bmin_ = Log2(double.MinValue);
            bmax_ = Log2(double.MaxValue);
            _a = axes.a;
            _b = axes.b;
            _f = (_a - _b) / _a;
            _fm1 = _b / _a;
            _e2 = ((_a - _b) * (_a + _b)) / (_a * _a);
            _e2m1 = (_b * _b) / (_a * _a);
            _e12 = ((_a - _b) * (_a + _b)) / (_b * _b);
            _e12p1 = (_a * _a) / (_b * _b);
            _n = (_a - _b) / (_a + _b);
            _e = Sqrt(Abs(_a - _b) * (_a + _b)) / _a;
            _e1 = Sqrt(Abs(_a - _b) * (_a + _b)) / _b;
            _n2 = _n * _n;
            _q = _e12p1 + (_f == 0 ? 1 : (_f > 0 ? Asinh(_e1) : Atan(_e)) / _e);

            if (!(IsFinite(_a) && _a > 0))
                throw new GeographicException("Equatorial radius is not positive");
            if (!(IsFinite(_b) && _b > 0))
                throw new GeographicException("Polar semi-axis is not positive");

            _c.AsSpan().Fill(double.NaN);
        }

        /// <summary>
        /// A global instantiation of <see cref="AuxLatitude"/> with the parameters for the WGS84 ellipsoid.
        /// </summary>
        public static readonly AuxLatitude WGS84 = new AuxLatitude(Constants.WGS84_a, Constants.WGS84_f);

        /// <summary>
        /// Construct and return an AuxLatitude object specified in terms of the semi-axes.
        /// </summary>
        /// <param name="a">Equatorial radius.</param>
        /// <param name="b">Polar semi-axis.</param>
        /// <returns></returns>
        /// <exception cref="GeographicException"/>
        public static AuxLatitude FromAxes(double a, double b)
        {
            return new AuxLatitude((a, b));
        }

        /// <summary>
        /// Use Clenshaw to sum a Fouier series.
        /// </summary>
        /// <param name="sinp">If <see langword="true"/> sum the sine series, else sum the cosine series.</param>
        /// <param name="szeta">sin(<i>zeta</i>).</param>
        /// <param name="czeta">cos(<i>zeta</i>).</param>
        /// <param name="c">The array of coefficients.</param>
        /// <param name="K">The number of coefficients.</param>
        /// <returns>
        /// sum(c[k] * sin((2*k+2) * zeta), i, 0, K-1);, 
        /// if <paramref name="sinp"/> is <see langword="true"/>;
        /// with <paramref name="sinp"/> <see langword="false"/>, replace sin by cos.
        /// </returns>
        public static double Clenshaw(bool sinp, double szeta, double czeta,
                         ReadOnlySpan<double> c, int K)
        {
            // Evaluate
            // y = sum(c[k] * sin( (2*k+2) * zeta), i, 0, K-1) if  sinp
            // y = sum(c[k] * cos( (2*k+2) * zeta), i, 0, K-1) if !sinp
            // Approx operation count = (K + 5) mult and (2 * K + 2) add
            int k = K;
            double u0 = 0, u1 = 0,        // accumulators for sum
              x = 2 * (czeta - szeta) * (czeta + szeta); // 2 * cos(2*zeta)
            for (; k > 0;)
            {
                double t = x * u0 - u1 + c[--k];
                u1 = u0; u0 = t;
            }
            // u0*f0(zeta) - u1*fm1(zeta)
            // f0 = sinp ? sin(2*zeta) : cos(2*zeta)
            // fm1 = sinp ? 0 : 1
            double f0 = sinp ? 2 * szeta * czeta : x / 2, fm1 = sinp ? 0 : 1;
            return f0 * u0 - fm1 * u1;
        }

        /// <summary>
        /// Convert between any two auxiliary latitudes specified as <see cref="AuxAngle"/>.
        /// </summary>
        /// <param name="auxin">Type of auxiliary latitude <i>zeta</i>.</param>
        /// <param name="auxout">Type of auxiliary latitude <i>eta</i>.</param>
        /// <param name="zeta">The input auxiliary latitude as an <see cref="AuxAngle"/>.</param>
        /// <param name="exact">If <see langword="true"/> use the exact equations instead of the Taylor series.</param>
        /// <returns>The output auxiliary latitude <i>eta</i> as an <see cref="AuxAngle"/>.</returns>
        /// <remarks>
        /// With <paramref name="exact"/> = <see langword="false"/>, the Fourier coefficients for a specific <paramref name="auxin"/>
        /// and <paramref name="auxout"/> are computed and saved on the first call; the saved
        /// coefficients are used on subsequent calls. The series method is
        /// accurate for abs(<i>f</i>) ≤ 1/150; for other <i>f</i>, the exact method
        /// should be used.
        /// </remarks>
        public AuxAngle Convert(
            AuxLatitudeType auxin,
            AuxLatitudeType auxout,
            AuxAngle zeta,
            bool exact = false)
        {
            int k = ind(auxout, auxin);
            if (k < 0) return AuxAngle.NaN;
            if (auxin == auxout) return zeta;
            if (exact)
            {
                if ((int)auxin < 3 && (int)auxout < 3)
                    // Need extra real because, since C++11, pow(float, int) returns double
                    return new AuxAngle(zeta.Y * Pow(_fm1, auxout - auxin), zeta.X);
                else
                    return ToAuxiliary(auxout, FromAuxiliary(auxin, zeta));
            }
            else
            {
                if (double.IsNaN(_c[Lmax * (k + 1) - 1])) fillcoeff(auxin, auxout, k);
                AuxAngle zetan = zeta.Normalized();
                double d = Clenshaw(true, zetan.Y, zetan.X, _c.AsSpan().Slice(Lmax * k), Lmax);
                zetan += AuxAngle.FromRadians(d);
                return zetan;
            }
        }

        /// <summary>
        /// Convert between any two auxiliary latitudes specified in degrees.
        /// </summary>
        /// <param name="auxin">Type of auxiliary latitude <i>zeta</i>.</param>
        /// <param name="auxout">Type of auxiliary latitude <i>eta</i>.</param>
        /// <param name="zeta">The input auxiliary latitude in degrees.</param>
        /// <param name="exact">If <see langword="true"/> use the exact equations instead of the Taylor series.</param>
        /// <returns>The output auxiliary latitude <i>eta</i> in degrees.</returns>
        /// <remarks>
        /// With <paramref name="exact"/> = <see langword="false"/>, the Fourier coefficients for a specific <paramref name="auxin"/>
        /// and <paramref name="auxout"/> are computed and saved on the first call; the saved
        /// coefficients are used on subsequent calls. The series method is
        /// accurate for abs(<i>f</i>) ≤ 1/150; for other <i>f</i>, the exact method
        /// should be used.
        /// </remarks>
        public double Convert(
            AuxLatitudeType auxin,
            AuxLatitudeType auxout,
            double zeta,
            bool exact = false)
        {
            AuxAngle zetaa = AuxAngle.FromDegrees(zeta);
            double m = Round((zeta - zetaa.Degrees) / TD);
            return TD * m + Convert(auxin, auxout, zetaa, exact).Degrees;
        }

        private unsafe AuxAngle ToAuxiliary(AuxLatitudeType auxout, AuxAngle phi, double* diff)
        {
            switch (auxout)
            {
                case AuxLatitudeType.Geographic: if (diff != null) *diff = 1; return phi;
                case AuxLatitudeType.Parametric: return Parametric(phi, diff);
                case AuxLatitudeType.Geocentric: return Geocentric(phi, diff);
                case AuxLatitudeType.Rectifying: return Rectifying(phi, diff);
                case AuxLatitudeType.Conformal: return Conformal(phi, diff);
                case AuxLatitudeType.Authalic: return Authalic(phi, diff);
                default:
                    if (diff != null) *diff = double.NaN;
                    return AuxAngle.NaN;
            }
        }

        /// <summary>
        /// Convert geographic latitude to an auxiliary latitude <i>eta</i>.
        /// </summary>
        /// <param name="auxout">The auxiliary latitude to be returned.</param>
        /// <param name="phi">The geographic latitude.</param>
        /// <param name="diff">the derivative d tan(<i>eta</i>) / d tan(<i>phi</i>).</param>
        /// <returns>The auxiliary latitude <i>eta</i>.</returns>
        /// <remarks>
        /// This uses the exact equations.
        /// </remarks>
        public unsafe AuxAngle ToAuxiliary(AuxLatitudeType auxout, AuxAngle phi, out double diff)
        {
            diff = 0;
            fixed (double* pdiff = &diff)
            {
                return ToAuxiliary(auxout, phi, pdiff);
            }
        }

        /// <summary>
        /// Convert geographic latitude to an auxiliary latitude <i>eta</i>.
        /// </summary>
        /// <param name="auxout">The auxiliary latitude to be returned.</param>
        /// <param name="phi">The geographic latitude.</param>
        /// <returns>The auxiliary latitude <i>eta</i>.</returns>
        /// <remarks>
        /// This uses the exact equations.
        /// </remarks>
        public unsafe AuxAngle ToAuxiliary(AuxLatitudeType auxout, AuxAngle phi)
        {
            return ToAuxiliary(auxout, phi, null);
        }

        private unsafe AuxAngle FromAuxiliary(AuxLatitudeType auxin, AuxAngle zeta,
                           int* niter)
        {
            int n = 0; if (niter != null) *niter = n;
            double tphi = _fm1;
            switch (auxin)
            {
                case AuxLatitudeType.Geographic: return zeta;
                // case PARAMETRIC:                   break;
                case AuxLatitudeType.Parametric: return new AuxAngle(zeta.Y / _fm1, zeta.X);
                // case GEOCENTRIC: tphi *= _fm1  ; break;
                case AuxLatitudeType.Geocentric: return new AuxAngle(zeta.Y / _e2m1, zeta.X);
                case AuxLatitudeType.Rectifying: tphi *= Sqrt(_fm1); break;
                case AuxLatitudeType.Conformal: tphi *= _fm1; break;
                case AuxLatitudeType.Authalic: tphi *= Cbrt(_fm1); break;
                default: return AuxAngle.NaN;
            }

            // Drop through to solution by Newton's method
            double tzeta = Abs(zeta.Tan), ltzeta = Log2(tzeta);
            if (!IsFinite(ltzeta)) return zeta;
            tphi = tzeta / tphi;
            double ltphi = Log2(tphi),
              bmin = Min(ltphi, bmin_), bmax = Max(ltphi, bmax_);
            for (int sign = 0, osign = 0, ntrip = 0; n < numit_;)
            {
                ++n;
                double diff;
                AuxAngle zeta1 = ToAuxiliary(auxin, new AuxAngle(tphi), &diff);
                double tzeta1 = zeta1.Tan, ltzeta1 = Log2(tzeta1);
                // Convert derivative from dtan(zeta)/dtan(phi) to
                // dlog(tan(zeta))/dlog(tan(phi))
                diff *= tphi / tzeta1;
                osign = sign;
                if (tzeta1 == tzeta)
                    break;
                else if (tzeta1 > tzeta)
                {
                    sign = 1;
                    bmax = ltphi;
                }
                else
                {
                    sign = -1;
                    bmin = ltphi;
                }
                double dltphi = -(ltzeta1 - ltzeta) / diff;
                ltphi += dltphi;
                tphi = Exp2(ltphi);
                if (!(Abs(dltphi) >= tol_))
                {
                    ++n;
                    // Final Newton iteration without the logs
                    zeta1 = ToAuxiliary(auxin, new AuxAngle(tphi), &diff);
                    tphi -= (zeta1.Tan - tzeta) / diff;
                    break;
                }
                if ((sign * osign < 0 && n - ntrip > 2) ||
                    ltphi >= bmax || ltphi <= bmin)
                {
                    sign = 0; ntrip = n;
                    ltphi = (bmin + bmax) / 2;
                    tphi = Exp2(ltphi);
                }
            }
            if (niter != null) *niter = n;
            return new AuxAngle(tphi).CopyQuadrant(zeta);
        }

        /// <summary>
        /// Convert an auxiliary latitude <paramref name="zeta"/> to geographic latitude.
        /// </summary>
        /// <param name="auxin">The type of auxiliary latitude <paramref name="zeta"/>.</param>
        /// <param name="zeta">The input auxiliary latitude.</param>
        /// <param name="niter">The number of iterations.</param>
        /// <returns>The geographic latitude <i>phi</i>.</returns>
        /// <remarks>
        /// This uses the exact equations.
        /// </remarks>
        public unsafe AuxAngle FromAuxiliary(AuxLatitudeType auxin, AuxAngle zeta,
                           out int niter)
        {
            niter = 0;
            fixed (int* p = &niter)
            {
                return FromAuxiliary(auxin, zeta, p);
            }
        }

        /// <summary>
        /// Convert an auxiliary latitude <paramref name="zeta"/> to geographic latitude.
        /// </summary>
        /// <param name="auxin">The type of auxiliary latitude <paramref name="zeta"/>.</param>
        /// <param name="zeta">The input auxiliary latitude.</param>
        /// <returns>The geographic latitude <i>phi</i>.</returns>
        /// <remarks>
        /// This uses the exact equations.
        /// </remarks>
        public unsafe AuxAngle FromAuxiliary(AuxLatitudeType auxin, AuxAngle zeta)
        {
            return FromAuxiliary(auxin, zeta, null);
        }

        /// <summary>
        /// Return the rectifying radius.
        /// </summary>
        /// <param name="exact">If <see langword="true"/> use the exact expression instead of the Taylor series.</param>
        /// <returns>The rectifying radius in the same units as <i>a</i>.</returns>
        public double RectifyingRadius(bool exact = false)
        {
            if (exact)
            {
                return EllipticFunction.RG(Sq(_a), Sq(_b)) * 4 / PI;
            }
            else
            {
                // Maxima code for these coefficients:
                // df[i]:=if i<0 then df[i+2]/(i+2) else i!!$
                // R(Lmax):=sum((df[2*j-3]/df[2*j])^2*n^(2*j),j,0,floor(Lmax/2))$
                // cf(Lmax):=block([t:R(Lmax)],
                //  t:makelist(coeff(t,n,2*(floor(Lmax/2)-j)),j,0,floor(Lmax/2)),
                //  map(lambda([x],num(x)/
                //         (if denom(x) = 1 then 1 else real(denom(x)))),t))$

                // GEOGRAPHICLIB_AUXLATITUDE_ORDER == 6
                Span<double> coeff = stackalloc double[] { 1 / 256.0, 1 / 64.0, 1 / 4.0, 1 };

                int m = Lmax / 2;
                return (_a + _b) / 2 * PolyVal(m, coeff, _n2);
            }
        }

        /// <summary>
        /// Return the authalic radius squared.
        /// </summary>
        /// <param name="exact">If <see langword="true"/> use the exact expression instead of the Taylor series.</param>
        /// <returns>The authalic radius squared in the same units as <i>a</i>.</returns>
        public double AuthalicRadiusSquared(bool exact = false)
        {
            if (exact)
            {
                return Sq(_b) * _q / 2;
            }
            else
            {
                // Using a * (a + b) / 2 as the multiplying factor leads to a rapidly
                // converging series in n.  Of course, using this series isn't really
                // necessary, since the exact expression is simple to evaluate.  However,
                // we do it for consistency with RectifyingRadius; and, presumably, the
                // roundoff error is smaller compared to that for the exact expression.
                //
                // Maxima code for these coefficients:
                // c2:subst(e=2*sqrt(n)/(1+n),
                //          (atanh(e)/e * (1-n)^2 + (1+n)^2)/(2*(1+n)))$
                // cf(Lmax):=block([t:expand(ratdisrep(taylor(c2,n,0,Lmax)))],
                //  t:makelist(coeff(t,n,Lmax-j),j,0,Lmax),
                //  map(lambda([x],num(x)/
                //         (if denom(x) = 1 then 1 else real(denom(x)))),t))$
                // N.B. Coeff of n^j = 1                     for j = 0
                //                     -1/3                  for j = 1
                //                     4*(2*j-5)!!/(2*j+1)!! for j > 1

                // GEOGRAPHICLIB_AUXLATITUDE_ORDER == 6
                Span<double> coeff = stackalloc double[] {
                    4/1287.0, 4/693.0, 4/315.0, 4/105.0, 4/15.0,
                    -1/3.0, 1
                  };

                int m = Lmax;
                return _a * (_a + _b) / 2 * PolyVal(m, coeff, _n);
            }
        }

        /// <summary>
        /// <i>a</i>, the equatorial radius of the ellipsoid (meters).
        /// </summary>
        public double EquatorialRadius => _a;

        /// <summary>
        /// <i>b</i>, the polar semi-axis of the ellipsoid (meters).
        /// </summary>
        public double PolarSemiAxis => _b;

        /// <summary>
        /// <i>f</i>, the flattening of the ellipsoid.
        /// </summary>
        public double Flattening => _f;

        private unsafe AuxAngle Parametric(AuxAngle phi, double* diff)
        {
            if (diff != null)
            {
                *diff = _fm1;
            }
            return new AuxAngle(phi.Y * _fm1, phi.X);
        }

        /// <summary>
        /// Convert geographic latitude to parametric latitude.
        /// </summary>
        /// <param name="phi">Geographic latitude.</param>
        /// <param name="diff">The derivative d tan(<i>beta</i>) / d tan(<i>phi</i>).</param>
        /// <returns><i>beta</i>, the parametric latitude.</returns>
        internal unsafe AuxAngle Parametric(AuxAngle phi, out double diff)
        {
            diff = 0;
            fixed (double* pdiff = &diff)
            {
                return Parametric(phi, pdiff);
            }
        }

        /// <summary>
        /// Convert geographic latitude to parametric latitude.
        /// </summary>
        /// <param name="phi">Geographic latitude.</param>
        /// <returns><i>beta</i>, the parametric latitude.</returns>
        internal unsafe AuxAngle Parametric(AuxAngle phi)
        {
            return Parametric(phi, null);
        }

        private unsafe AuxAngle Geocentric(AuxAngle phi, double* diff)
        {
            if (diff != null)
            {
                *diff = _e2m1;
            }
            return new AuxAngle(phi.Y * _e2m1, phi.X);
        }

        /// <summary>
        /// Convert geographic latitude to geocentric latitude.
        /// </summary>
        /// <param name="phi">Geographic latitude.</param>
        /// <param name="diff">The derivative d tan(<i>theta</i>) / d tan(<i>phi</i>).</param>
        /// <returns><i>theta</i>, the geocentric latitude.</returns>
        internal unsafe AuxAngle Geocentric(AuxAngle phi, out double diff)
        {
            diff = 0;
            fixed (double* pdiff = &diff)
            {
                return Geocentric(phi, pdiff);
            }
        }

        /// <summary>
        /// Convert geographic latitude to geocentric latitude.
        /// </summary>
        /// <param name="phi">Geographic latitude.</param>
        /// <returns><i>theta</i>, the geocentric latitude.</returns>
        internal unsafe AuxAngle Geocentric(AuxAngle phi)
        {
            return Geocentric(phi, null);
        }

        private unsafe AuxAngle Rectifying(AuxAngle phi, double* diff)
        {
            AuxAngle beta = Parametric(phi).Normalized();
            double sbeta = Abs(beta.Y), cbeta = Abs(beta.X);
            double a = 1, b = _fm1, ka = _e2, kb = -_e12, ka1 = _e2m1, kb1 = _e12p1,
              smu, cmu, mr;
            if (_f < 0)
            {
                Swap(ref a, ref b); Swap(ref ka, ref kb);
                Swap(ref ka1, ref kb1); Swap(ref sbeta, ref cbeta);
            }
            // now a,b = larger/smaller semiaxis
            // beta now measured from larger semiaxis
            // kb,ka = modulus-squared for distance from beta = 0,pi/2
            // NB kb <= 0; 0 <= ka <= 1
            // sa = b*E(beta,sqrt(kb)), sb = a*E(beta',sqrt(ka))
            //    1 - ka * (1 - sb2) = 1 -ka + ka*sb2
            double
              sb2 = sbeta * sbeta,
              cb2 = cbeta * cbeta,
              db2 = 1 - kb * sb2,
              da2 = ka1 + ka * sb2,
              // DLMF Eq. 19.25.9
              sa = b * sbeta * (EllipticFunction.RF(cb2, db2, 1) -
                                 kb * sb2 * EllipticFunction.RD(cb2, db2, 1) / 3),
              // DLMF Eq. 19.25.10 with complementary angles
              sb = a * cbeta * (ka1 * EllipticFunction.RF(sb2, da2, 1)
                                 + ka * ka1 * cb2 * EllipticFunction.RD(sb2, 1, da2) / 3
                                 + ka * sbeta / Sqrt(da2));
            // sa + sb  = 2*EllipticFunction::RG(a*a, b*b) = a*E(e) = b*E(i*e')
            // mr = a*E(e)*(2/pi) = b*E(i*e')*(2/pi)
            mr = (2 * (sa + sb)) / PI;
            smu = Sin(sa / mr);
            cmu = Sin(sb / mr);
            if (_f < 0) { Swap(ref smu, ref cmu); Swap(ref a, ref b); }
            // mu is normalized
            AuxAngle mu = new AuxAngle(smu, cmu).CopyQuadrant(phi);

            if (diff != null)
            {
                double cphi = phi.Normalized().X, tphi = phi.Tan;
                if (!double.IsInfinity(tphi))
                {
                    cmu = mu.X; cbeta = beta.X;
                    *diff = _fm1 * b / mr * Sq(cbeta / cmu) * (cbeta / cphi);
                }
                else
                    *diff = _fm1 * mr / a;
            }

            return mu;
        }

        /// <summary>
        /// Convert geographic latitude to rectifying latitude.
        /// </summary>
        /// <param name="phi">Geographic latitude.</param>
        /// <param name="diff">The derivative d tan(<i>mu</i>) / d tan(<i>phi</i>).</param>
        /// <returns><i>mu</i>, the rectifying latitude.</returns>
        internal unsafe AuxAngle Rectifying(AuxAngle phi, out double diff)
        {
            diff = 0;
            fixed (double* pdiff = &diff)
            {
                return Rectifying(phi, pdiff);
            }
        }

        /// <summary>
        /// Convert geographic latitude to rectifying latitude.
        /// </summary>
        /// <param name="phi">Geographic latitude.</param>
        /// <returns><i>mu</i>, the rectifying latitude.</returns>
        internal unsafe AuxAngle Rectifying(AuxAngle phi)
        {
            return Rectifying(phi, null);
        }

        /// <summary>
        /// Convert geographic latitude to rectifying latitude.
        /// </summary>
        /// <param name="phi">Geographic conformal.</param>
        /// <param name="diff">The derivative d tan(<i>chi</i>) / d tan(<i>phi</i>).</param>
        /// <returns><i>chi</i>, the conformal latitude.</returns>
        private unsafe AuxAngle Conformal(AuxAngle phi, double* diff)
        {
            double tphi = Abs(phi.Tan), tchi = tphi;
            if (!(!IsFinite(tphi) || tphi == 0 || _f == 0))
            {
                double scphi = sc(tphi),
                  sig = Sinh(_e2 * atanhee(tphi)),
                  scsig = sc(sig);
                if (_f <= 0)
                {
                    tchi = tphi * scsig - sig * scphi;
                }
                else
                {
                    // The general expression for tchi is
                    //   tphi * scsig - sig * scphi
                    // This involves cancellation if f > 0, so change to
                    //   (tphi - sig) * (tphi + sig) / (tphi * scsig + sig * scphi)
                    // To control overflow, write as (sigtphi = sig / tphi)
                    //   (tphi - sig) * (1 + sigtphi) / (scsig + sigtphi * scphi)
                    double sigtphi = sig / tphi, tphimsig;
                    if (sig < tphi / 2)
                        tphimsig = tphi - sig;
                    else
                    {
                        // Still have possibly dangerous cancellation in tphi - sig.
                        //
                        // Write tphi - sig = (1 - e) * Dg(1, e)
                        //   Dg(x, y) = (g(x) - g(y)) / (x - y)
                        //   g(x) = sinh(x * atanh(sphi * x))
                        // Note sinh(atanh(sphi)) = tphi
                        // Turn the crank on divided differences, substitute
                        //   sphi = tphi/sc(tphi)
                        //   atanh(x) = asinh(x/sqrt(1-x^2))
                        double em1 = _e2m1 / (1 + _e),              // 1 - e
                          atanhs = Asinh(tphi),                // atanh(sphi)
                          scbeta = sc(_fm1 * tphi),            // sec(beta)
                          scphibeta = sc(tphi) / scbeta,       // sec(phi)/sec(beta)
                          atanhes = Asinh(_e * tphi / scbeta), // atanh(e * sphi)
                          t1 = (atanhs - _e * atanhes) / 2,
                          t2 = Asinh(em1 * (tphi * scphibeta)) / em1,
                          Dg = Cosh((atanhs + _e * atanhes) / 2) * (Sinh(t1) / t1)
                          * ((atanhs + atanhes) / 2 + (1 + _e) / 2 * t2);
                        tphimsig = em1 * Dg;  // tphi - sig
                    }
                    tchi = tphimsig * (1 + sigtphi) / (scsig + sigtphi * scphi);
                }
            }
            AuxAngle chi = new AuxAngle(tchi).CopyQuadrant(phi);

            if (diff != null)
            {
                if (!double.IsInfinity(tphi))
                {
                    double cchi = chi.Normalized().X,
                      cphi = phi.Normalized().X,
                      cbeta = Parametric(phi, out _).Normalized().X;
                    *diff = _e2m1 * (cbeta / cchi) * (cbeta / cphi);
                }
                else
                {
                    double ss = _f > 0 ? Sinh(_e * Asinh(_e1)) : Sinh(-_e * Atan(_e));
                    *diff = _f > 0 ? 1 / (sc(ss) + ss) : sc(ss) - ss;
                }
            }
            return chi;
        }

        /// <summary>
        /// Convert geographic latitude to rectifying latitude.
        /// </summary>
        /// <param name="phi">Geographic conformal.</param>
        /// <param name="diff">The derivative d tan(<i>chi</i>) / d tan(<i>phi</i>).</param>
        /// <returns><i>chi</i>, the conformal latitude.</returns>
        internal unsafe AuxAngle Conformal(AuxAngle phi, out double diff)
        {
            diff = 0;
            fixed (double* pdiff = &diff)
            {
                return Conformal(phi, pdiff);
            }
        }

        /// <summary>
        /// Convert geographic latitude to rectifying latitude.
        /// </summary>
        /// <param name="phi">Geographic conformal.</param>
        /// <returns><i>chi</i>, the conformal latitude.</returns>
        internal unsafe AuxAngle Conformal(AuxAngle phi)
        {
            return Conformal(phi, null);
        }

        private unsafe AuxAngle Authalic(AuxAngle phi, double* diff)
        {
            double tphi = Abs(phi.Tan);
            AuxAngle xi = phi, phin = phi.Normalized();
            if (!(!IsFinite(tphi) || tphi == 0 || _f == 0))
            {
                double qv = q(tphi),
                  Dqp = Dq(tphi),
                  Dqm = (_q + qv) / (1 + Abs(phin.Y)); // Dq(-tphi)
                xi = new AuxAngle(CopySign(qv, phi.Y), phin.X * Sqrt(Dqp * Dqm));
            }
            if (diff != null)
            {
                if (!double.IsNaN(tphi))
                {
                    double cbeta = Parametric(phi).Normalized().X,
                      cxi = xi.Normalized().X;
                    *diff =
                      (2 / _q) * Sq(cbeta / cxi) * (cbeta / cxi) * (cbeta / phin.X);
                }
                else
                    *diff = _e2m1 * Sqrt(_q / 2);
            }
            return xi;
        }

        /// <summary>
        /// Convert geographic latitude to authalic latitude.
        /// </summary>
        /// <param name="phi">Geographic conformal.</param>
        /// <param name="diff">The derivative d tan(<i>xi</i>) / d tan(<i>phi</i>).</param>
        /// <returns><i>xi</i>, the authalic latitude.</returns>
        internal unsafe AuxAngle Authalic(AuxAngle phi, out double diff)
        {
            diff = 0;
            fixed (double* pdiff = &diff)
            {
                return Authalic(phi, pdiff);
            }
        }

        /// <summary>
        /// Convert geographic latitude to authalic latitude.
        /// </summary>
        /// <param name="phi">Geographic conformal.</param>
        /// <returns><i>xi</i>, the authalic latitude.</returns>
        internal unsafe AuxAngle Authalic(AuxAngle phi)
        {
            return Authalic(phi, null);
        }

        /// <summary>
        /// 1d index into AUXNUMBER x AUXNUMBER data
        /// </summary>
        /// <param name="auxout"></param>
        /// <param name="auxin"></param>
        /// <returns></returns>
        internal static int ind(AuxLatitudeType auxout, AuxLatitudeType auxin)
        {
            return (auxout >= 0 && (int)auxout < AUXNUMBER &&
                    auxin >= 0 && (int)auxin < AUXNUMBER) ?
              AUXNUMBER * (int)auxout + (int)auxin : -1;
        }

        /// <summary>
        /// the function sqrt(1 + tphi^2), convert tan to sec
        /// </summary>
        /// <param name="tphi"></param>
        /// <returns></returns>
        internal static double sc(double tphi) => Hypot(1, tphi);

        /// <summary>
        /// the function tphi / sqrt(1 + tphi^2), convert tan to sin
        /// </summary>
        /// <param name="tphi"></param>
        /// <returns></returns>
        internal static double sn(double tphi) => double.IsInfinity(tphi) ? CopySign(1, tphi) : tphi / sc(tphi);

        /// <summary>
        /// Populate [_c[Lmax * k], _c[Lmax * (k + 1)])
        /// </summary>
        /// <param name="auxin"></param>
        /// <param name="auxout"></param>
        /// <param name="k"></param>
        internal void fillcoeff(AuxLatitudeType auxin, AuxLatitudeType auxout, int k)
        {
            if (k < 0) return;          // auxout or auxin out of range
            if (auxout == auxin)
                _c.AsSpan().Slice(Lmax * k, Lmax).Fill(0);
            else
            {
                int o = s_ptrs[k];
                double d = _n;
                if (auxin <= AuxLatitudeType.Rectifying && auxout <= AuxLatitudeType.Rectifying)
                {
                    for (int l = 0; l < Lmax; ++l)
                    {
                        int m = (Lmax - l - 1) / 2; // order of polynomial in n^2
                        _c[Lmax * k + l] = d * PolyVal(m, s_coeffs.AsSpan().Slice(o), _n2);
                        o += m + 1;
                        d *= _n;
                    }
                }
                else
                {
                    for (int l = 0; l < Lmax; ++l)
                    {
                        int m = (Lmax - l - 1); // order of polynomial in n
                        _c[Lmax * k + l] = d * PolyVal(m, s_coeffs.AsSpan().Slice(o), _n);
                        o += m + 1;
                        d *= _n;
                    }
                }
                // assert (o == ptrs[AUXNUMBER * auxout + auxin + 1])
            }
        }

        /// <summary>
        /// the function atanh(e * sphi)/e; works for e^2 = 0 and e^2 &lt; 0
        /// </summary>
        /// <param name="tphi"></param>
        /// <returns></returns>
        private double atanhee(double tphi)
        {
            double s = _f <= 0 ? sn(tphi) : sn(_fm1 * tphi);
            return _f == 0 ? s :
              // atanh(e * sphi) = asinh(e' * sbeta)
              (_f < 0 ? Atan(_e * s) : Asinh(_e1 * s)) / _e;
        }

        /// <summary>
        /// the function atanh(e * sphi)/e + sphi / (1 - (e * sphi)^2);
        /// </summary>
        /// <param name="tphi"></param>
        /// <returns></returns>
        private double q(double tphi)
        {
            double scbeta = sc(_fm1 * tphi);
            return atanhee(tphi) + (tphi / scbeta) * (sc(tphi) / scbeta);
        }

        /// <summary>
        /// The divided difference of (q(1) - q(sphi)) / (1 - sphi)
        /// </summary>
        /// <param name="tphi"></param>
        /// <returns></returns>
        private double Dq(double tphi)
        {
            double scphi = sc(tphi), sphi = sn(tphi),
                   // d = (1 - sphi) can underflow to zero for large tphi
                   d = tphi > 0 ? 1 / (scphi * scphi * (1 + sphi)) : 1 - sphi;
            if (tphi <= 0)
                // This branch is not reached; this case is open-coded in Authalic.
                return (_q - q(tphi)) / d;
            else if (d == 0)
                return 2 / Sq(_e2m1);
            else
            {
                // General expression for Dq(1, sphi) is
                // atanh(e * d / (1 - e2 * sphi)) / (e * d) +
                //   (1 + e2 * sphi) / ((1 - e2 * sphi * sphi) * e2m1);
                // atanh( e * d / (1 - e2 * sphi))
                // = atanh( e * d * scphi/(scphi - e2 * tphi))
                // =
                double scbeta = sc(_fm1 * tphi);
                return (_f == 0 ? 1 :
                        (_f > 0 ? Asinh(_e1 * d * scphi / scbeta) :
                         Atan(_e * d / (1 - _e2 * sphi))) / (_e * d)) +
                  (_f > 0 ?
                   ((scphi + _e2 * tphi) / (_e2m1 * scbeta)) * (scphi / scbeta) :
                  (1 + _e2 * sphi) / ((1 - _e2 * sphi * sphi) * _e2m1));
            }
        }
    }

    /// <summary>
    /// Represents different types of auxiliary latitudes.
    /// </summary>
    public enum AuxLatitudeType
    {
        /// <summary>
        /// Geographic latitude, <i>phi</i>, φ.
        /// </summary>
        Geographic = 0,

        /// <summary>
        /// Parametric latitude, <i>beta</i>, β.
        /// </summary>
        Parametric = 1,

        /// <summary>
        /// Geocentric latitude, <i>theta</i>, θ.
        /// </summary>
        Geocentric = 2,

        /// <summary>
        /// Rectifying latitude, <i>mu</i>, µ.
        /// </summary>
        Rectifying = 3,

        /// <summary>
        /// Conformal latitude, <i>chi</i>, χ.
        /// </summary>
        Conformal = 4,

        /// <summary>
        /// Authalic latitude, <i>xi</i>, ξ.
        /// </summary>
        Authalic = 5,

        /// <summary>
        /// An alias for <see cref="Geographic"/>.
        /// </summary>
        Phi = Geographic,

        /// <summary>
        /// An alias for <see cref="Parametric"/>.
        /// </summary>
        Beta = Parametric,

        /// <summary>
        /// An alias for <see cref="Geocentric"/>.
        /// </summary>
        Theta = Geocentric,

        /// <summary>
        /// An alias for <see cref="Rectifying"/>.
        /// </summary>
        Mu = Rectifying,

        /// <summary>
        /// An alias for <see cref="Conformal"/>.
        /// </summary>
        Chi = Conformal,

        /// <summary>
        /// An alias for <see cref="Authalic"/>.
        /// </summary>
        Xi = Authalic,

        /// <summary>
        /// An alias for <see cref="Geographic"/>.
        /// </summary>
        Common = Geographic,

        /// <summary>
        /// An alias for <see cref="Geographic"/>.
        /// </summary>
        Geodetic = Geographic,

        /// <summary>
        /// An alias for <see cref="Parametric"/>.
        /// </summary>
        Reduced = Parametric
    }
}
