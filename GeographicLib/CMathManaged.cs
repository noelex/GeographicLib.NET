using static System.Math;

namespace GeographicLib
{
    partial class CMathManaged : CMath
    {
        private const double
                // SPLIT = 4294967297, // (0x1p32 + 1)
                SPLIT = 134217729, // (0x1p27 + 1)
                _0x1p64 = 1.8446744073709552E+19,
                _0x1p700 = 5.2601359015483735E+210,
                _0x1p_700 = 1.90109156629516E-211;

        /// <summary>
        /// Computes the remainder of two integer values, 
        /// and stores an integer value with the sign and approximate magnitude of the quotient in a location that's specified in a parameter.
        /// </summary>
        /// <param name="x">The numerator.</param>
        /// <param name="y">The denominator.</param>
        /// <param name="quo">A pointer to an integer to store a value that has the sign and approximate magnitude of the quotient.</param>
        /// <returns>
        /// Returns the floating-point remainder of <paramref name="x"/> / <paramref name="y"/>.
        /// If the value of <paramref name="y"/> is 0.0, this method returns a quiet <see cref="double.NaN"/>.</returns>
        public override double Remquo(double x, double y, out int quo)
        {
            Bit64 ux = x, uy = y;

            int ex = (int)(ux.UInt64 >> 52 & 0x7ff);
            int ey = (int)(uy.UInt64 >> 52 & 0x7ff);
            int sx = (int)(ux.UInt64 >> 63);
            int sy = (int)(uy.UInt64 >> 63);
            uint q;
            ulong i;
            ulong uxi = ux.UInt64;

            quo = 0;
            if (uy.UInt64 << 1 == 0 || double.IsNaN(y) || ex == 0x7ff)
                return (x * y) / (x * y);
            if (ux.UInt64 << 1 == 0)
                return x;

            /* normalize x and y */
            if (ex == 0)
            {
                for (i = uxi << 12; i >> 63 == 0; ex--, i <<= 1) ;
                uxi <<= -ex + 1;
            }
            else
            {
                unchecked
                {
                    uxi &= (ulong)-1L >> 12;
                    uxi |= 1L << 52;
                }
            }
            if (ey == 0)
            {
                for (i = uy.UInt64 << 12; i >> 63 == 0; ey--, i <<= 1) ;
                uy.UInt64 <<= -ey + 1;
            }
            else
            {
                unchecked
                {
                    uy.UInt64 &= (ulong)-1L >> 12;
                    uy.UInt64 |= 1L << 52;
                }
            }

            q = 0;
            if (ex < ey)
            {
                if (ex + 1 == ey)
                    goto end;
                return x;
            }

            /* x mod y */
            for (; ex > ey; ex--)
            {
                i = uxi - uy.UInt64;
                if (i >> 63 == 0)
                {
                    uxi = i;
                    q++;
                }
                uxi <<= 1;
                q <<= 1;
            }
            i = uxi - uy.UInt64;
            if (i >> 63 == 0)
            {
                uxi = i;
                q++;
            }
            if (uxi == 0)
                ex = -60;
            else
                for (; uxi >> 52 == 0; uxi <<= 1, ex--) ;
                end:
            /* scale result and decide between |x| and |x|-|y| */
            if (ex > 0)
            {
                uxi -= 1L << 52;
                uxi |= (ulong)ex << 52;
            }
            else
            {
                uxi >>= -ex + 1;
            }
            ux.UInt64 = uxi;
            x = ux.Double;
            if (sy != 0)
                y = -y;
            if (ex == ey || (ex + 1 == ey && (2 * x > y || (2 * x == y && (q % 2) != 0))))
            {
                x -= y;
                q++;
            }
            q &= 0x7fffffff;
            quo = (sx ^ sy) != 0 ? -(int)q : (int)q;
            return sx != 0 ? -x : x;
        }

        /// <summary>
        /// Computes the square root of the sum of the squares of x and y.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns>Hypotenuse of a right-angled triangle computed as √(x^2+y^2).</returns>
        public override double Hypot(double x, double y)
        {
            void sq(out double hi, out double lo, double a)
            {
                double xh, xl, xc;

                xc = a * SPLIT;
                xh = a - xc + xc;
                xl = a - xh;
                hi = a * a;
                lo = xh * xh - hi + 2 * xh * xl + xl * xl;
            }

            Bit64 ux = x, uy = y;

            int ex, ey;
            double hx, lx, hy, ly, z;

            unchecked
            {
                /* arrange |x| >= |y| */
                ux.UInt64 &= (ulong)-1L >> 1;
                uy.UInt64 &= (ulong)-1L >> 1;
            }

            if (ux.UInt64 < uy.UInt64)
            {
                MathEx.Swap(ref ux, ref uy);
            }

            /* special cases */
            ex = (int)(ux.UInt64 >> 52);
            ey = (int)(uy.UInt64 >> 52);
            x = ux.Double;
            y = uy.Double;
            /* note: hypot(inf,nan) == inf */
            if (ey == 0x7ff)
                return y;
            if (ex == 0x7ff || uy.UInt64 == 0)
                return x;
            /* note: hypot(x,y) ~= x + y*y/x/2 with inexact for small y/x */
            /* 64 difference is enough for ld80 double_t */
            if (ex - ey > 64)
                return x + y;

            /* precise sqrt argument in nearest rounding mode without overflow */
            /* xh*xh must not overflow and xl*xl must not underflow in sq */
            z = 1;
            if (ex > 0x3ff + 510)
            {
                z = _0x1p700;
                x *= _0x1p_700;
                y *= _0x1p_700;
            }
            else if (ey < 0x3ff - 450)
            {
                z = _0x1p_700;
                x *= _0x1p700;
                y *= _0x1p700;
            }
            sq(out hx, out lx, x);
            sq(out hy, out ly, y);
            return z * Sqrt(ly + lx + hy + hx);
        }

        private const double _0x1p1023 = 8.98846567431158E+307;

        /// <summary>
        /// Compute log(1+x) without losing precision for small values of <paramref name="x"/>.
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public override double Log1p(double x)
        {
            const double
                ln2_hi = 6.93147180369123816490e-01,  /* 3fe62e42 fee00000 */
                ln2_lo = 1.90821492927058770002e-10,  /* 3dea39ef 35793c76 */
                Lg1 = 6.666666666666735130e-01,  /* 3FE55555 55555593 */
                Lg2 = 3.999999999940941908e-01,  /* 3FD99999 9997FA04 */
                Lg3 = 2.857142874366239149e-01,  /* 3FD24924 94229359 */
                Lg4 = 2.222219843214978396e-01,  /* 3FCC71C5 1D8E78AF */
                Lg5 = 1.818357216161805012e-01,  /* 3FC74664 96CB03DE */
                Lg6 = 1.531383769920937332e-01,  /* 3FC39A09 D078C69F */
                Lg7 = 1.479819860511658591e-01;  /* 3FC2F112 DF3E5244 */

            Bit64 u = x;

            double hfsq, f = 0, c = 0, s, z, R, w, t1, t2, dk;
            uint hx, hu;
            int k;

            hx = (uint)(u.UInt64 >> 32);
            k = 1;
            if (hx < 0x3fda827a || (hx >> 31) != 0)
            {  /* 1+x < sqrt(2)+ */
                if (hx >= 0xbff00000)
                {  /* x <= -1.0 */
                    if (x == -1)
                        return x / 0.0; /* log1p(-1) = -inf */
                    return (x - x) / 0.0;     /* log1p(x<-1) = NaN */
                }
                if (hx << 1 < 0x3ca00000 << 1)
                {  /* |x| < 2**-53 */
                    /* underflow if subnormal */
                    if ((hx & 0x7ff00000) == 0)
                        _ = (float)x;// FORCE_EVAL((float)x);
                    return x;
                }
                if (hx <= 0xbfd2bec4)
                {  /* sqrt(2)/2- <= 1+x < sqrt(2)+ */
                    k = 0;
                    c = 0;
                    f = x;
                }
            }
            else if (hx >= 0x7ff00000)
                return x;
            if (k != 0)
            {
                u.Double = 1 + x;
                hu = (uint)(u.UInt64 >> 32);
                hu += 0x3ff00000 - 0x3fe6a09e;
                k = (int)(hu >> 20) - 0x3ff;
                /* correction term ~ log(1+x)-log(u), avoid underflow in c/u */
                if (k < 54)
                {
                    c = k >= 2 ? 1 - (u.Double - x) : x - (u.Double - 1);
                    c /= u.Double;
                }
                else
                    c = 0;
                /* reduce u into [sqrt(2)/2, sqrt(2)] */
                hu = (hu & 0x000fffff) + 0x3fe6a09e;
                u.UInt64 = (ulong)hu << 32 | (u.UInt64 & 0xffffffff);
                f = u.Double - 1;
            }
            hfsq = 0.5 * f * f;
            s = f / (2.0 + f);
            z = s * s;
            w = z * z;
            t1 = w * (Lg2 + w * (Lg4 + w * Lg6));
            t2 = z * (Lg1 + w * (Lg3 + w * (Lg5 + w * Lg7)));
            R = t2 + t1;
            dk = k;
            return s * (hfsq + R) + (dk * ln2_lo + c) - hfsq + f + dk * ln2_hi;
        }

        /// <summary>
        /// Compute exp(x) - 1 without loss of precision for small values of x.
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public override double Expm1(double x)
        {
            const double
                o_threshold = 7.09782712893383973096e+02, /* 0x40862E42, 0xFEFA39EF */
                ln2_hi = 6.93147180369123816490e-01, /* 0x3fe62e42, 0xfee00000 */
                ln2_lo = 1.90821492927058770002e-10, /* 0x3dea39ef, 0x35793c76 */
                invln2 = 1.44269504088896338700e+00, /* 0x3ff71547, 0x652b82fe */
                /* Scaled Q's: Qn_here = 2**n * Qn_above, for R(2*z) where z = hxs = x*x/2: */
                Q1 = -3.33333333333331316428e-02, /* BFA11111 111110F4 */
                Q2 = 1.58730158725481460165e-03, /* 3F5A01A0 19FE5585 */
                Q3 = -7.93650757867487942473e-05, /* BF14CE19 9EAADBB7 */
                Q4 = 4.00821782732936239552e-06, /* 3ED0CFCA 86E65239 */
                Q5 = -2.01099218183624371326e-07; /* BE8AFDB7 6E09C32D */

            double y, hi, lo, c = 0, t, e, hxs, hfx, r1, twopk;
            Bit64 u = x;

            uint hx = (uint)(u.UInt64 >> 32 & 0x7fffffff);
            int k, sign = (int)(u.UInt64 >> 63);

            /* filter out huge and non-finite argument */
            if (hx >= 0x4043687A)
            {  /* if |x|>=56*ln2 */
                if (double.IsNaN(x))
                    return x;
                if (sign != 0)
                    return -1;
                if (x > o_threshold)
                {
                    x *= _0x1p1023;
                    return x;
                }
            }

            /* argument reduction */
            if (hx > 0x3fd62e42)
            {  /* if  |x| > 0.5 ln2 */
                if (hx < 0x3FF0A2B2)
                {  /* and |x| < 1.5 ln2 */
                    if (sign == 0)
                    {
                        hi = x - ln2_hi;
                        lo = ln2_lo;
                        k = 1;
                    }
                    else
                    {
                        hi = x + ln2_hi;
                        lo = -ln2_lo;
                        k = -1;
                    }
                }
                else
                {
                    k = (int)(invln2 * x + (sign != 0 ? -0.5 : 0.5));
                    t = k;
                    hi = x - t * ln2_hi;  /* t*ln2_hi is exact here */
                    lo = t * ln2_lo;
                }
                x = hi - lo;
                c = (hi - x) - lo;
            }
            else if (hx < 0x3c900000)
            {  /* |x| < 2**-54, return x */
                if (hx < 0x00100000)
                    _ = (float)x;//FORCE_EVAL((float)x);
                return x;
            }
            else
                k = 0;

            /* x is now in primary range */
            hfx = 0.5 * x;
            hxs = x * hfx;
            r1 = 1.0 + hxs * (Q1 + hxs * (Q2 + hxs * (Q3 + hxs * (Q4 + hxs * Q5))));
            t = 3.0 - r1 * hfx;
            e = hxs * ((r1 - t) / (6.0 - x * t));
            if (k == 0)   /* c is 0 */
                return x - (x * e - hxs);
            e = x * (e - c) - c;
            e -= hxs;
            /* exp(x) ~ 2^k (x_reduced - e + 1) */
            if (k == -1)
                return 0.5 * (x - e) - 0.5;
            if (k == 1)
            {
                if (x < -0.25)
                    return -2.0 * (e - (x + 0.5));
                return 1.0 + 2.0 * (x - e);
            }
            u.UInt64 = (ulong)(0x3ff + k) << 52;  /* 2^k */
            twopk = u.Double;
            if (k < 0 || k > 56)
            {  /* suffice to return exp(x)-1 */
                y = x - e + 1.0;
                if (k == 1024)
                    y = y * 2.0 * _0x1p1023;
                else
                    y = y * twopk;
                return y - 1.0;
            }
            u.UInt64 = (ulong)(0x3ff - k) << 52;  /* 2^-k */
            if (k < 20)
                y = (x - e + (1 - u.Double)) * twopk;
            else
                y = (x - (e + u.Double) + 1) * twopk;
            return y;
        }

        public override double Frexp(double x, out int e)
        {
            Bit64 y = x;

            e = 0;
            int ee = (int)(y.UInt64 >> 52 & 0x7ff);

            if (ee == 0)
            {
                if (x != 0)
                {
                    x = Frexp(x * _0x1p64, out e);
                    e -= 64;
                }
                else e = 0;
                return x;
            }
            else if (ee == 0x7ff)
            {
                return x;
            }

            e = ee - 0x3fe;
            y.UInt64 &= 0x800ffffffffffffful;
            y.UInt64 |= 0x3fe0000000000000ul;
            return y.Double;
        }

#if !NET5_0_OR_GREATER

        private const int
            LOG2_TABLE_BITS = 6,
            LOG2_POLY_ORDER = 7,
            LOG2_POLY1_ORDER = 11,
            N = 1 << LOG2_TABLE_BITS;

        private const double
            _0x1p_1022 = 2.2250738585072014E-308, // 0x1p-1022
            _0x1p53 = 9007199254740992, // 0x1p53
            _0x1p63 = 9.2233720368547758E+18;

        private static readonly double
            _0x0_ffffff8p_63 = (double)(Bit64)(4323455642007240704L);

        private static readonly log2_data __log2_data = new log2_data
        {
            invln2hi = (double)(Bit64)(4609176140020449280L),
            invln2lo = (double)(Bit64)(4460540536611119616L),
            poly = new[]
            {
                (double)(Bit64)(-4618699496460942535L),
                (double)(Bit64)(4602334714382648510L),
                (double)(Bit64)(-4623203096100589573L),
                (double)(Bit64)(4598869476581264456L),
                (double)(Bit64)(-4625540194963385668L),
                (double)(Bit64)(4596594284514769198L),
            },
            poly1 = new[]
            {
                (double)(Bit64)(-4618699496460942594L),
                (double)(Bit64)(4602334714382648311L),
                (double)(Bit64)(-4623203096088314817L),
                (double)(Bit64)(4598869476596800484L),
                (double)(Bit64)(-4625540922078752795L),
                (double)(Bit64)(4596593529635268406L),
                (double)(Bit64)(-4627706728612059237L),
                (double)(Bit64)(4594943554432801610L),
                (double)(Bit64)(-4628985851833206724L),
                (double)(Bit64)(4593868635039788646L),
            },
            tab = new[]
            {
                new log2_data.tabitem{invc=(double)(Bit64)(4609192499661090040L),logc=(double)(Bit64)(-4620401435072417792L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4609119721358544033L),logc=(double)(Bit64)(-4620547443313422336L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4609048552263439498L),logc=(double)(Bit64)(-4620691827473188864L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4608978935752977541L),logc=(double)(Bit64)(-4620976044666830848L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4608910827810035960L),logc=(double)(Bit64)(-4621258533714272256L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4608844174476640876L),logc=(double)(Bit64)(-4621537994585174016L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4608778932670517176L),logc=(double)(Bit64)(-4621814478722457600L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4608715056657886750L),logc=(double)(Bit64)(-4622088054713425920L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4608652505995016245L),logc=(double)(Bit64)(-4622358774421803008L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4608591236480216541L),logc=(double)(Bit64)(-4622626711648313344L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4608531212716768208L),logc=(double)(Bit64)(-4622891907264966656L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4608472393775796206L),logc=(double)(Bit64)(-4623154431617540096L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4608414747454330217L),logc=(double)(Bit64)(-4623414321238929408L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4608358236140386643L),logc=(double)(Bit64)(-4623671641629847552L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4608302825730857126L),logc=(double)(Bit64)(-4623926447329796096L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4608248488588498998L),logc=(double)(Bit64)(-4624178767387516928L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4608195190521789770L),logc=(double)(Bit64)(-4624428665301803008L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4608142902552467893L),logc=(double)(Bit64)(-4624676184581849088L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4608091596839446450L),logc=(double)(Bit64)(-4624921367194574848L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4608041245183688113L),logc=(double)(Bit64)(-4625164260603383808L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607991821672542679L),logc=(double)(Bit64)(-4625612992611196928L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607943301500489090L),logc=(double)(Bit64)(-4626089859604447232L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607895658995436368L),logc=(double)(Bit64)(-4626562396589875200L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607848870546001381L),logc=(double)(Bit64)(-4627030681651781632L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607802913757468081L),logc=(double)(Bit64)(-4627494786956263424L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607757767841591096L),logc=(double)(Bit64)(-4627954774287458304L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607713409864934365L),logc=(double)(Bit64)(-4628410733127335936L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607669821755832624L),logc=(double)(Bit64)(-4628862708716888064L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607626980202848768L),logc=(double)(Bit64)(-4629310806041968640L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607584869060826940L),logc=(double)(Bit64)(-4629809704699330560L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607543467839055042L),logc=(double)(Bit64)(-4630690701708247040L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607502757543314897L),logc=(double)(Bit64)(-4631564337986371584L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607462725802730023L),logc=(double)(Bit64)(-4632430631954546688L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607423351295589007L),logc=(double)(Bit64)(-4633289804019056640L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607384620866568409L),logc=(double)(Bit64)(-4634141906601590784L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607346516590018701L),logc=(double)(Bit64)(-4635770193422516224L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607309024952925059L),logc=(double)(Bit64)(-4637446934582853632L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607272131607993180L),logc=(double)(Bit64)(-4639512833354366976L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607235821368544930L),logc=(double)(Bit64)(-4642813028705239040L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607200080338058785L),logc=(double)(Bit64)(-4650211837073686528L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4607112595412526955L),logc=(double)(Bit64)(4577625706679238656L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4606976147961868676L),logc=(double)(Bit64)(4584977561914671104L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4606843801959665791L),logc=(double)(Bit64)(4588127867468742656L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4606715378043149032L),logc=(double)(Bit64)(4590199577454788608L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4606590705194581478L),logc=(double)(Bit64)(4591728374053568512L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4606469619618625715L),logc=(double)(Bit64)(4593235018675142656L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4606351967316207939L),logc=(double)(Bit64)(4594195890746384384L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4606237606996243811L),logc=(double)(Bit64)(4594927996701155328L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4606126401917471569L),logc=(double)(Bit64)(4595649931110391808L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4606018221995330636L),logc=(double)(Bit64)(4596361981185310720L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4605912946732701130L),logc=(double)(Bit64)(4597064405355782144L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4605810462065308559L),logc=(double)(Bit64)(4597757451025965056L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4605710654407004217L),logc=(double)(Bit64)(4598308306428592128L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4605613423158833629L),logc=(double)(Bit64)(4598645833741746176L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4605518667591971211L),logc=(double)(Bit64)(4598979039133097984L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4605426298098548397L),logc=(double)(Bit64)(4599308018491162624L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4605336222619703554L),logc=(double)(Bit64)(4599632888650723328L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4605248357336792306L),logc=(double)(Bit64)(4599953748864847872L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4605162622515894769L),logc=(double)(Bit64)(4600270694594048000L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4605078941046797933L),logc=(double)(Bit64)(4600583822264602624L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4604997241283328404L),logc=(double)(Bit64)(4600893218024939520L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4604917450612805666L),logc=(double)(Bit64)(4601198981460250624L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4604839505954084298L),logc=(double)(Bit64)(4601501185019535360L)},
                new log2_data.tabitem{invc=(double)(Bit64)(4604763342788034444L),logc=(double)(Bit64)(4601799915281006592L)},

            },
#if !__FP_FAST_FMA
            tab2 = new[]
            {
                new log2_data.tab2item{chi=(double)(Bit64)(4604402853719116430L),clo=(double)(Bit64)(4362022846676121093L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4604473222222001318L),clo=(double)(Bit64)(4351039760133978025L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4604543589872431813L),clo=(double)(Bit64)(4360186820181135546L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4604613960589203972L),clo=(double)(Bit64)(4362017200139375398L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4604684327854904547L),clo=(double)(Bit64)(4362577641042051761L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4604754697484045840L),clo=(double)(Bit64)(-4867549565258857663L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4604825066178698197L),clo=(double)(Bit64)(4357242425499200374L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4604895435458996646L),clo=(double)(Bit64)(-4872322070862224713L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4604965803051496146L),clo=(double)(Bit64)(-4862759197468145361L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4605036172874022083L),clo=(double)(Bit64)(4354409044815026719L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4605106540704983888L),clo=(double)(Bit64)(-4859942193201572355L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4605176910554010982L),clo=(double)(Bit64)(-4861067402868720827L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4605247277844338501L),clo=(double)(Bit64)(-4866145770030620612L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4605317646212964305L),clo=(double)(Bit64)(-4866421422011946659L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4605388016825550695L),clo=(double)(Bit64)(4355254363668381187L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4605458384160370085L),clo=(double)(Bit64)(4363565660153595394L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4605528752627056047L),clo=(double)(Bit64)(-4870066486530691458L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4605599121485922197L),clo=(double)(Bit64)(4361311456321974891L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4605669489900424491L),clo=(double)(Bit64)(-4863178645741759522L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4605739858968511466L),clo=(double)(Bit64)(4348246094361144892L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4605810228037762381L),clo=(double)(Bit64)(-4859671117886424932L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4605880596196038234L),clo=(double)(Bit64)(4347041365419758135L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4605950965127898740L),clo=(double)(Bit64)(4356097777077636608L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4606021334850119855L),clo=(double)(Bit64)(4356400052040473536L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4606091704796635033L),clo=(double)(Bit64)(4353835862923526601L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4606162073095370670L),clo=(double)(Bit64)(-4865005117326593834L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4606232442403421578L),clo=(double)(Bit64)(-4865564075919651809L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4606302808848468317L),clo=(double)(Bit64)(4362010240734639575L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4606373178132256463L),clo=(double)(Bit64)(-4875098182039641666L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4606443545414737043L),clo=(double)(Bit64)(-4865290399984216589L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4606513914071649036L),clo=(double)(Bit64)(4361457311920661046L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4606584286181359126L),clo=(double)(Bit64)(-4866982593563404276L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4606654653268379537L),clo=(double)(Bit64)(4362158653603537632L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4606725023436613336L),clo=(double)(Bit64)(-4864437948897232902L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4606795391376932876L),clo=(double)(Bit64)(-4860364199508469743L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4606865761279265660L),clo=(double)(Bit64)(-4859524565416184076L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4606936130211511756L),clo=(double)(Bit64)(-4859522701981363313L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4607006497577911562L),clo=(double)(Bit64)(-4881777750598746148L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4607076865290327347L),clo=(double)(Bit64)(-4875115961130309456L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4607147233707525769L),clo=(double)(Bit64)(4361728657611394670L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4607217603241911143L),clo=(double)(Bit64)(4367224947017401985L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4607287971444676200L),clo=(double)(Bit64)(4361448205140792699L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4607358340836992855L),clo=(double)(Bit64)(-4875178626333279343L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4607428709847089100L),clo=(double)(Bit64)(-4858959367304253420L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4607499078026115829L),clo=(double)(Bit64)(-4860515690030777662L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4607569446445093293L),clo=(double)(Bit64)(-4855585107559847953L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4607639815933671137L),clo=(double)(Bit64)(4367365781677602303L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4607710184813904456L),clo=(double)(Bit64)(-4859239979544846944L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4607780553393856110L),clo=(double)(Bit64)(-4861356726666069181L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4607850922517715089L),clo=(double)(Bit64)(4363141597397045634L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4607921291312686810L),clo=(double)(Bit64)(-4868388502340990011L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4607991658823892537L),clo=(double)(Bit64)(4367220851589949177L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4608062027965722754L),clo=(double)(Bit64)(4358905217570067028L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4608132396583914336L),clo=(double)(Bit64)(-4855439679068229141L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4608202766363307724L),clo=(double)(Bit64)(4356653090751095927L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4608273134405792310L),clo=(double)(Bit64)(4367504049820239923L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4608343503126235135L),clo=(double)(Bit64)(4358931040529366312L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4608413872044932273L),clo=(double)(Bit64)(4360056423802313126L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4608484240640089993L),clo=(double)(Bit64)(-4861079047304659758L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4608554609382081523L),clo=(double)(Bit64)(-4866120542784629975L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4608624977183407887L),clo=(double)(Bit64)(4366384804694283484L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4608695346758398831L),clo=(double)(Bit64)(-4861906384145340409L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4608765715311097786L),clo=(double)(Bit64)(4361727028155808514L)},
                new log2_data.tab2item{chi=(double)(Bit64)(4608836083952399489L),clo=(double)(Bit64)(-4858961166217451202L)},
            }
#endif
        };

        private struct log2_data
        {
            public double invln2hi;
            public double invln2lo;
            public double[] poly;
            public double[] poly1;

            public tabitem[] tab;
            
            public struct tabitem
            {
                public double invc, logc;
            }

#if !__FP_FAST_FMA
            public tab2item[] tab2;

            public struct tab2item
            {
                public double chi, clo;
            }
#endif
        }

        private struct num { public ulong m; public int e; public int sign; };

        /// <summary>
        /// Returns x * 2^n computed efficiently.
        /// </summary>
        /// <param name="x">A double-precision floating-point number that specifies the base value.</param>
        /// <param name="n">A 32-bit integer that specifies the power.</param>
        /// <returns>x * 2^n computed efficiently.</returns>
        public override double ScaleB(double x, int n)
        {
            Bit64 u = x;

            double y = x;

            if (n > 1023)
            {
                y *= _0x1p1023;
                n -= 1023;
                if (n > 1023)
                {
                    y *= _0x1p1023;
                    n -= 1023;
                    if (n > 1023)
                        n = 1023;
                }
            }
            else if (n < -1022)
            {
                /* make sure final n < -53 to avoid double
                   rounding in the subnormal range */
                y *= _0x1p_1022 * _0x1p53;
                n += 1022 - 53;
                if (n < -1022)
                {
                    y *= _0x1p_1022 * _0x1p53;
                    n += 1022 - 53;
                    if (n < -1022)
                        n = -1022;
                }
            }
            u.UInt64 = (ulong)(0x3ff + n) << 52;
            x = y * u.Double;
            return x;
        }

        /// <summary>
        /// Returns a value with the magnitude of <paramref name="x"/> and the sign of <paramref name="y"/>.
        /// </summary>
        /// <param name="x">A number whose magnitude is used in the result.</param>
        /// <param name="y">A number whose sign is the used in the result.</param>
        /// <returns>A value with the magnitude of <paramref name="x"/> and the sign of <paramref name="y"/>.</returns>
        public override double CopySign(double x, double y)
        {
            Bit64 ux = x, uy = y;

            unchecked
            {
                ux.UInt64 &= (ulong)-1L / 2;
                ux.UInt64 |= uy.UInt64 & 1UL << 63;
            }
            return ux.Double;
        }

        public override double FusedMultiplyAdd(double x, double y, double z)
        {
            const int ZEROINFNAN = 0x7ff - 0x3ff - 52 - 1;

            num normalize(double x_)
            {
                var ix = (ulong)(Bit64)x_;
                int e_ = (int)(ix >> 52);
                int sign_ = e_ & 0x800;
                e_ &= 0x7ff;
                if (e_ == 0)
                {
                    ix = (ulong)(Bit64)(x_ * _0x1p63);
                    e_ = (int)(ix >> 52 & 0x7ff);
                    e_ = e_ != 0 ? e_ - 63 : 0x800;
                }
                ix &= (1ul << 52) - 1;
                ix |= 1ul << 52;
                ix <<= 1;
                e_ -= 0x3ff + 52 + 1;

                return new num { m = ix, e = e_, sign = sign_ };
            }

            void mul(out ulong hi, out ulong lo, ulong x_, ulong y_)
            {
                ulong t1, t2, t3;
                ulong xlo = (uint)x_, xhi = x_ >> 32;
                ulong ylo = (uint)y_, yhi = y_ >> 32;

                t1 = xlo * ylo;
                t2 = xlo * yhi + xhi * ylo;
                t3 = xhi * yhi;
                lo = t1 + (t2 << 32);
                hi = t3 + (t2 >> 32) + (t1 > lo ? 1UL : 0);
            }

            int a_clz_64(ulong x_)
            {
                x_ |= x_ >> 1;
                x_ |= x_ >> 2;
                x_ |= x_ >> 4;
                x_ |= x_ >> 8;
                x_ |= x_ >> 16;
                x_ |= x_ >> 32;
                x_ -= (x_ >> 1) & 0x5555555555555555;
                x_ = (x_ & 0x3333333333333333) + ((x_ >> 2) & 0x3333333333333333);
                x_ = (x_ + (x_ >> 4)) & 0x0f0f0f0f0f0f0f0f;
                return 64 - (int)((x_ * 0x0101010101010101) >> 56);
            }

            /* normalize so top 10bits and last bit are 0 */
            num nx, ny, nz;
            nx = normalize(x);
            ny = normalize(y);
            nz = normalize(z);

            if (nx.e >= ZEROINFNAN || ny.e >= ZEROINFNAN)
                return x * y + z;
            if (nz.e >= ZEROINFNAN)
            {
                if (nz.e > ZEROINFNAN) /* z==0 */
                    return x * y + z;
                return z;
            }

            /* mul: r = x*y */
            ulong zhi, zlo;
            mul(out var rhi, out var rlo, nx.m, ny.m);
            /* either top 20 or 21 bits of rhi and last 2 bits of rlo are 0 */

            /* align exponents */
            int e = nx.e + ny.e;
            int d = nz.e - e;
            /* shift bits z<<=kz, r>>=kr, so kz+kr == d, set e = e+kr (== ez-kz) */
            if (d > 0)
            {
                if (d < 64)
                {
                    zlo = nz.m << d;
                    zhi = nz.m >> 64 - d;
                }
                else
                {
                    zlo = 0;
                    zhi = nz.m;
                    e = nz.e - 64;
                    d -= 64;
                    if (d == 0)
                    {
                    }
                    else if (d < 64)
                    {
                        rlo = rhi << 64 - d | rlo >> d | ((rlo << 64 - d) == 0 ? 0UL : 1);
                        rhi = rhi >> d;
                    }
                    else
                    {
                        rlo = 1;
                        rhi = 0;
                    }
                }
            }
            else
            {
                zhi = 0;
                d = -d;
                if (d == 0)
                {
                    zlo = nz.m;
                }
                else if (d < 64)
                {
                    zlo = nz.m >> d | ((nz.m << 64 - d) == 0 ? 0UL : 1);
                }
                else
                {
                    zlo = 1;
                }
            }

            /* add */
            int sign = nx.sign ^ ny.sign;
            int samesign = (sign ^ nz.sign) == 0 ? 1 : 0;
            int nonzero = 1;
            if (samesign != 0)
            {
                /* r += z */
                rlo += zlo;
                rhi += zhi + (rlo < zlo ? 1UL : 0);
            }
            else
            {
                /* r -= z */
                ulong t = rlo;
                rlo -= zlo;
                rhi = rhi - zhi - (t < rlo ? 1UL : 0);
                if ((rhi >> 63) != 0)
                {
                    rlo = (ulong)-(long)rlo;
                    rhi = (ulong)-(long)rhi - (rlo == 0 ? 0UL : 1);
                    sign = sign == 0 ? 1 : 0;
                }
                nonzero = rhi == 0 ? 0 : 1;
            }

            /* set rhi to top 63bit of the result (last bit is sticky) */
            if (nonzero != 0)
            {
                e += 64;
                d = a_clz_64(rhi) - 1;
                /* note: d > 0 */
                rhi = rhi << d | rlo >> 64 - d | ((rlo << d) == 0 ? 0UL : 1);
            }
            else if (rlo != 0)
            {
                d = a_clz_64(rlo) - 1;
                if (d < 0)
                    rhi = rlo >> 1 | (rlo & 1);
                else
                    rhi = rlo << d;
            }
            else
            {
                /* exact +-0 */
                return x * y + z;
            }
            e -= d;

            /* convert to double */
            ulong i = rhi; /* i is in [1<<62,(1<<63)-1] */
            if (sign != 0)
                i = (ulong)-(long)i;
            double r = i; /* |r| is in [0x1p62,0x1p63] */

            if (e < -1022 - 62)
            {
                /* result is subnormal before rounding */
                if (e == -1022 - 63)
                {
                    double c = _0x1p63;
                    if (sign != 0)
                        c = -c;
                    if (r == c)
                    {
                        /* min normal after rounding, underflow depends
                           on arch behaviour which can be imitated by
                           a double to float conversion */
                        float fltmin = (float)(_0x0_ffffff8p_63 * MathEx.FLT_MIN * r);
                        return MathEx.DBL_MIN / MathEx.FLT_MIN * fltmin;
                    }
                    /* one bit is lost when scaled, add another top bit to
                       only round once at conversion if it is inexact */
                    if ((rhi << 53) != 0)
                    {
                        i = rhi >> 1 | (rhi & 1) | 1ul << 62;
                        if (sign != 0)
                            i = (ulong)-(long)i;
                        r = i;
                        r = 2 * r - c; /* remove top bit */

                        /* raise underflow portably, such that it
                           cannot be optimized away */
                        {
                            double tiny = MathEx.DBL_MIN / MathEx.FLT_MIN * r;
                            r += (double)(tiny * tiny) * (r - r);
                        }
                    }
                }
                else
                {
                    /* only round once when scaled */
                    d = 10;
                    i = (rhi >> d | ((rlo << 64 - d) == 0 ? 0UL : 1)) << d;
                    if (sign != 0)
                        i = (ulong)-(long)i;
                    r = i;
                }
            }
            return ScaleB(r, e);
        }

        public override double Log2(double x)
        {
            const ulong
                _0x1p52 = 4503599627370496,
                OFF = 0x3fe6000000000000,
                LO = 4606800540372828160, // asuint64(1.0 - 0x1.5b51p-5)
                HI = 4607381812656734208; // asuint64(1.0 + 0x1.6ab2p-5)
            double z, r, r2, r4, y, invc, logc, kd, hi, lo, t1, t2, t3, p;
            ulong ix, iz, tmp;
            uint top;
            int k, i;

            ulong asuint64(Bit64 x_) =>  (ulong)x_;

            double asdouble(Bit64 x_) => (double)x_;

            uint top16(double x_) => (uint)(asuint64(x_) >> 48);

            ix = asuint64(x);
            top = top16(x);

#if !__FP_FAST_FMA
            double rhi, rlo;
            var T2 = __log2_data.tab2;
#endif
            var T = __log2_data.tab;
            var B = __log2_data.poly1;
            var A = __log2_data.poly;
            var InvLn2hi = __log2_data.invln2hi;
            var InvLn2lo = __log2_data.invln2lo;

            if (ix - LO < HI - LO)
            {
                /* Handle close to 1.0 inputs separately.  */
                /* Fix sign of zero with downward rounding when x==1.  */
                if (ix == asuint64(1.0))
                    return 0;
                r = x - 1.0;

#if __FP_FAST_FMA
                hi = r * InvLn2hi;
                lo = r * InvLn2lo + FusedMultiplyAdd(r, InvLn2hi, -hi);

#else
                rhi = asdouble(asuint64(r) & unchecked((ulong)-1L) << 32);
                rlo = r - rhi;
                hi = rhi * InvLn2hi;
                lo = rlo * InvLn2hi + r * InvLn2lo;
#endif

                r2 = r * r; /* rounding error: 0x1p-62.  */
                r4 = r2 * r2;
                /* Worst-case error is less than 0.54 ULP (0.55 ULP without fma).  */
                p = r2 * (B[0] + r * B[1]);
                y = hi + p;
                lo += hi - y + p;
                lo += r4 * (B[2] + r * B[3] + r2 * (B[4] + r * B[5]) +
                        r4 * (B[6] + r * B[7] + r2 * (B[8] + r * B[9])));
                y += lo;
                return y;
            }
            if (top - 0x0010 >= 0x7ff0 - 0x0010)
            {
                /* x < 0x1p-1022 or inf or nan.  */
                if (ix * 2 == 0)
                    return double.NegativeInfinity;
                if (ix == asuint64(double.PositiveInfinity)) /* log(inf) == inf.  */
                    return x;
                if (((top & 0x8000U) != 0) || (top & 0x7ff0U) == 0x7ff0U)
                    return double.NaN;
                /* x is subnormal, normalize it.  */
                ix = asuint64(x * _0x1p52);
                ix -= 52UL << 52;
            }

            /* x = 2^k z; where z is in range [OFF,2*OFF) and exact.
               The range is split into N subintervals.
               The ith subinterval contains z and c is near its center.  */
            tmp = ix - OFF;
            i = (int)(tmp >> (52 - LOG2_TABLE_BITS)) % N;
            k = (int)((long)tmp >> 52); /* arithmetic shift */
            iz = ix - (tmp & 0xfffUL << 52);
            invc = T[i].invc;
            logc = T[i].logc;
            z = asdouble(iz);
            kd = k;

            /* log2(x) = log2(z/c) + log2(c) + k.  */
            /* r ~= z/c - 1, |r| < 1/(2*N).  */

#if __FP_FAST_FMA
            /* rounding error: 0x1p-55/N.  */
            r = FusedMultiplyAdd(z, invc, -1.0);
            t1 = r * InvLn2hi;
            t2 = r * InvLn2lo + FusedMultiplyAdd(r, InvLn2hi, -t1);

#else
            /* rounding error: 0x1p-55/N + 0x1p-65.  */
            r = (z - T2[i].chi - T2[i].clo) * invc;
            rhi = asdouble(asuint64(r) & unchecked((ulong)-1L) << 32);
            rlo = r - rhi;
            t1 = rhi * InvLn2hi;
            t2 = rlo * InvLn2hi + r * InvLn2lo;
#endif

            /* hi + lo = r/ln2 + log2(c) + k.  */
            t3 = kd + logc;
            hi = t3 + t1;
            lo = t3 - hi + t1 + t2;

            /* log2(r+1) = r/ln2 + r^2*poly(r).  */
            /* Evaluation is optimized assuming superscalar pipelined execution.  */
            r2 = r * r; /* rounding error: 0x1p-54/N^2.  */
            r4 = r2 * r2;
            /* Worst-case error if |y| > 0x1p-4: 0.547 ULP (0.550 ULP without fma).
               ~ 0.5 + 2/N/ln2 + abs-poly-error*0x1p56 ULP (+ 0.003 ULP without fma).  */
            p = A[0] + r * A[1] + r2 * (A[2] + r * A[3]) + r4 * (A[4] + r * A[5]);
            y = lo + r2 * p + hi;
            return y;
        }

#endif

#if NETSTANDARD2_0
        /// <summary>
        /// Returns the angle whose hyperbolic tangent is the specified number.
        /// </summary>
        /// <param name="x">A number representing a hyperbolic tangent, where d must be greater than or equal to -1, but less than or equal to 1.</param>
        /// <returns>
        /// An angle, θ, measured in radians, such that -∞ &lt; θ &lt; -1, or 1 &lt; θ &lt; ∞.
        /// -or- <see cref="double.NaN"/> if <paramref name="x"/> &lt; -1 or <paramref name="x"/> > 1 or <paramref name="x"/> equals <see cref="double.NaN"/>.
        /// </returns>
        public override double Atanh(double x)
        {
            Bit64 u = x;

            uint e = (uint)(u.UInt64 >> 52 & 0x7ff);
            uint s = (uint)(u.UInt64 >> 63);
            double y;

            /* |x| */
            u.UInt64 &= unchecked((ulong)-1L / 2);
            y = u.Double;

            if (e < 0x3ff - 1)
            {
                if (e < 0x3ff - 32)
                {
                    /* handle underflow */
                    if (e == 0)
                        _ = (float)y; // FORCE_EVAL((float)y)
                }
                else
                {
                    /* |x| < 0.5, up to 1.7ulp error */
                    y = 0.5 * Log1p(2 * y + 2 * y * y / (1 - y));
                }
            }
            else
            {
                /* avoid overflow */
                y = 0.5 * Log1p(2 * (y / (1 - y)));
            }
            return s != 0 ? -y : y;
        }

        /// <summary>
        /// Returns the angle whose hyperbolic sine is the specified number.
        /// </summary>
        /// <param name="x">A number representing a hyperbolic sine, where d must be greater than or equal
        /// to <see cref="double.NegativeInfinity"/>, but less than or equal to <see cref="double.PositiveInfinity"/>.
        /// </param>
        /// <returns>
        /// An angle, θ, measured in radians, such that -∞ &lt; θ ≤ -1, or 1 ≤ θ &lt; ∞. -or- <see cref="double.NaN"/>
        /// if <paramref name="x"/> equals <see cref="double.NaN"/>.
        /// </returns>
        public override double Asinh(double x)
        {
            Bit64 u = x;

            uint e = (uint)(u.UInt64 >> 52 & 0x7ff);
            uint s = (uint)(u.UInt64 >> 63);

            /* |x| */
            u.UInt64 &= unchecked((ulong)-1L / 2);
            x = u.Double;

            if (e >= 0x3ff + 26)
            {
                /* |x| >= 0x1p26 or inf or nan */
                x = Log(x) + 0.693147180559945309417232121458176568;
            }
            else if (e >= 0x3ff + 1)
            {
                /* |x| >= 2 */
                x = Log(2 * x + 1 / (Sqrt(x * x + 1) + x));
            }
            else if (e >= 0x3ff - 26)
            {
                /* |x| >= 0x1p-26, up to 1.6ulp error in [0.125,0.5] */
                x = Log1p(x + x * x / (Sqrt(x * x + 1) + 1));
            }
            else
            {
                /* |x| < 0x1p-26, raise inexact if x != 0 */
                _ = (float)x + 1.3292279957849159E+36f;//FORCE_EVAL(x + 0x1p120f);
            }
            return s != 0 ? -x : x;
        }

        /// <summary>
        /// Returns the cube root of a specified number.
        /// </summary>
        /// <param name="x">The number whose cube root is to be found.</param>
        /// <returns>The cube root of <paramref name="x"/>. -or- <see cref="double.NaN"/> if <paramref name="x"/> equals <see cref="double.NaN"/>.</returns>
        public override double Cbrt(double x)
        {
            const uint
            B1 = 715094163, /* B1 = (1023-1023/3-0.03306235651)*2**20 */
            B2 = 696219795; /* B2 = (1023-1023/3-54/3-0.03306235651)*2**20 */

            /* |1/cbrt(x) - p(x)| < 2**-23.5 (~[-7.93e-8, 7.929e-8]). */
            const double
            P0 = 1.87595182427177009643,  /* 0x3ffe03e6, 0x0f61e692 */
            P1 = -1.88497979543377169875,  /* 0xbffe28e0, 0x92f02420 */
            P2 = 1.621429720105354466140, /* 0x3ff9f160, 0x4a49d6c2 */
            P3 = -0.758397934778766047437, /* 0xbfe844cb, 0xbee751d9 */
            P4 = 0.145996192886612446982; /* 0x3fc2b000, 0xd4e4edd7 */

            const double _0x1p54 = 18014398509481984d;

            Bit64 u = x;
            double r, s, t, w;
            uint hx = (uint)(u.UInt64 >> 32 & 0x7fffffff);

            if (hx >= 0x7ff00000)  /* cbrt(NaN,INF) is itself */
                return x + x;

            /*
             * Rough cbrt to 5 bits:
             *    cbrt(2**e*(1+m) ~= 2**(e/3)*(1+(e%3+m)/3)
             * where e is integral and >= 0, m is real and in [0, 1), and "/" and
             * "%" are integer division and modulus with rounding towards minus
             * infinity.  The RHS is always >= the LHS and has a maximum relative
             * error of about 1 in 16.  Adding a bias of -0.03306235651 to the
             * (e%3+m)/3 term reduces the error to about 1 in 32. With the IEEE
             * floating point representation, for finite positive normal values,
             * ordinary integer divison of the value in bits magically gives
             * almost exactly the RHS of the above provided we first subtract the
             * exponent bias (1023 for doubles) and later add it back.  We do the
             * subtraction virtually to keep e >= 0 so that ordinary integer
             * division rounds towards minus infinity; this is also efficient.
             */
            if (hx < 0x00100000)
            { /* zero or subnormal? */
                u.Double = x * _0x1p54;
                hx = (uint)(u.UInt64 >> 32 & 0x7fffffff);
                if (hx == 0)
                    return x;  /* cbrt(0) is itself */
                hx = hx / 3 + B2;
            }
            else
                hx = hx / 3 + B1;
            u.UInt64 &= 1UL << 63;
            u.UInt64 |= (ulong)hx << 32;
            t = u.Double;

            /*
             * New cbrt to 23 bits:
             *    cbrt(x) = t*cbrt(x/t**3) ~= t*P(t**3/x)
             * where P(r) is a polynomial of degree 4 that approximates 1/cbrt(r)
             * to within 2**-23.5 when |r - 1| < 1/10.  The rough approximation
             * has produced t such than |t/cbrt(x) - 1| ~< 1/32, and cubing this
             * gives us bounds for r = t**3/x.
             *
             * Try to optimize for parallel evaluation as in __tanf.c.
             */
            r = (t * t) * (t / x);
            t = t * ((P0 + r * (P1 + r * P2)) + ((r * r) * r) * (P3 + r * P4));

            /*
             * Round t away from zero to 23 bits (sloppily except for ensuring that
             * the result is larger in magnitude than cbrt(x) but not much more than
             * 2 23-bit ulps larger).  With rounding towards zero, the error bound
             * would be ~5/6 instead of ~4/6.  With a maximum error of 2 23-bit ulps
             * in the rounded t, the infinite-precision error in the Newton
             * approximation barely affects third digit in the final error
             * 0.667; the error in the rounded t can be up to about 3 23-bit ulps
             * before the final error is larger than 0.667 ulps.
             */
            u.Double = t;
            u.UInt64 = (u.UInt64 + 0x80000000) & 0xffffffffc0000000UL;
            t = u.Double;

            /* one step Newton iteration to 53 bits with error < 0.667 ulps */
            s = t * t;         /* t*t is exact */
            r = x / s;         /* error <= 0.5 ulps; |r| < |t| */
            w = t + t;         /* t+t is exact */
            r = (r - t) / (w + r); /* r-t is exact; w+r ~= 3*t */
            t = t + t * r;       /* error <= 0.5 + 0.5/3 + epsilon */
            return t;
        }
#endif
    }
}
