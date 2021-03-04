#if CMATH_MANAGED
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

using static System.Math;

namespace GeographicLib
{
    static partial class MathEx
    {
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
        public static double Remquo(double x, double y, out int quo)
        {
            Span<double>
                u_f = stackalloc[] { x, y };

            var u_i = MemoryMarshal.Cast<double, ulong>(u_f);

            int ex = (int)(u_i[0] >> 52 & 0x7ff);
            int ey = (int)(u_i[1] >> 52 & 0x7ff);
            int sx = (int)(u_i[0] >> 63);
            int sy = (int)(u_i[1] >> 63);
            uint q;
            ulong i;
            ulong uxi = u_i[0];

            quo = 0;
            if (u_i[1] << 1 == 0 || double.IsNaN(y) || ex == 0x7ff)
                return (x * y) / (x * y);
            if (u_i[0] << 1 == 0)
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
                for (i = u_i[1] << 12; i >> 63 == 0; ey--, i <<= 1) ;
                u_i[1] <<= -ey + 1;
            }
            else
            {
                unchecked
                {
                    u_i[1] &= (ulong)-1L >> 12;
                    u_i[1] |= 1L << 52;
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
                i = uxi - u_i[1];
                if (i >> 63 == 0)
                {
                    uxi = i;
                    q++;
                }
                uxi <<= 1;
                q <<= 1;
            }
            i = uxi - u_i[1];
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
            u_i[0] = uxi;
            x = u_f[0];
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
        public static double Hypot(double x, double y)
        {
            const double
                // SPLIT = 4294967297, // (0x1p32 + 1)
                SPLIT = 134217729, // (0x1p27 + 1)
                _0x1p700 = 5.2601359015483735E+210,
                _0x1p_700 = 1.90109156629516E-211;

            void sq(out double hi, out double lo, double a)
            {
                double xh, xl, xc;

                xc = a * SPLIT;
                xh = a - xc + xc;
                xl = a - xh;
                hi = a * a;
                lo = xh * xh - hi + 2 * xh * xl + xl * xl;
            }

            Span<double>
                u_f = stackalloc[] { x, y };

            var u_i = MemoryMarshal.Cast<double, ulong>(u_f);

            int ex, ey;
            double hx, lx, hy, ly, z, ut;

            unchecked
            {
                /* arrange |x| >= |y| */
                u_i[0] &= (ulong)-1L >> 1;
                u_i[1] &= (ulong)-1L >> 1;
            }

            if (u_i[0] < u_i[1])
            {
                ut = u_f[0];
                u_f[0] = u_f[1];
                u_f[1] = ut;
            }

            /* special cases */
            ex = (int)(u_i[0] >> 52);
            ey = (int)(u_i[1] >> 52);
            x = u_f[0];
            y = u_f[1];
            /* note: hypot(inf,nan) == inf */
            if (ey == 0x7ff)
                return y;
            if (ex == 0x7ff || u_i[1] == 0)
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

        /// <summary>
        /// Compute log(1+x) without losing precision for small values of <paramref name="x"/>.
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double Log1p(double x)
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

            Span<double> u_f = stackalloc[] { x };
            var u_i = MemoryMarshal.Cast<double, ulong>(u_f);

            double hfsq, f = 0, c = 0, s, z, R, w, t1, t2, dk;
            uint hx, hu;
            int k;

            hx = (uint)(u_i[0] >> 32);
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
                u_f[0] = 1 + x;
                hu = (uint)(u_i[0] >> 32);
                hu += 0x3ff00000 - 0x3fe6a09e;
                k = (int)(hu >> 20) - 0x3ff;
                /* correction term ~ log(1+x)-log(u), avoid underflow in c/u */
                if (k < 54)
                {
                    c = k >= 2 ? 1 - (u_f[0] - x) : x - (u_f[0] - 1);
                    c /= u_f[0];
                }
                else
                    c = 0;
                /* reduce u into [sqrt(2)/2, sqrt(2)] */
                hu = (hu & 0x000fffff) + 0x3fe6a09e;
                u_i[0] = (ulong)hu << 32 | (u_i[0] & 0xffffffff);
                f = u_f[0] - 1;
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
        public static double Expm1(double x)
        {
            const double
                _0x1p1023 = 8.98846567431158E+307,
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
            Span<double> u_f = stackalloc[] { x };
            var u_i = MemoryMarshal.Cast<double, ulong>(u_f);

            uint hx = (uint)(u_i[0] >> 32 & 0x7fffffff);
            int k, sign = (int)(u_i[0] >> 63);

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
            u_i[0] = (ulong)(0x3ff + k) << 52;  /* 2^k */
            twopk = u_f[0];
            if (k < 0 || k > 56)
            {  /* suffice to return exp(x)-1 */
                y = x - e + 1.0;
                if (k == 1024)
                    y = y * 2.0 * _0x1p1023;
                else
                    y = y * twopk;
                return y - 1.0;
            }
            u_i[0] = (ulong)(0x3ff - k) << 52;  /* 2^-k */
            if (k < 20)
                y = (x - e + (1 - u_f[0])) * twopk;
            else
                y = (x - (e + u_f[0]) + 1) * twopk;
            return y;
        }

#if !NET5_0
        /// <summary>
        /// Returns x * 2^n computed efficiently.
        /// </summary>
        /// <param name="x">A double-precision floating-point number that specifies the base value.</param>
        /// <param name="n">A 32-bit integer that specifies the power.</param>
        /// <returns>x * 2^n computed efficiently.</returns>
        public static double ScaleB(double x, int n)
        {
            const double
                _0x1p1023 = 8.98846567431158E+307, // 0x1p1023
                _0x1p_1022 = 2.2250738585072014E-308, // 0x1p-1022
                _0x1p53 = 9007199254740992; // 0x1p53

            Span<double>
                u_f = stackalloc[] { x };

            var u_i = MemoryMarshal.Cast<double, ulong>(u_f);

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
            u_i[0] = (ulong)(0x3ff + n) << 52;
            x = y * u_f[0];
            return x;
        }

        /// <summary>
        /// Returns a value with the magnitude of <paramref name="x"/> and the sign of <paramref name="y"/>.
        /// </summary>
        /// <param name="x">A number whose magnitude is used in the result.</param>
        /// <param name="y">A number whose sign is the used in the result.</param>
        /// <returns>A value with the magnitude of <paramref name="x"/> and the sign of <paramref name="y"/>.</returns>
        public static double CopySign(double x, double y)
        {
            Span<double>
                u_f = stackalloc[] { x, y };

            var u_i = MemoryMarshal.Cast<double, ulong>(u_f);

            unchecked
            {
                u_i[0] &= (ulong)-1L / 2;
                u_i[0] |= u_i[1] & 1UL << 63;
            }
            return u_f[0];
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
        public static double Atanh(double x)
        {
            Span<double> u_f = stackalloc[] { x };
            var u_i = MemoryMarshal.Cast<double, ulong>(u_f);

            uint e = (uint)(u_i[0] >> 52 & 0x7ff);
            uint s = (uint)(u_i[0] >> 63);
            double y;

            /* |x| */
            u_i[0] &= unchecked((ulong)-1L / 2);
            y = u_f[0];

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
        public static double Asinh(double x)
        {
            Span<double> u_f = stackalloc[] { x };
            var u_i = MemoryMarshal.Cast<double, ulong>(u_f);

            uint e = (uint)(u_i[0] >> 52 & 0x7ff);
            uint s = (uint)(u_i[0] >> 63);

            /* |x| */
            u_i[0] &= unchecked((ulong)-1L / 2);
            x = u_f[0];

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
        public static double Cbrt(double x)
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

            Span<double> u_f = stackalloc[] { x };
            var u_i = MemoryMarshal.Cast<double, ulong>(u_f);
            double r, s, t, w;
            uint hx = (uint)(u_i[0] >> 32 & 0x7fffffff);

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
                u_f[0] = x * _0x1p54;
                hx = (uint)(u_i[0] >> 32 & 0x7fffffff);
                if (hx == 0)
                    return x;  /* cbrt(0) is itself */
                hx = hx / 3 + B2;
            }
            else
                hx = hx / 3 + B1;
            u_i[0] &= 1UL << 63;
            u_i[0] |= (ulong)hx << 32;
            t = u_f[0];

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
            u_f[0] = t;
            u_i[0] = (u_i[0] + 0x80000000) & 0xffffffffc0000000UL;
            t = u_f[0];

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
#endif
