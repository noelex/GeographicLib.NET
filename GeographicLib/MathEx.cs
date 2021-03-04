using System;
using System.Collections.Generic;
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;

using static System.Math;
using static GeographicLib.Macros;

namespace GeographicLib
{
    /// <summary>
    /// Defines extended mathematical methods.
    /// </summary>
    public static partial class MathEx
    {
        /// <summary>
        /// The number of radians in a degree
        /// </summary>
        public const double Degree = PI / 180;

        /// <summary>
        /// 
        /// </summary>
        internal const int FLT_RADIX = 2;

        /// <summary>
        /// Minimum positive, normal value of a <see cref="Double"/> (<c>2.2250738585072014e-308</c>).
        /// </summary>
        internal const double DBL_MIN = 2.2250738585072014e-308;

        /// <summary>
        /// Gives an upper bound on the relative error due to rounding in floating point arithmetic. 
        /// Not to be confused with <see cref="double.Epsilon"/>.
        /// </summary>
        internal const double DBL_EPSILON = 2.2204460492503131e-16;

        /// <summary>
        /// Maximum binary exponent of a <see cref="double"/>.
        /// </summary>
        internal const int DBL_MAX_EXP = 1024;

        /// <summary>
        /// Number of digits in the radix specified by FLT_RADIX in the floating-point significand.
        /// </summary>
        internal const int DBL_MANT_DIG = 53;

#if NETSTANDARD2_0
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static bool IsFinite(double x) => unchecked((ulong)(BitConverter.DoubleToInt64Bits(x) & -1L >> 1)) < 0x7ffUL << 52;
#else
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static bool IsFinite(double x) => double.IsFinite(x);
#endif

        /// <summary>
        /// Multiplies a floating-point number by an integral power of two.
        /// </summary>
        /// <param name="number">Floating-point value.</param>
        /// <param name="exp">Integer exponent.</param>
        /// <returns>
        /// Return the value of x * 2^<paramref name="exp"/> if successful.
        /// On overflow, and depending on the sign of <paramref name="number"/>,
        /// <see cref="double.NegativeInfinity"/> or <see cref="double.PositiveInfinity"/> is returned.
        /// </returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Ldexp(double number, int exp) => ScaleB(number, exp);

        /// <summary>
        /// Swap value of two variables.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="x"></param>
        /// <param name="y"></param>
        public static void Swap<T>(ref T x, ref T y) where T : struct
        {
            var a = x;
            x = y;
            y = a;
        }

        /// <summary>
        /// Normalize a two-vector.
        /// </summary>
        /// <param name="x">on output set to <i>x</i>/Hypot(<i>x</i>, <i>y</i>).</param>
        /// <param name="y">on output set to <i>y</i>/hypot(<i>x</i>, <i>y</i>).</param>
        public static void Norm(ref double x, ref double y)
        {
            var h = Hypot(x, y);
            x /= h;
            y /= h;
        }

        /// <summary>
        /// Computes the square root of the sum of the squares of 1 and x.
        /// </summary>
        /// <param name="x"></param>
        /// <returns>Hypotenuse of a right-angled triangle computed as √(1^2+x^2).</returns>
        public static double Hypot(double x) => Hypot(1, x);

        /// <summary>
        /// Square a number.
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double Sq(double x) => x * x;

        /// <summary>
        /// Evaluate the sine and cosine function with the argument in degrees.
        /// </summary>
        /// <param name="x">x in degrees.</param>
        /// <param name="sinx">sin(<i>x</i>).</param>
        /// <param name="cosx">cos(<i>x</i>).</param>
        /// <remarks>
        /// The results obey exactly the elementary properties of the trigonometric
        /// functions, e.g., sin 9° = cos 81° = - sin 123456789°.
        /// If x = -0, then sinx = -0; this is the only case where
        /// -0 is returned.
        /// </remarks>
        public static void SinCosd(double x, out double sinx, out double cosx)
        {
            // In order to minimize round-off errors, this function exactly reduces
            // the argument to the range [-45, 45] before converting it to radians.
            double r;
            int q;

            r = Remquo(x, 90, out q);
            r *= Degree;

            double s = Sin(r), c = Cos(r);

            switch ((uint)q & 3U)
            {
                case 0U: sinx = s; cosx = c; break;
                case 1U: sinx = c; cosx = -s; break;
                case 2U: sinx = -s; cosx = -c; break;
                default: sinx = -c; cosx = s; break; // case 3U
            }

            // TODO: Handle signed zero.
            // Set sign of 0 results.  -0 only produced for sin(-0)
            if (x != 0) { sinx += 0d; cosx += 0d; }
        }

        /// <summary>
        /// Evaluate the sine function with the argument in degrees.
        /// </summary>
        /// <param name="x">in degrees.</param>
        /// <returns>sin(<i>x</i>)</returns>
        public static double Sind(double x)
        {
            // See sincosd
            double r; int q;
            r = Remquo(x, 90d, out q); // now abs(r) <= 45
            r *= Degree;

            var p = (uint)q;

            r = (p & 1U) != 0 ? Cos(r) : Sin(r);
            if ((p & 2U) != 0) r = -r;

            if (x != 0) r += 0d;

            return r;
        }

        /// <summary>
        /// Evaluate the cosine function with the argument in degrees.
        /// </summary>
        /// <param name="x">in degrees.</param>
        /// <returns>cos(<i>x</i>)</returns>
        public static double Cosd(double x)
        {
            // See sincosd
            double r; int q;
            r = Remquo(x, 90d, out q); // now abs(r) <= 45
            r *= Degree;

            var p = (uint)(q + 1);
            r = (p & 1U) != 0 ? Cos(r) : Sin(r);
            if ((p & 2U) != 0) r = -r;
            return 0d + r;
        }

        /// <summary>
        /// Evaluate the tangent function with the argument in degrees.
        /// </summary>
        /// <param name="x">in degrees.</param>
        /// <returns>tan(<i>x</i>).</returns>
        /// <remarks>
        /// If <paramref name="x"/> = ±90°, then a suitably large (but finite) value is returned.
        /// </remarks>
        public static double Tand(double x)
        {
            const double overflow = 1 / (DBL_EPSILON * DBL_EPSILON);

            double s, c;
            SinCosd(x, out s, out c);
            return c != 0 ? s / c : (s < 0 ? -overflow : overflow);
        }

        /// <summary>
        /// Evaluate the atan2 function with the result in degrees.
        /// </summary>
        /// <param name="y"></param>
        /// <param name="x"></param>
        /// <returns>atan2(<i>y</i>, <i>x</i>) in degrees.</returns>
        /// <remarks>
        /// The result is in the range (-180°, 180°].  N.B.,
        /// atan2d(±0, -1) = +180°; atan2d(-ε,-1) = -180°, for ε positive and tiny;
        /// atan2d(±0, +1) = ±0°.
        /// </remarks>
        public static double Atan2d(double y, double x)
        {
            // In order to minimize round-off errors, this function rearranges the
            // arguments so that result of atan2 is in the range [-pi/4, pi/4] before
            // converting it to degrees and mapping the result to the correct
            // quadrant.
            int q = 0;
            if (Abs(y) > Abs(x)) { Swap(ref x, ref y); q = 2; }
            if (x < 0) { x = -x; ++q; }
            // here x >= 0 and x >= abs(y), so angle is in [-pi/4, pi/4]
            var ang = Atan2(y, x) / Degree;

            // Note that atan2d(-0.0, 1.0) will return -0.  However, we expect that
            // atan2d will not be called with y = -0.  If need be, include
            //
            //   case 0: ang = 0 + ang; break;
            //
            // and handle mpfr as in AngRound.
            switch (q)
            {
                case 0: return 0 + ang;
                case 1: return (y >= 0 ? 180 : -180) - ang;
                case 2: return 90 - ang;
                case 3: return -90 + ang;
                default: throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Evaluate <i>e</i> atanh(<i>e x</i>)
        /// </summary>
        /// <param name="x"></param>
        /// <param name="es">The signed eccentricity = sign(<i>e^2</i>)sqrt(|<i>e^2</i>|)</param>
        /// <returns><i>e</i> atanh(<i>e x</i>)</returns>
        /// <remarks>
        /// If <i>e^2</i> is negative (<i>e</i> is imaginary), the expression is evaluated in terms of atan.
        /// </remarks>
        public static double EAtanhE(double x, double es) => es > 0d ? es * Atanh(es * x) : -es * Atan(es * x);

        /// <summary>
        /// tanχ in terms of tanφ
        /// </summary>
        /// <param name="tau">τ = tanφ</param>
        /// <param name="es">The signed eccentricity = sign(<i>e^2</i>)sqrt(|<i>e^2</i>|)</param>
        /// <returns>τ′ = tanχ</returns>
        /// <remarks>
        /// See Eqs. (7--9) of
        /// C.F.F.Karney,
        /// <a href = "https://doi.org/10.1007/s00190-011-0445-3" >
        /// Transverse Mercator with an accuracy of a few nanometers,</a>
        /// J.Geodesy 85(8), 475--485 (Aug. 2011)
        /// (preprint
        /// <a href = "https://arxiv.org/abs/1002.1417" > arXiv:1002.1417</a>).
        /// </remarks>
        public static double Taupf(double tau, double es)
        {
            if (IsFinite(tau))
            {
                double tau1 = Hypot(1d, tau),
                  sig = Sinh(EAtanhE(tau / tau1, es));
                return Hypot(1d, sig) * tau - sig * tau1;
            }
            else
                return tau;
        }

        /// <summary>
        /// tanφ in terms of tanχ
        /// </summary>
        /// <param name="taup">τ′ = tanχ</param>
        /// <param name="es">The signed eccentricity = sign(<i>e^2</i>)sqrt(|<i>e^2</i>|)</param>
        /// <returns>τ = tanφ</returns>
        /// <remarks>
        /// See Eqs. (19--21) of
        /// C.F.F.Karney,
        /// <a href = "https://doi.org/10.1007/s00190-011-0445-3" >
        /// Transverse Mercator with an accuracy of a few nanometers,</a>
        /// J.Geodesy 85(8), 475--485 (Aug. 2011)
        /// (preprint
        /// <a href = "https://arxiv.org/abs/1002.1417" > arXiv:1002.1417</a>).
        /// </remarks>
        public static double Tauf(double taup, double es)
        {
            const int numit = 5;

            // min iterations = 1, max iterations = 2; mean = 1.95
            var tol = Sqrt(DBL_EPSILON) / 10;
            var taumax = 2 / Sqrt(DBL_EPSILON);

            double e2m = 1d - Sq(es),
              // To lowest order in e^2, taup = (1 - e^2) * tau = _e2m * tau; so use
              // tau = taup/e2m as a starting guess. Only 1 iteration is needed for
              // |lat| < 3.35 deg, otherwise 2 iterations are needed.  If, instead, tau
              // = taup is used the mean number of iterations increases to 1.999 (2
              // iterations are needed except near tau = 0).
              //
              // For large tau, taup = exp(-es*atanh(es)) * tau.  Use this as for the
              // initial guess for |taup| > 70 (approx |phi| > 89deg).  Then for
              // sufficiently large tau (such that sqrt(1+tau^2) = |tau|), we can exit
              // with the intial guess and avoid overflow problems.  This also reduces
              // the mean number of iterations slightly from 1.963 to 1.954.
              tau = Abs(taup) > 70 ? taup * Exp(EAtanhE(1d, es)) : taup / e2m,
              stol = tol * Max(1d, Abs(taup));

            if (!(Abs(tau) < taumax)) return tau; // handles +/-inf and nan

            for (int i = 0; i < numit || GEOGRAPHICLIB_PANIC; ++i)
            {
                double taupa = Taupf(tau, es),
                  dtau = (taup - taupa) * (1 + e2m * Sq(tau)) /
                  (e2m * Hypot(1d, tau) * Hypot(1d, taupa));
                tau += dtau;
                if (!(Abs(dtau) >= stol))
                    break;
            }

            return tau;
        }

        /// <summary>
        /// Evaluate the atan function with the result in degrees.
        /// </summary>
        /// <param name="x">in degrees</param>
        /// <returns>atan(<i>x</i>) in degrees.</returns>
        public static double Atand(double x) => Atan2d(x, 1d);

        /// <summary>
        /// The error-free sum of two numbers.
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <param name="t">the exact error given by (<paramref name="u"/> + <paramref name="v"/>) - <i>s</i>.</param>
        /// <returns><i>s</i> = round(<paramref name="u"/> + <paramref name="v"/>)</returns>
        /// <remarks>
        /// See D. E. Knuth, TAOCP, Vol 2, 4.2.2, Theorem B.  (Note that <paramref name="t"/> can be the same as one of the first two arguments.)
        /// </remarks>
        public static double Sum(double u, double v, out double t)
        {
            var s = u + v;
            var up = s - v;
            var vpp = s - up;

            up -= u;
            vpp -= v;
            t = -(up + vpp);

            // u + v =       s      + t
            //       = round(u + v) + t
            return s;
        }

        /// <summary>
        /// Normalize an angle.
        /// </summary>
        /// <param name="x">The angle in degrees.</param>
        /// <returns>The angle reduced to the range (-180°, 180°]</returns>
        public static double AngNormalize(double x)
        {
            x = IEEERemainder(x, 360d);
            return x != -180 ? x : 180;
        }

        /// <summary>
        /// The exact difference of two angles reduced to (-180°, 180°].
        /// </summary>
        /// <param name="x">The first angle in degrees.</param>
        /// <param name="y">The second angle in degrees.</param>
        /// <param name="e">The error term in degrees.</param>
        /// <returns><i>d</i>, the truncated value of <paramref name="y"/> - <paramref name="x"/>. </returns>
        /// <remarks>
        /// This computes <i>z</i> = <i>y</i> - <i>x</i> exactly, reduced to
        /// (-180°, 180°]; and then sets <i>z</i> = <i>d</i> + <i>e</i> where <i>d</i>
        /// is the nearest representable number to <i>z</i> and <i>e</i> is the truncation
        /// error.If <i>d</i> = -180, then <i>e</i> > 0; If <i>d</i> = 180, then <i>e</i> &lt;= 0.
        /// </remarks>
        public static double AngDiff(double x, double y, out double e)
        {
            double t, d = AngNormalize(Sum(IEEERemainder(-x, 360d),
                                           IEEERemainder(y, 360d), out t));

            // Here y - x = d + t (mod 360), exactly, where d is in (-180,180] and
            // abs(t) <= eps (eps = 2^-45 for doubles).  The only case where the
            // addition of t takes the result outside the range (-180,180] is d = 180
            // and t > 0.  The case, d = -180 + eps, t = -eps, can't happen, since
            // sum would have returned the exact result in such a case (i.e., given t
            // = 0).
            return Sum(d == 180 && t > 0 ? -180 : d, t, out e);
        }

        /// <summary>
        /// The exact difference of two angles reduced to (-180°, 180°].
        /// </summary>
        /// <param name="x">The first angle in degrees.</param>
        /// <param name="y">The second angle in degrees.</param>
        /// <returns><i>d</i>, the truncated value of <paramref name="y"/> - <paramref name="x"/>. </returns>
        /// <remarks>
        /// This computes <i>z</i> = <i>y</i> - <i>x</i> exactly, reduced to
        /// (-180°, 180°]; and then sets <i>z</i> = <i>d</i> + <i>e</i> where <i>d</i>
        /// is the nearest representable number to <i>z</i> and <i>e</i> is the truncation
        /// error.If <i>d</i> = -180, then <i>e</i> > 0; If <i>d</i> = 180, then <i>e</i> &lt;= 0.
        /// </remarks>
        public static double AngDiff(double x, double y) => AngDiff(x, y, out _);

        /// <summary>
        /// Normalize a latitude.
        /// </summary>
        /// <param name="x">The angle in degrees.</param>
        /// <returns><paramref name="x"/> if it is in the range [-90°, 90°], otherwise <see cref="double.NaN"/>.</returns>
        public static double LatFix(double x) => Abs(x) > 90 ? double.NaN : x;

        /// <summary>
        /// Round specified angle value.
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double AngRound(double x)
        {
            const double z = 1 / 16d;

            if (x is 0) return 0;
            var y = Abs(x);

            // The compiler mustn't "simplify" z - (z - y) to y
            y = y < z ? z - (z - y) : y;
            return x < 0 ? -y : y;
        }

        /// <summary>
        /// Evaluate a polynomial.
        /// </summary>
        /// <param name="N">the order of the polynomial.</param>
        /// <param name="p">the coefficient array (of size <paramref name="N"/> + 1).</param>
        /// <param name="x">the variable.</param>
        /// <returns>the value of the polynomial.</returns>
        /// <remarks>
        /// The evaluation uses Horner's method.
        /// </remarks>
        public static double PolyVal(int N, ReadOnlySpan<double> p, double x)
        {
            var i = 0;
            double y = N < 0 ? 0 : p[i++];
            while (--N >= 0) y = y * x + p[i++];
            return y;
        }

        /// <summary>
        /// Evaluate a polynomial.
        /// </summary>
        /// <param name="N">the order of the polynomial.</param>
        /// <param name="p">the coefficient array (of size <paramref name="N"/> + 1).</param>
        /// <param name="x">the variable.</param>
        /// <returns>the value of the polynomial.</returns>
        /// <remarks>
        /// The evaluation uses Horner's method.
        /// </remarks>
        public static double PolyVal(int N, ReadOnlyMemory<double> p, double x) => PolyVal(N, p.Span, x);
    }
}
