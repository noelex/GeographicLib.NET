using System;
using System.Runtime.CompilerServices;
using static GeographicLib.Macros;
using static System.Math;

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

        internal const float FLT_MIN = 1.17549435082228750797e-38F;

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

        /// <summary>
        /// # of decimal digits of precision
        /// </summary>
        internal const int DBL_DIG = 10;

#if NETSTANDARD2_0
        /// <summary>
        /// Gets or sets a value representing that whether <see cref="MathEx"/> should use managed implementations of C mathematical functions
        /// when there's no corresponding implementation provided by .NET runtime.
        /// </summary>
        /// <remarks>
        /// The following functions have managed implementation:
        /// <list type="bullet">
        /// <item><see cref="Atanh(double)"/></item>
        /// <item><see cref="Asinh(double)"/></item>
        /// <item><see cref="Cbrt(double)"/></item>
        /// <item><see cref="ScaleB(double, int)"/></item>
        /// <item><see cref="CopySign(double, double)"/></item>
        /// <item><see cref="Log2(double)"/></item>
        /// <item><see cref="Expm1(double)"/></item>
        /// <item><see cref="Log1p(double)"/></item>
        /// <item><see cref="Hypot(double, double)"/></item>
        /// <item><see cref="Remquo(double, double, out int)"/></item>
        /// <item><see cref="Frexp(double, out int)"/></item>
        /// </list>
        /// When set to <see langword="true"/>, <see cref="MathEx"/> will use managed implementations when the above functions are called.
        /// <para>
        /// When set to <see langword="false"/>, <see cref="MathEx"/> will use platform dependent implementations provided by system C runtime
        /// library when the above functions are called.
        /// </para>
        /// </remarks>
#elif NETSTANDARD2_1
        /// <summary>
        /// Gets or sets a value representing that whether <see cref="MathEx"/> should use managed implementations of C mathematical functions
        /// when there's no corresponding implementation provided by .NET runtime.
        /// </summary>
        /// <remarks>
        /// The following functions have managed implementation:
        /// <list type="bullet">
        /// <item><see cref="ScaleB(double, int)"/></item>
        /// <item><see cref="CopySign(double, double)"/></item>
        /// <item><see cref="Log2(double)"/></item>
        /// <item><see cref="Expm1(double)"/></item>
        /// <item><see cref="Log1p(double)"/></item>
        /// <item><see cref="Hypot(double, double)"/></item>
        /// <item><see cref="Remquo(double, double, out int)"/></item>
        /// <item><see cref="Frexp(double, out int)"/></item>
        /// </list>
        /// When set to <see langword="true"/>, <see cref="MathEx"/> will use managed implementations when the above functions are called.
        /// <para>
        /// When set to <see langword="false"/>, <see cref="MathEx"/> will use platform dependent implementations provided by system C runtime
        /// library when the above functions are called.
        /// </para>
        /// </remarks>
#else
        /// <summary>
        /// Gets or sets a value representing that whether <see cref="MathEx"/> should use managed implementations of C mathematical functions
        /// when there's no corresponding implementation provided by .NET runtime.
        /// </summary>
        /// <remarks>
        /// The following functions have managed implementation:
        /// <list type="bullet">
        /// <item><see cref="Expm1(double)"/></item>
        /// <item><see cref="Log1p(double)"/></item>
        /// <item><see cref="Hypot(double, double)"/></item>
        /// <item><see cref="Remquo(double, double, out int)"/></item>
        /// <item><see cref="Frexp(double, out int)"/></item>
        /// </list>
        /// When set to <see langword="true"/>, <see cref="MathEx"/> will use managed implementations when the above functions are called.
        /// <para>
        /// When set to <see langword="false"/>, <see cref="MathEx"/> will use platform dependent implementations provided by system C runtime
        /// library when the above functions are called.
        /// </para>
        /// </remarks>
#endif
        public static bool UseManagedCMath { get => CMath.UseManagedImplementation; set => CMath.UseManagedImplementation = value; }

#if NETSTANDARD2_0
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static bool IsFinite(double x) => ((ulong)BitConverter.DoubleToInt64Bits(x) & unchecked((ulong)-1L) >> 1) < 0x7ffUL << 52;

        /// <summary>
        /// Returns the angle whose hyperbolic tangent is the specified number.
        /// </summary>
        /// <param name="x">A number representing a hyperbolic tangent, where d must be greater than or equal to -1, but less than or equal to 1.</param>
        /// <returns>
        /// An angle, θ, measured in radians, such that -∞ &lt; θ &lt; -1, or 1 &lt; θ &lt; ∞.
        /// -or- <see cref="double.NaN"/> if <paramref name="x"/> &lt; -1 or <paramref name="x"/> > 1 or <paramref name="x"/> equals <see cref="double.NaN"/>.
        /// </returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)] 
        public static double Atanh(double x) => CMath.Instance.Atanh(x);

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
        [MethodImpl(MethodImplOptions.AggressiveInlining)] 
        public static double Asinh(double x) => CMath.Instance.Asinh(x);

        /// <summary>
        /// Returns the cube root of a specified number.
        /// </summary>
        /// <param name="x">The number whose cube root is to be found.</param>
        /// <returns>The cube root of <paramref name="x"/>. -or- <see cref="double.NaN"/> if <paramref name="x"/> equals <see cref="double.NaN"/>.</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)] 
        public static double Cbrt(double x) => CMath.Instance.Cbrt(x);
#else
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static bool IsFinite(double x) => double.IsFinite(x);
#endif

#if !NET5_0_OR_GREATER
        /// <summary>
        /// Returns x * 2^n computed efficiently.
        /// </summary>
        /// <param name="x">A double-precision floating-point number that specifies the base value.</param>
        /// <param name="n">A 32-bit integer that specifies the power.</param>
        /// <returns>x * 2^n computed efficiently.</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)] 
        public static double ScaleB(double x, int n) => CMath.Instance.ScaleB(x, n);

        /// <summary>
        /// Returns a value with the magnitude of <paramref name="x"/> and the sign of <paramref name="y"/>.
        /// </summary>
        /// <param name="x">A number whose magnitude is used in the result.</param>
        /// <param name="y">A number whose sign is the used in the result.</param>
        /// <returns>A value with the magnitude of <paramref name="x"/> and the sign of <paramref name="y"/>.</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)] 
        public static double CopySign(double x, double y) => CMath.Instance.CopySign(x, y);

        /// <summary>
        /// Returns (x * y) + z, rounded as one ternary operation.
        /// </summary>
        /// <param name="x">The number to be multiplied with y.</param>
        /// <param name="y">The number to be multiplied with x.</param>
        /// <param name="z">The number to be added to the result of x multiplied by y.</param>
        /// <returns>(x * y) + z, rounded as one ternary operation.</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double FusedMultiplyAdd(double x, double y, double z) => CMath.Instance.FusedMultiplyAdd(x, y, z);

        /// <summary>
        /// Returns the base 2 logarithm of a specified number.
        /// </summary>
        /// <param name="x">A number whose logarithm is to be found.</param>
        /// <returns>
        /// One of the values in the following table.
        /// <list type="table">
        /// <item><paramref name="x"/> parameter – Return value</item>
        /// <item>Positive – The base 2 log of <paramref name="x"/>; that is, log 2<paramref name="x"/>.</item>
        /// <item>Zero - <see cref="double.NegativeInfinity"/></item>
        /// <item>Negative - <see cref="double.NaN"/></item>
        /// <item>Equal to <see cref="double.NaN"/> - <see cref="double.NaN"/></item>
        /// <item>Equal to <see cref="double.PositiveInfinity"/> - <see cref="double.PositiveInfinity"/></item>
        /// </list>
        /// </returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Log2(double x) => CMath.Instance.Log2(x);
#endif

        /// <summary>
        /// Compute exp(x) - 1 without loss of precision for small values of x.
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Expm1(double x) => CMath.Instance.Expm1(x);

        /// <summary>
        /// Computes the square root of the sum of the squares of x and y.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns>Hypotenuse of a right-angled triangle computed as √(x^2+y^2).</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Hypot(double x, double y) => CMath.Instance.Hypot(x, y);

        /// <summary>
        /// Compute log(1+x) without losing precision for small values of <paramref name="x"/>.
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Log1p(double x) => CMath.Instance.Log1p(x);

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
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Remquo(double x, double y, out int quo) => CMath.Instance.Remquo(x, y, out quo);

        /// <summary>
        /// Decomposes given floating point value <paramref name="x"/> into a normalized fraction and an integral power of two.
        /// </summary>
        /// <param name="x">Floating point value.</param>
        /// <param name="e">Pointer to integer value to store the exponent to.</param>
        /// <returns>
        /// <para>
        /// If <paramref name="x"/> is zero, returns zero and stores zero in <paramref name="e"/>.
        /// </para>
        /// <para>
        /// Otherwise (if <paramref name="x"/> is not zero), if no errors occur,
        /// returns the value x in the range (-1;-0.5], [0.5; 1) and stores an integer value in <paramref name="e"/> such that
        /// x×2(<paramref name="e"/>)=<paramref name="x"/>
        /// </para>
        /// <para>
        /// If the value to be stored in <paramref name="e"/> is outside the range of int, the behavior is unspecified.
        /// </para>
        /// <para>
        /// If <paramref name="x"/> is not a floating-point number, the behavior is unspecified.
        /// </para>
        /// </returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Frexp(double x, out int e) => CMath.Instance.Frexp(x, out e);

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
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
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
        /// <para>
        /// If <i>x</i> = -0 or a negative multiple of 180°, then sin<i>x</i> = -0;
        /// this is the only case where -0 is returned.
        /// </para>
        /// </remarks>
        public static void SinCosd(double x, out double sinx, out double cosx)
        {
            // In order to minimize round-off errors, this function exactly reduces
            // the argument to the range [-45, 45] before converting it to radians.
            var r = Remquo(x, 90, out var q);
            r *= Degree;

#if NET6_0_OR_GREATER
            // TODO: Math.SinCos is failing test case PlanimeterTests.TestComputeWithPoints on Windows,
            // procuding error ranging from 0.0002 to 0.02 meters, which is much greater than expected error 1e-8.
            // Needs further investigation.
            // var (s, c) = SinCos(r);
            var (s, c) = (Sin(r), Cos(r));
#else
            double s = Sin(r), c = Cos(r);
#endif

            switch ((uint)q & 3U)
            {
                case 0U: sinx = s; cosx = c; break;
                case 1U: sinx = c; cosx = -s; break;
                case 2U: sinx = -s; cosx = -c; break;
                default: sinx = -c; cosx = s; break; // case 3U
            }

            // http://www.open-std.org/jtc1/sc22/wg14/www/docs/n1950.pdf
            cosx += 0d;                              // special values from F.10.1.12
            if (sinx == 0) sinx = CopySign(sinx, x); // special values from F.10.1.13
        }

        /// <summary>
        /// Evaluate the sine function with the argument in degrees.
        /// </summary>
        /// <param name="x">in degrees.</param>
        /// <returns>sin(<i>x</i>)</returns>
        /// <remarks>
        /// The result is +0 for <paramref name="x"/> = +0 and positive multiples of 180°. The
        /// result is -0 for <paramref name="x"/> = -0 and negative multiples of 180°.
        /// </remarks>
        public static double Sind(double x)
        {
            // See sincosd
            double r;
            r = Remquo(x, 90d, out var q); // now abs(r) <= 45
            r *= Degree;

            var p = (uint)q;

            r = (p & 1U) != 0 ? Cos(r) : Sin(r);
            if ((p & 2U) != 0) r = -r;

            if (r == 0) r = CopySign(r, x);

            return r;
        }

        /// <summary>
        /// Evaluate the cosine function with the argument in degrees.
        /// </summary>
        /// <param name="x">in degrees.</param>
        /// <returns>cos(<i>x</i>)</returns>
        /// <remarks>
        /// The result is +0 for <paramref name="x"/> an odd multiple of 90°.
        /// </remarks>
        public static double Cosd(double x)
        {
            // See sincosd
            double r;
            r = Remquo(x, 90d, out var q); // now abs(r) <= 45
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
        /// If <paramref name="x"/> is an odd multiple of 90°, then a suitably large (but finite) value is returned.
        /// </remarks>
        public static double Tand(double x)
        {
            const double overflow = 1 / (DBL_EPSILON * DBL_EPSILON);

            SinCosd(x, out var s, out var c);
            // http://www.open-std.org/jtc1/sc22/wg14/www/docs/n1950.pdf
            var r = s / c;  // special values from F.10.1.14
                            // With C++17 this becomes clamp(s / c, -overflow, overflow);
                            // Use max/min here (instead of fmax/fmin) to preserve NaN
            return Min(Max(r, -overflow), overflow);
        }

        /// <summary>
        /// Evaluate the atan2 function with the result in degrees.
        /// </summary>
        /// <param name="y"></param>
        /// <param name="x"></param>
        /// <returns>atan2(<i>y</i>, <i>x</i>) in degrees.</returns>
        /// <remarks>
        /// The result is in the range [-180°, 180°].  N.B.,
        /// atan2d(±0, -1) = ±180°.
        /// </remarks>
        public static double Atan2d(double y, double x)
        {
            // In order to minimize round-off errors, this function rearranges the
            // arguments so that result of atan2 is in the range [-pi/4, pi/4] before
            // converting it to degrees and mapping the result to the correct
            // quadrant.
            int q = 0;
            if (Abs(y) > Abs(x)) { Swap(ref x, ref y); q = 2; }
            if (SignBit(x)) { x = -x; ++q; }
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
                case 1: ang = (SignBit(y) ? -180 : 180) - ang; break;
                case 2: ang = 90 - ang; break;
                case 3: ang = -90 + ang; break;
                default: break;
            }

            return ang;
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
            // if s = 0, then t = 0 and give t the same sign as s
            up -= u;
            vpp -= v;
            t = s != 0 ? 0d - (up + vpp) : s;

            // u + v =       s      + t
            //       = round(u + v) + t
            return s;
        }

        /// <summary>
        /// Normalize an angle.
        /// </summary>
        /// <param name="x">The angle in degrees.</param>
        /// <returns>The angle reduced to the range [-180°, 180°]</returns>
        public static double AngNormalize(double x)
        {
            var y = IEEERemainder(x, 360d);
            return Abs(y) == 180 ? CopySign(180, x) : y;
        }

        /// <summary>
        /// The exact difference of two angles reduced to [-180°, 180°].
        /// </summary>
        /// <param name="x">The first angle in degrees.</param>
        /// <param name="y">The second angle in degrees.</param>
        /// <param name="e">The error term in degrees.</param>
        /// <returns><i>d</i>, the truncated value of <paramref name="y"/> - <paramref name="x"/>. </returns>
        /// <remarks>
        /// This computes <i>z</i> = <i>y</i> - <i>x</i> exactly, reduced to
        /// [-180°, 180°]; and then sets <i>z</i> = <i>d</i> + <i>e</i> where <i>d</i>
        /// is the nearest representable number to <i>z</i> and <i>e</i> is the truncation
        /// error. If <i>z</i> = ±0° or ±180°, then the sign of <i>d</i> is given
        /// by the sign of <i>y</i> - <i>x</i>. The maximum absolute value of <i>e</i> is 2^-26
        /// (for doubles).
        /// </remarks>
        public static double AngDiff(double x, double y, out double e)
        {
            // Use remainder instead of AngNormalize, since we treat boundary cases
            // later taking account of the error
            var d = IEEERemainder(Sum(IEEERemainder(-x, 360),
                                IEEERemainder(y, 360), out e), 360);
            // This second sum can only change d if abs(d) < 128, so don't need to
            // apply remainder yet again.
            d = Sum(d, e, out e);
            // Fix the sign if d = -180, 0, 180.
            if (d == 0 || Abs(d) == 180)
                // If e == 0, take sign from y - x
                // else (e != 0, implies d = +/-180), d and e must have opposite signs
                d = CopySign(d, e == 0 ? y - x : -e);
            return d;
        }

        /// <summary>
        /// The exact difference of two angles reduced to [-180°, 180°].
        /// </summary>
        /// <param name="x">The first angle in degrees.</param>
        /// <param name="y">The second angle in degrees.</param>
        /// <returns><i>d</i>, the truncated value of <paramref name="y"/> - <paramref name="x"/>. </returns>
        /// <remarks>
        /// This computes <i>z</i> = <i>y</i> - <i>x</i> exactly, reducing
        /// it to [-180°, 180°] and rounding the result.
        /// </remarks>
        public static double AngDiff(double x, double y) => AngDiff(x, y, out _);

        /// <summary>
        /// Normalize a latitude.
        /// </summary>
        /// <param name="x">The angle in degrees.</param>
        /// <returns><paramref name="x"/> if it is in the range [-90°, 90°], otherwise <see cref="double.NaN"/>.</returns>
        public static double LatFix(double x) => Abs(x) > 90 ? double.NaN : x;

        /// <summary>
        /// Coarsen a value close to zero.
        /// </summary>
        /// <param name="x"></param>
        /// <returns>the coarsened value.</returns>
        /// <remarks>
        /// The makes the smallest gap in <i>x</i> = 1/16 - nextafter(1/16, 0) = 1/2^57
        /// for doubles = 0.8 pm on the earth if <i>x</i> is an angle
        /// in degrees.  (This is about 2000 times more resolution than we get with
        /// angles around 90°.)  We use this to avoid having to deal with near
        /// singular cases when <i>x</i> is non-zero but tiny (e.g.,
        /// 10^-200.  The sign of ±0 is preserved.
        /// </remarks>
        public static double AngRound(double x)
        {
            const double z = 1 / 16d;

            var y = Abs(x);

            // The compiler mustn't "simplify" z - (z - y) to y
            y = y < z ? z - (z - y) : y;
            return CopySign(y, x);
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

        private static readonly long SignMask = BitConverter.DoubleToInt64Bits(-0.0) ^ BitConverter.DoubleToInt64Bits(+0.0);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static bool SignBit(double v) => (BitConverter.DoubleToInt64Bits(v) & SignMask) == SignMask;
    }
}
