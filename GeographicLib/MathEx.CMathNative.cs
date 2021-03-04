#if !CMATH_MANAGED && !NET5_0
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib
{
    abstract class CMath
    {
        public abstract double Expm1(double x);

        public abstract double ScaleB(double number, int exp);

        public abstract double Remquo(double x, double y, out int quo);

        public abstract double Hypot(double x, double y);

        public abstract double Log1p(double x);

        public abstract double CopySign(double x, double y);

#if NETSTANDARD2_0
        public abstract double Atanh(double x);

        public abstract double Asinh(double x);

        public abstract double Cbrt(double x);
#endif

        public static CMath Current { get; } = OperatingSystem.IsWindows() ? (CMath)new CMathWindows() : new CMathUnix();

        private class CMathWindows : CMath
        {
            private const string CRuntimeLibrary = "ucrtbase";

            [DllImport(CRuntimeLibrary)]
            public static extern double expm1(double x);

            [DllImport(CRuntimeLibrary)]
            public static extern double scalbn(double number, int exp);

            [DllImport(CRuntimeLibrary)]
            public static extern double remquo(double x, double y, out int quo);

            [DllImport(CRuntimeLibrary)]
            public static extern double hypot(double x, double y);

            [DllImport(CRuntimeLibrary)]
            public static extern double log1p(double x);

            [DllImport(CRuntimeLibrary)]
            public static extern double copysign(double x, double y);

            public override double Expm1(double x) => expm1(x);

            public override double ScaleB(double number, int exp) => scalbn(number, exp);

            public override double Remquo(double x, double y, out int quo) => remquo(x, y, out quo);

            public override double Hypot(double x, double y) => hypot(x, y);

            public override double Log1p(double x) => log1p(x);

            public override double CopySign(double x, double y) => copysign(x, y);

#if NETSTANDARD2_0
            [DllImport(CRuntimeLibrary)]
            public static extern double atanh(double x);

            [DllImport(CRuntimeLibrary)]
            public static extern double asinh(double x);

            [DllImport(CRuntimeLibrary)]
            public static extern double cbrt(double x);

            public override double Atanh(double x) => atanh(x);

            public override double Asinh(double x) => asinh(x);

            public override double Cbrt(double x) => cbrt(x);
#endif
        }

        private class CMathUnix : CMath
        {
            private const string CRuntimeLibrary = "libm.so.6";

            [DllImport(CRuntimeLibrary)]
            public static extern double expm1(double x);

            [DllImport(CRuntimeLibrary)]
            public static extern double scalbn(double number, int exp);

            [DllImport(CRuntimeLibrary)]
            public static extern double remquo(double x, double y, out int quo);

            [DllImport(CRuntimeLibrary)]
            public static extern double hypot(double x, double y);

            [DllImport(CRuntimeLibrary)]
            public static extern double log1p(double x);

            [DllImport(CRuntimeLibrary)]
            public static extern double copysign(double x, double y);

            public override double Expm1(double x) => expm1(x);

            public override double ScaleB(double number, int exp) => scalbn(number, exp);

            public override double Remquo(double x, double y, out int quo) => remquo(x, y, out quo);

            public override double Hypot(double x, double y) => hypot(x, y);

            public override double Log1p(double x) => log1p(x);

            public override double CopySign(double x, double y) => copysign(x, y);

#if NETSTANDARD2_0
            [DllImport(CRuntimeLibrary)]
            public static extern double atanh(double x);

            [DllImport(CRuntimeLibrary)]
            public static extern double asinh(double x);

            [DllImport(CRuntimeLibrary)]
            public static extern double cbrt(double x);

            public override double Atanh(double x) => atanh(x);

            public override double Asinh(double x) => asinh(x);

            public override double Cbrt(double x) => cbrt(x);
#endif
        }
    }

    static partial class MathEx
    {
        /// <summary>
        /// Compute exp(x) - 1 without loss of precision for small values of x.
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double Expm1(double x) => CMath.Current.Expm1(x);

        /// <summary>
        /// Computes the square root of the sum of the squares of x and y.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns>Hypotenuse of a right-angled triangle computed as √(x^2+y^2).</returns>
        public static double Hypot(double x, double y) => CMath.Current.Hypot(x, y);

        /// <summary>
        /// Compute log(1+x) without losing precision for small values of <paramref name="x"/>.
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double Log1p(double x) => CMath.Current.Log1p(x);

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
        public static double Remquo(double x, double y, out int quo) => CMath.Current.Remquo(x, y, out quo);

        /// <summary>
        /// Returns x * 2^n computed efficiently.
        /// </summary>
        /// <param name="x">A double-precision floating-point number that specifies the base value.</param>
        /// <param name="n">A 32-bit integer that specifies the power.</param>
        /// <returns>x * 2^n computed efficiently.</returns>
        public static double ScaleB(double x, int n) => CMath.Current.ScaleB(x, n);

        /// <summary>
        /// Returns a value with the magnitude of <paramref name="x"/> and the sign of <paramref name="y"/>.
        /// </summary>
        /// <param name="x">A number whose magnitude is used in the result.</param>
        /// <param name="y">A number whose sign is the used in the result.</param>
        /// <returns>A value with the magnitude of <paramref name="x"/> and the sign of <paramref name="y"/>.</returns>
        public static double CopySign(double x, double y) => CMath.Current.CopySign(x, y);

#if NETSTANDARD2_0
        /// <summary>
        /// Returns the angle whose hyperbolic tangent is the specified number.
        /// </summary>
        /// <param name="x">A number representing a hyperbolic tangent, where d must be greater than or equal to -1, but less than or equal to 1.</param>
        /// <returns>
        /// An angle, θ, measured in radians, such that -∞ &lt; θ &lt; -1, or 1 &lt; θ &lt; ∞.
        /// -or- <see cref="double.NaN"/> if <paramref name="x"/> &lt; -1 or <paramref name="x"/> > 1 or <paramref name="x"/> equals <see cref="double.NaN"/>.
        /// </returns>
        public static double Atanh(double x) => CMath.Current.Atanh(x);

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
        public static double Asinh(double x) => CMath.Current.Asinh(x);

        /// <summary>
        /// Returns the cube root of a specified number.
        /// </summary>
        /// <param name="x">The number whose cube root is to be found.</param>
        /// <returns>The cube root of <paramref name="x"/>. -or- <see cref="double.NaN"/> if <paramref name="x"/> equals <see cref="double.NaN"/>.</returns>
        public static double Cbrt(double x) => CMath.Current.Cbrt(x);
#endif
    }
}
#endif