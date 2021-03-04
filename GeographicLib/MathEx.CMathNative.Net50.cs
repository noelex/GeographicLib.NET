#if !CMATH_MANAGED && NET5_0
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
    static partial class MathEx
    {
        [ModuleInitializer]
        internal static void Init()
        {
            NativeLibrary.SetDllImportResolver(Assembly.GetExecutingAssembly(), (libraryName, assembly, searchPath) =>
            {
                // Use urcrtbase.dll as C runtime library on Windows.
                if (libraryName == "m")
                {
                    if (OperatingSystem.IsWindows())
                    {
                        return NativeLibrary.Load("ucrtbase.dll", assembly, searchPath);
                    }
                    else if (OperatingSystem.IsMacOS())
                    {
                        return NativeLibrary.Load("libSystem.B.dylib", assembly, searchPath);
                    }
                    else
                    {
                        return NativeLibrary.Load("libm.so.6", assembly, searchPath);
                    }
                }

                return IntPtr.Zero;
            });
        }

        /// <summary>
        /// Compute exp(x) - 1 without loss of precision for small values of x.
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        [DllImport("m", EntryPoint = "expm1")]
        public static extern double Expm1(double x);

        /// <summary>
        /// Computes the square root of the sum of the squares of x and y.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns>Hypotenuse of a right-angled triangle computed as √(x^2+y^2).</returns>
        [DllImport("m", EntryPoint = "hypot")]
        public static extern double Hypot(double x, double y);

        /// <summary>
        /// Compute log(1+x) without losing precision for small values of <paramref name="x"/>.
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        [DllImport("m", EntryPoint = "log1p")]
        public static extern double Log1p(double x);

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
        [DllImport("m", EntryPoint = "remquo")]
        public static extern double Remquo(double x, double y, out int quo);
    }
}
#endif