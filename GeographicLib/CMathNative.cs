﻿using System;
using System.Reflection;
using System.Runtime.InteropServices;

namespace GeographicLib
{
#if NET5_0_OR_GREATER
    class CMathNative : CMath
    {
        static CMathNative()
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

        [SuppressGCTransition, DllImport("m")]
        private static extern double hypot(double x, double y);

        [SuppressGCTransition, DllImport("m")]
        private static extern double remquo(double x, double y, out int quo);

        [SuppressGCTransition, DllImport("m")]
        private static extern double expm1(double x);

        [SuppressGCTransition, DllImport("m")]
        private static extern double log1p(double x);

        [SuppressGCTransition, DllImport("m")]
        private static extern double exp2(double x);

        [SuppressGCTransition, DllImport("m")]
        private static extern double frexp(double x, out int e);

        public override double Hypot(double x, double y) => hypot(x, y);

        public override double Remquo(double x, double y, out int quo) => remquo(x, y, out quo);

        public override double Expm1(double x) => expm1(x);

        public override double Log1p(double x) => log1p(x);

        public override double Exp2(double x) => exp2(x);

        public override double Frexp(double x, out int e) => frexp(x, out e);

        public static CMath Create() => new CMathNative();
    }
#else
    class CMathNative
    {
        public static CMath Create()
        {
            if (OperatingSystem.IsWindows())
            {
                return new CMathWindows();
            }
            else if (OperatingSystem.IsMacOS())
            {
                return new CMathOSX();
            }

            return new CMathLinux();
        }

        private class CMathWindows : CMath
        {
            private const string CRuntimeLibrary = "ucrtbase.dll";

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
            public static extern double exp2(double x);

            [DllImport(CRuntimeLibrary)]
            public static extern double copysign(double x, double y);

            [DllImport(CRuntimeLibrary)]
            public static extern double fma(double x, double y, double z);

            [DllImport(CRuntimeLibrary)]
            private static extern double frexp(double x, out int e);

            [DllImport(CRuntimeLibrary)]
            public static extern double log2(double x);

            public override double Expm1(double x) => expm1(x);

            public override double ScaleB(double number, int exp) => scalbn(number, exp);

            public override double Remquo(double x, double y, out int quo) => remquo(x, y, out quo);

            public override double Hypot(double x, double y) => hypot(x, y);

            public override double Log1p(double x) => log1p(x);

            public override double CopySign(double x, double y) => copysign(x, y);

            public override double FusedMultiplyAdd(double x, double y, double z) => fma(x, y, z);

            public override double Frexp(double x, out int e) => frexp(x, out e);

            public override double Log2(double x) => log2(x);

            public override double Exp2(double x) => exp2(x);

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

        private class CMathLinux : CMath
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
            public static extern double exp2(double x);

            [DllImport(CRuntimeLibrary)]
            public static extern double copysign(double x, double y);

            [DllImport(CRuntimeLibrary)]
            public static extern double fma(double x, double y, double z);

            [DllImport(CRuntimeLibrary)]
            private static extern double frexp(double x, out int e);

            [DllImport(CRuntimeLibrary)]
            public static extern double log2(double x);

            public override double Expm1(double x) => expm1(x);

            public override double ScaleB(double number, int exp) => scalbn(number, exp);

            public override double Remquo(double x, double y, out int quo) => remquo(x, y, out quo);

            public override double Hypot(double x, double y) => hypot(x, y);

            public override double Log1p(double x) => log1p(x);

            public override double CopySign(double x, double y) => copysign(x, y);

            public override double FusedMultiplyAdd(double x, double y, double z) => fma(x, y, z);

            public override double Frexp(double x, out int e) => frexp(x, out e);

            public override double Log2(double x) => log2(x);

            public override double Exp2(double x) => exp2(x);

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

        private class CMathOSX : CMath
        {
            private const string CRuntimeLibrary = "libSystem.B.dylib";

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

            [DllImport(CRuntimeLibrary)]
            public static extern double fma(double x, double y, double z);

            [DllImport(CRuntimeLibrary)]
            private static extern double frexp(double x, out int e);

            [DllImport(CRuntimeLibrary)]
            public static extern double log2(double x);

            [DllImport(CRuntimeLibrary)]
            public static extern double exp2(double x);

            public override double Expm1(double x) => expm1(x);

            public override double ScaleB(double number, int exp) => scalbn(number, exp);

            public override double Remquo(double x, double y, out int quo) => remquo(x, y, out quo);

            public override double Hypot(double x, double y) => hypot(x, y);

            public override double Log1p(double x) => log1p(x);

            public override double CopySign(double x, double y) => copysign(x, y);

            public override double FusedMultiplyAdd(double x, double y, double z) => fma(x, y, z);

            public override double Frexp(double x, out int e) => frexp(x, out e);

            public override double Log2(double x) => log2(x);

            public override double Exp2(double x) => exp2(x);

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
#endif
}
