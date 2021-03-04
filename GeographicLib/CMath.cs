using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib
{
    abstract class CMath
    {
        private static readonly CMath managed = new CMathManaged(), native = CMathNative.Create();

        public abstract double Expm1(double x);

        public abstract double Remquo(double x, double y, out int quo);

        public abstract double Hypot(double x, double y);

        public abstract double Log1p(double x);

#if !NET5_0
        public abstract double CopySign(double x, double y);

        public abstract double ScaleB(double number, int exp);
#endif

#if NETSTANDARD2_0
        public abstract double Atanh(double x);

        public abstract double Asinh(double x);

        public abstract double Cbrt(double x);
#endif
        public static CMath Instance { get; private set; } = managed;

        public static bool UseManagedImplementation
        {
            get => Instance == managed;
            set { Instance = value ? managed : native; }
        }
    }
}
