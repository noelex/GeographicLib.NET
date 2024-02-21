namespace GeographicLib
{
    abstract class CMath
    {
        private static readonly CMath managed = new CMathManaged(), native = CMathNative.Create();

        public abstract double Remquo(double x, double y, out int quo);

        public abstract double Hypot(double x, double y);

        public abstract double Frexp(double x, out int e);

        public abstract double Log1p(double x);

        public abstract double Expm1(double x);

        public abstract double Exp2(double x);

#if !NET5_0_OR_GREATER
        public abstract double CopySign(double x, double y);

        public abstract double ScaleB(double number, int exp);

        public abstract double FusedMultiplyAdd(double x, double y, double z);

        public abstract double Log2(double x);
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
