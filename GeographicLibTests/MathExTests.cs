using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using static GeographicLib.Tests.MathExTestData;

namespace GeographicLib.Tests
{
    [TestClass]
    public class MathExTests
    {
        [DataTestMethod]
        [DynamicData("Sum", typeof(MathExTestData))]
        public void TestSum_Sign(double u, double v, double r)
        {
            var r1 = MathEx.Sum(u, v, out _);
            Assert.That.EqualsExactly(r, r1);
        }

        [DataTestMethod]
        [DynamicData("Sind", typeof(MathExTestData))]
        public void TestSind_Sign(double x, double r)
        {
            var r1 = MathEx.Sind(x);
            Assert.That.EqualsExactly(r, r1);
        }

        [DataTestMethod]
        [DynamicData("Cosd", typeof(MathExTestData))]
        public void TestCosd_Sign(double x, double r)
        {
            var r1 = MathEx.Cosd(x);
            Assert.That.EqualsExactly(r, r1);
        }

        [DataTestMethod]
        [DynamicData("Tand", typeof(MathExTestData))]
        public void TestTand_Sign(double x, double r)
        {
            var r1 = MathEx.Tand(x);
            Assert.That.EqualsExactly(r, r1);
        }

        [DataTestMethod]
        [DynamicData("Atan2d", typeof(MathExTestData))]
        public void TestAtan2d_Sign(double y, double x, double r)
        {
            var r1 = MathEx.Atan2d(y, x);
            Assert.That.EqualsExactly(r, r1);
        }

        [DataTestMethod]
        [DynamicData("AngRound", typeof(MathExTestData))]
        public void TestAngRound_Sign(double x, double r)
        {
            var r1 = MathEx.AngRound(x);
            Assert.That.EqualsExactly(r, r1);
        }

        [DataTestMethod]
        [DynamicData("AngNormalize", typeof(MathExTestData))]
        public void TestAngNormalize_Sign(double x, double r)
        {
            var r1 = MathEx.AngNormalize(x);
            Assert.That.EqualsExactly(r, r1);
        }

        [DataTestMethod]
        [DynamicData("AngDiff", typeof(MathExTestData))]
        public void TestAngDiff_Sign(double x, double y, double r)
        {
            var r1 = MathEx.AngDiff(x, y);
            Assert.That.EqualsExactly(r, r1);
        }

        [DataTestMethod]
        [DynamicData("SinCosd", typeof(MathExTestData))]
        public void TestSinCosd_Sign(double x, double sinx, double cosx)
        {
            MathEx.SinCosd(x, out var sinx1, out var cosx1);
            Assert.That.EqualsExactly(sinx, sinx1);
            Assert.That.EqualsExactly(cosx, cosx1);
        }

        [TestMethod]
        public void TestSinCosd_Accuracy()
        {
            var sin9 = Math.Sin(9 * MathEx.Degree);

            MathEx.SinCosd(9, out var dsin9, out _);
            MathEx.SinCosd(81, out _, out var dcos81);
            MathEx.SinCosd(123456789, out var dsin123456789, out _);

            Assert.AreEqual(sin9, dsin9);
            Assert.AreEqual(dsin9, dcos81);
            Assert.AreEqual(dsin9, -dsin123456789);
        }

        [TestMethod]
        public void TestAtan2d_Accuracy()
        {
            var s = 7e-16;
            var (a, b) = (MathEx.Atan2d(s, -1), 180 - MathEx.Atan2d(s, 1));

            Assert.That.EqualsExactly(a, b);
        }

        [TestMethod]
        public void TestAngDiff_Accuracy()
        {
            double x = 138 + 128 * eps, y = -164;
            var r = MathEx.AngDiff(x, y);

            Assert.That.EqualsExactly(58 - 128 * eps, r);
        }

        [DataTestMethod]
        [DynamicData("Hypot", typeof(MathExTestData))]
        public void TestHypot(long x, long y, long expected)
        {
            var result = BitConverter.DoubleToInt64Bits(MathEx.Hypot(BitConverter.Int64BitsToDouble(x), BitConverter.Int64BitsToDouble(y)));
            var delta = Math.Abs(expected - result);
            if (delta > 1) Assert.Fail($"result={result}, expected={expected}, delta={delta}");
            else Assert.IsTrue(true);
        }

        [TestMethod]
        public void TestHypot_MonotonicityBug()
        {
            // hypot for Visual Studio (A=win32) fails monotonicity, e.g., with
            //   x  = 0.6102683302836215
            //   y1 = 0.7906090004346522
            //   y2 = y1 + 1e-16
            // the test
            //   hypot(x, y2) >= hypot(x, y1)
            // fails.  See also
            //   https://bugs.python.org/issue43088

            var x = 0.6102683302836215;
            var y1 = 0.7906090004346522;
            var y2 = y1 + 1e-16;

            // This may fail in 32-bit applications on Windows with UseManagedCMath set to false.
            Assert.IsTrue(MathEx.Hypot(x, y2) >= MathEx.Hypot(x, y1));
        }

        [TestMethod]
        public void TestHypot_GH75651()
        {
            var x = MathEx.Hypot(1, +1e20); // 1E+20 okA
            Assert.IsTrue(x > 0);
            x = MathEx.Hypot(1, -1e10); // 1E+10 ok
            Assert.IsTrue(x > 0);
            x = MathEx.Hypot(1, -1e20); // -1E+20 wrong
            Assert.IsTrue(x > 0);
        }

        [DataTestMethod]
        [DynamicData("Remquo", typeof(MathExTestData))]
        public void TestRemquo(long x, long y, long d, int q)
        {
            var r = MathEx.Remquo(BitConverter.Int64BitsToDouble(x), BitConverter.Int64BitsToDouble(y), out var rq);

            Assert.AreEqual(q, rq);
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(r), d);
        }

        [DataTestMethod]
        [DynamicData("Frexp", typeof(MathExTestData))]
        public void TestFrexp(long x, long expected, int e)
        {
            var result = BitConverter.DoubleToInt64Bits(MathEx.Frexp(BitConverter.Int64BitsToDouble(x), out var ae));
            Assert.AreEqual(e, ae);
            Assert.AreEqual(expected, result);
        }

        [DataTestMethod]
        [DynamicData("Log1p", typeof(MathExTestData))]
        public void TestLog1p(long x, long expected)
        {
            var result = BitConverter.DoubleToInt64Bits(MathEx.Log1p(BitConverter.Int64BitsToDouble(x)));
            var delta = Math.Abs(expected - result);
            if (delta > 1) Assert.Fail($"result={result}, expected={expected}, delta={delta}");
            else Assert.IsTrue(true);
            // Assert.AreEqual(expected, BitConverter.DoubleToInt64Bits(MathEx.Log1p(BitConverter.Int64BitsToDouble(x))));
        }

        [DataTestMethod]
        [DynamicData("Expm1", typeof(MathExTestData))]
        public void TestExpm1(long x, long expected)
        {
            var result = BitConverter.DoubleToInt64Bits(MathEx.Expm1(BitConverter.Int64BitsToDouble(x)));
            var delta = Math.Abs(expected - result);
            if (delta > 2) Assert.Fail($"result={result}, expected={expected}, delta={delta}");
            else Assert.IsTrue(true);
            // Assert.AreEqual(expected, BitConverter.DoubleToInt64Bits(MathEx.Log1p(BitConverter.Int64BitsToDouble(x))));
        }

        [DataTestMethod]
        [DataRow(4, 16)]
        [DataRow(0.5, 1.4142135623730951)]
        [DataRow(-4, 0.0625)]
        [DataRow(-0, 1.0)]
        [DataRow(0, 1.0)]
        [DataRow(double.NegativeInfinity, 0.0)]
        [DataRow(1024.0, double.PositiveInfinity)]
        [DataRow(-1075, 0)]
        [DataRow(double.PositiveInfinity, double.PositiveInfinity)]
        [DataRow(double.NaN, double.NaN)]
        [DynamicData("Exp2_Tiny", typeof(MathExTestData))]
        [DynamicData("Exp2_SmallPositive", typeof(MathExTestData))]
        [DynamicData("Exp2_SmallNegative", typeof(MathExTestData))]
        [DynamicData("Exp2_MediumPositive", typeof(MathExTestData))]
        [DynamicData("Exp2_MediumNegative", typeof(MathExTestData))]
        public void TestExp2(double input, double expected)
        {
            Assert.AreEqual(expected, MathEx.Exp2(input));
        }

#if !NET5_0_OR_GREATER
        [DataTestMethod]
        [DynamicData("CopySign", typeof(MathExTestData))]
        public void TestCopySign(long x, long y, long expected)
        {
            Assert.AreEqual(BitConverter.Int64BitsToDouble(expected), 
                MathEx.CopySign(BitConverter.Int64BitsToDouble(x), BitConverter.Int64BitsToDouble(y)));
        }

        [DataTestMethod]
        [DynamicData("ScaleB", typeof(MathExTestData))]
        public void TestScaleB(long x, int n, long expected)
        {
            Assert.AreEqual(BitConverter.Int64BitsToDouble(expected), MathEx.ScaleB(BitConverter.Int64BitsToDouble(x),n));
        }

        [DataTestMethod]
        [DynamicData("FusedMultiplyAdd", typeof(MathExTestData))]
        public void TestFusedMultiplyAdd(long x, long y,long z, long expected)
        {
            var result = BitConverter.DoubleToInt64Bits(MathEx.FusedMultiplyAdd(
                BitConverter.Int64BitsToDouble(x), BitConverter.Int64BitsToDouble(y), BitConverter.Int64BitsToDouble(z)));
            var delta = Math.Abs(expected - result);

            if (delta > 1) Assert.Fail($"result={result}, expected={expected}, delta={delta}");
            else Assert.IsTrue(true);
        }

        [DataTestMethod]
        [DynamicData("Log2", typeof(MathExTestData))]
        public void TestLog2(long x, long expected)
        {
            var result = BitConverter.DoubleToInt64Bits(MathEx.Log2(BitConverter.Int64BitsToDouble(x)));
            Assert.AreEqual(expected, result);
        }
#endif

#if NETCOREAPP2_1
        [DataTestMethod]
        [DynamicData("Cbrt", typeof(MathExTestData))]
        public void TestCbrt(long x, long expected)
        {
            var result = BitConverter.DoubleToInt64Bits(MathEx.Cbrt(BitConverter.Int64BitsToDouble(x)));
            var delta = Math.Abs(expected - result);
            if (delta > 1) Assert.Fail($"result={result}, expected={expected}, delta={delta}");
            else Assert.IsTrue(true);
        }

        [DataTestMethod]
        [DynamicData("Atanh", typeof(MathExTestData))]
        public void TestAtanh(long x, long expected)
        {
            var result = BitConverter.DoubleToInt64Bits(MathEx.Atanh(BitConverter.Int64BitsToDouble(x)));
            var delta = Math.Abs(expected - result);
            if (delta > 2) Assert.Fail($"result={result}, expected={expected}, delta={delta}");
            else Assert.IsTrue(true);
        }

        [DataTestMethod]
        [DynamicData("Asinh", typeof(MathExTestData))]
        public void TestAsinh(long x, long expected)
        {
            var result = BitConverter.DoubleToInt64Bits(MathEx.Asinh(BitConverter.Int64BitsToDouble(x)));
            var delta = Math.Abs(expected - result);
            if (delta > 1) Assert.Fail($"result={result}, expected={expected}, delta={delta}");
            else Assert.IsTrue(true);
        }
#endif
    }
}