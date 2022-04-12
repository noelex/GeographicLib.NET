using Microsoft.VisualStudio.TestTools.UnitTesting;
using GeographicLib;
using System;
using System.Collections.Generic;
using System.Text;
using System.Runtime.InteropServices;

namespace GeographicLib.Tests
{
    [TestClass]
    public class MathExTests
    {
        [DataTestMethod]
        [DataRow(+0.0, -0.0, +1.0)]
        [DataRow(-0.0, -0.0, -1.0)]
        [DataRow(+0.0, +0.0, +0.0)]
        [DataRow(-0.0, +0.0, -0.0)]
        [DataRow(+0.0, -1.0, +1.0)]
        [DataRow(-0.0, -1.0, -1.0)]
        [DataRow(+0.0, +1.0, +0.0)]
        [DataRow(-0.0, +1.0, -0.0)]
        [DataRow(-1.0, +0.0, -0.5)]
        [DataRow(-1.0, -0.0, -0.5)]
        [DataRow(+1.0, +0.0, +0.5)]
        [DataRow(+1.0, -0.0, +0.5)]
        [DataRow(+1.0, double.NegativeInfinity, +1.0)]
        [DataRow(-1.0, double.NegativeInfinity, -1.0)]
        [DataRow(+1.0, double.PositiveInfinity, +0.0)]
        [DataRow(-1.0, double.PositiveInfinity, -0.0)]
        [DataRow(double.PositiveInfinity, +1.0, +0.5)]
        [DataRow(double.PositiveInfinity, -1.0, +0.5)]
        [DataRow(double.NegativeInfinity, +1.0, -0.5)]
        [DataRow(double.NegativeInfinity, -1.0, -0.5)]
        [DataRow(double.PositiveInfinity, double.NegativeInfinity, +0.75)]
        [DataRow(double.NegativeInfinity, double.NegativeInfinity, -0.75)]
        [DataRow(double.PositiveInfinity, double.PositiveInfinity, +0.25)]
        [DataRow(double.NegativeInfinity, double.PositiveInfinity, -0.25)]
        [DataRow(double.NaN, +1.0, double.NaN)]
        [DataRow(+1.0, double.NaN, double.NaN)]
        public void TestAtan(double y, double x, double r)
        {
            r *= 180;
            var r1 = MathEx.Atan2d(y, x);

            if (double.IsNaN(r))
            {
                Assert.IsTrue(double.IsNaN(r1));
            }
            else
            {
                Assert.AreEqual(r, r1);
                Assert.AreEqual(Math.Sign(r), Math.Sign(r1));
            }
        }

        [TestMethod]
        public void TestSinCosd()
        {
            var sin9 = Math.Sin(9 * MathEx.Degree);

            double dsin9, dcos81, dsin123456789;
            MathEx.SinCosd(9, out dsin9, out _);
            MathEx.SinCosd(81, out _, out dcos81);
            MathEx.SinCosd(123456789, out dsin123456789, out _);

            Assert.AreEqual(sin9, dsin9);
            Assert.AreEqual(dsin9, dcos81);
            Assert.AreEqual(dsin9, -dsin123456789);
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

        [DataTestMethod]
        [DynamicData("Remquo", typeof(MathExTestData))]
        public void TestRemquo(long x, long y, long d, int q)
        {
            var r = MathEx.Remquo(BitConverter.Int64BitsToDouble(x), BitConverter.Int64BitsToDouble(y), out var rq);

            Assert.AreEqual(q, rq);
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(r), d);
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
        [DynamicData("Frexp", typeof(MathExTestData))]
        public void TestFrexp(long x, long expected, int e)
        {
            var result = BitConverter.DoubleToInt64Bits(MathEx.Frexp(BitConverter.Int64BitsToDouble(x),out var ae ));
            Assert.AreEqual(e, ae);
            Assert.AreEqual(expected, result);
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