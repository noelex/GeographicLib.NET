using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static GeographicLib.Tests.MathExTestData;

namespace GeographicLib.Tests
{
    [TestClass]
    public class DMSTests
    {
        [DataTestMethod]
        [DataRow(33.3, 44.4, "33d18'00.0\"N 044d24'00.0\"E")]
        [DataRow(0.0, 123.456, "00d00'00.0\"N 123d27'21.6\"E")]
        public void TestEncode(double lat, double lon, string dms)
        {
            var latStr = DMS.Encode(lat, 5, HemisphereIndicator.Latitude);
            var lonStr = DMS.Encode(lon, 5, HemisphereIndicator.Longitude);
            Assert.AreEqual(dms, $"{latStr} {lonStr}");
        }

        [DataTestMethod]
        [DataRow(33.3, 44.4, "33d18'00.0\"N 044d24'00.0\"E")]
        [DataRow(33.3, 44.4, "33°18'00.0\"N 044°24'00.0\"E")]
        [DataRow(0.0, 123.456, "00d00'00.0\"N 123d27'21.6\"E")]
        public void TestDecode(double lat, double lon, string dms)
        {
            var components = dms.Split(' ');
            var (latVal, lonVal) = DMS.Decode(components[0], components[1]);
            Assert.AreEqual(lat, latVal, 0.001);
            Assert.AreEqual(lon, lonVal, 0.001);
        }

        [DataTestMethod]
        [DataRow("+0", +0.0)]
        [DataRow("-0", -0.0)]
        [DataRow("nan", nan)]
        [DataRow("+inf", +inf)]
        [DataRow("inf", +inf)]
        [DataRow("-inf", -inf)]
        [DataRow("+0N", +0.0)]
        [DataRow("-0N", -0.0)]
        [DataRow("+0S", -0.0)]
        [DataRow("-0S", +0.0)]
        public void TestSign(string input, double r)
        {
            var (r1, _) = DMS.Decode(input);
            Assert.That.EqualsExactly(r, r1);
        }

        [DataTestMethod]
        [DataRow(nan, TrailingUnit.Degree, 0, "nan")]
        [DataRow(-inf, TrailingUnit.Degree, 0, "-inf")]
        [DataRow(+inf, TrailingUnit.Degree, 0, "inf")]
        [DataRow(-3.5, TrailingUnit.Degree, 0, "-4")]
        [DataRow(-2.5, TrailingUnit.Degree, 0, "-2")]
        [DataRow(-1.5, TrailingUnit.Degree, 0, "-2")]
        [DataRow(-0.5, TrailingUnit.Degree, 0, "-0")]
        [DataRow(-0.0, TrailingUnit.Degree, 0, "-0")]
        [DataRow(+0.0, TrailingUnit.Degree, 0, "0")]
        [DataRow(+0.5, TrailingUnit.Degree, 0, "0")]
        [DataRow(+1.5, TrailingUnit.Degree, 0, "2")]
        [DataRow(+2.5, TrailingUnit.Degree, 0, "2")]
        [DataRow(+3.5, TrailingUnit.Degree, 0, "4")]
        [DataRow(-1.75, TrailingUnit.Degree, 1, "-1.8")]
        [DataRow(-1.25, TrailingUnit.Degree, 1, "-1.2")]
        [DataRow(-0.75, TrailingUnit.Degree, 1, "-0.8")]
        [DataRow(-0.25, TrailingUnit.Degree, 1, "-0.2")]
        [DataRow(-0.0, TrailingUnit.Degree, 1, "-0.0")]
        [DataRow(+0.0, TrailingUnit.Degree, 1, "0.0")]
        [DataRow(+0.25, TrailingUnit.Degree, 1, "0.2")]
        [DataRow(+0.75, TrailingUnit.Degree, 1, "0.8")]
        [DataRow(+1.25, TrailingUnit.Degree, 1, "1.2")]
        [DataRow(+1.75, TrailingUnit.Degree, 1, "1.8")]
        [DataRow(1e20, TrailingUnit.Degree, 0, "100000000000000000000")]
        [DataRow(1e21, TrailingUnit.Degree, 0, "1000000000000000000000")]
        public void TestEncode_RoundToEven(double angle, TrailingUnit trailing, int prec, string r)
        {
            var r1 = DMS.Encode(angle, trailing, prec);
            Assert.AreEqual(r, r1);
        }
    }
}
