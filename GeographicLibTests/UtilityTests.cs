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
    public class UtilityTest
    {
        [DataTestMethod]
        [DataRow("+0",   +0.0)]
        [DataRow("-0",   -0.0)]
        [DataRow("nan",  nan)]
        [DataRow("+inf", +inf)]
        [DataRow("inf",  +inf)]
        [DataRow("-inf", -inf)]
        public void TestParseDouble(string input, double r)
        {
            var r1 = Utility.ParseDouble(input);
            Assert.That.EqualsExactly(r, r1);
        }

        [DataTestMethod]
        [DataRow(nan,  "nan")]
        [DataRow(-inf, "-inf")]
        [DataRow(-3.5, "-4")]
        [DataRow(-2.5, "-2")]
        [DataRow(-1.5, "-2")]
        [DataRow(-0.5, "-0")]
        [DataRow(-0.0, "-0")]
        [DataRow(+0.0, "0")]
        [DataRow(+0.5, "0")]
        [DataRow(+1.5, "2")]
        [DataRow(+2.5, "2")]
        [DataRow(+3.5, "4")]
        [DataRow(+inf, "inf")]
        public void TestToFixedstring(double value, string r)
        {
            var r1 = value.ToFixedString(0);
            Assert.AreEqual(r, r1);
        }
    }
}
