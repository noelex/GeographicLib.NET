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
        [DataRow(nan, 0, "nan")]
        [DataRow(-inf, 0, "-inf")]
        [DataRow(+inf, 0, "inf")]
        [DataRow(-3.5, 0, "-4")]
        [DataRow(-2.5, 0, "-2")]
        [DataRow(-1.5, 0, "-2")]
        [DataRow(-0.5, 0, "-0")]
        [DataRow(-0.0, 0, "-0")]
        [DataRow(+0.0, 0, "0")]
        [DataRow(+0.5, 0, "0")]
        [DataRow(+1.5, 0, "2")]
        [DataRow(+2.5, 0, "2")]
        [DataRow(+3.5, 0, "4")]
        [DataRow(-1.75, 1, "-1.8")]
        [DataRow(-1.25, 1, "-1.2")]
        [DataRow(-0.75, 1, "-0.8")]
        [DataRow(-0.25, 1, "-0.2")]
        [DataRow(-0.0, 1, "-0.0")]
        [DataRow(+0.0, 1, "0.0")]
        [DataRow(+0.25, 1, "0.2")]
        [DataRow(+0.75, 1, "0.8")]
        [DataRow(+1.25, 1, "1.2")]
        [DataRow(+1.75, 1, "1.8")]
        public void TestToFixedstring(double value,int prec, string r)
        {
            var r1 = value.ToFixedString(prec);
            Assert.AreEqual(r, r1);
        }
    }
}
