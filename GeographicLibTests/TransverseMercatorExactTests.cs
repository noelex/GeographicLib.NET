using GeographicLib.Projections;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib.Tests
{
    [TestClass]
    public class TransverseMercatorExactTests
    {
        [DataTestMethod]
        [DynamicData("ExactData", typeof(TransverseMercatorTestData))]
        public void TestForward(double lat, double lon, double x, double y, double gamma, double k)
        {
            var (_x, _y) = TransverseMercatorExact.UTM.Forward(0, lat, lon, out var _gamma, out var _k);
            Assert.AreEqual(x, _x, TransverseMercatorTestData.Tolerance);
            Assert.AreEqual(y, _y, TransverseMercatorTestData.Tolerance);
            Assert.AreEqual(gamma, _gamma, TransverseMercatorTestData.Tolerance);
            Assert.AreEqual(k, _k, TransverseMercatorTestData.Tolerance);
        }

        //[DataTestMethod]
        //[DynamicData("ExactData", typeof(TransverseMercatorTestData))]
        //public void TestReverse(double lat, double lon, double x, double y, double gamma, double k)
        //{
        //    var (_lat, _lon) = TransverseMercatorExact.UTM.Reverse(0, x, y, out var _gamma, out var _k);
        //    Assert.AreEqual(lat, _lat, TransverseMercatorTestData.ToleranceReverse);
        //    Assert.AreEqual(lon, _lon, TransverseMercatorTestData.ToleranceReverse);
        //    Assert.AreEqual(gamma, _gamma, TransverseMercatorTestData.ToleranceReverse);
        //    Assert.AreEqual(k, _k, TransverseMercatorTestData.ToleranceReverse);
        //}
    }
}
