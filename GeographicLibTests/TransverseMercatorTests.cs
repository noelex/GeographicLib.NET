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
    public class TransverseMercatorTests
    {
        private const double tolerance = 1e-8;

        [DataTestMethod]
        [DynamicData("Data", typeof(TransverseMercatorTestData))]
        public void TestForward(double lat, double lon, double x, double y, double gamma, double k)
        {
            var (_x, _y) = TransverseMercator.UTM.Forward(0, lat, lon, out var _gamma, out var _k);
            Assert.AreEqual(x, _x, tolerance);
            Assert.AreEqual(y, _y, tolerance);
            Assert.AreEqual(gamma, _gamma, tolerance);
            Assert.AreEqual(k, _k, tolerance);
        }

        //[DataTestMethod]
        //[DynamicData("Data", typeof(TransverseMercatorTestData))]
        //public void TestReverse(double lat, double lon, double x, double y, double gamma, double k)
        //{
        //    var (_lat, _lon) = TransverseMercator.UTM.Reverse(0, x, y, out var _gamma, out var _k);
        //    Assert.AreEqual(lat, _lat, tolerance);
        //    Assert.AreEqual(lon, _lon, tolerance);
        //    Assert.AreEqual(gamma, _gamma, 1e-8);
        //    Assert.AreEqual(k, _k, tolerance);
        //}
    }
}
