using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib.Tests
{
    [TestClass]
    public class RhumbTests
    {
        [DataTestMethod]
        [DataRow("40:38:23N 073:46:44W", "01:21:33N 103:59:22E", 103.582833003410954, 18523563.0423774272, 45921660958919.070, false)]
        [DataRow("33.3 44.4", "12.34 56.78", 151.483174375421527, 2641839.7239669343, 3408569853610.518, false)]
        [DataRow("40:38:23N 073:46:44W", "01:21:33N 103:59:22E", 103.582833003410954, 18523563.0423774309, 45921660958919.070, true)]
        [DataRow("33.3 44.4", "12.34 56.78", 151.483174375421527, 2641839.7239669338, 3408569853610.518, true)]
        public void TestInverse(string p1, string p2, double azi12, double s12, double S12, bool exact)
        {
            var c1=new GeoCoords(p1);
            var c2 = new GeoCoords(p2);

            new Rhumb(Ellipsoid.WGS84, exact)
                .Inverse(c1.Latitude, c1.Longitude, c2.Latitude, c2.Longitude, out var _s12, out var _azi12, out var _S12);

            Assert.AreEqual(azi12, _azi12, 1e-8);
            Assert.AreEqual(s12, _s12, 1e-8);
            Assert.AreEqual(S12, _S12, 0.001);
        }

        [DataTestMethod]
        [DataRow("40:38:23N 073:46:44W", "01:21:33N 103:59:22E", 103.582833003410954, 18523563.0423774272, 45921660958919.102, false)]
        [DataRow("33.3 44.4", "12.34 56.78", 151.483174375421527, 2641839.7239669343, 3408569853610.518, false)]
        [DataRow("40:38:23N 073:46:44W", "01:21:33N 103:59:22E", 103.582833003410954, 18523563.0423774309, 45921660958919.102, true)]
        [DataRow("33.3 44.4", "12.34 56.78", 151.483174375421527, 2641839.7239669338, 3408569853610.517, true)]
        public void TestDirect(string p1, string p2, double azi12, double s12, double S12, bool exact)
        {
            var c1 = new GeoCoords(p1);
            var c2 = new GeoCoords(p2);

            new Rhumb(Ellipsoid.WGS84, exact)
                .Direct(c1.Latitude, c1.Longitude, azi12, s12, out var _lat2, out var _lon2, out var _S12);

            Assert.AreEqual(c2.Latitude, _lat2, 1e-8);
            Assert.AreEqual(c2.Longitude, _lon2, 1e-8);

            Assert.AreEqual(S12, _S12, 0.01);
        }
    }
}
