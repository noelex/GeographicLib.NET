using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Text;
using GeographicLib;
using System.Linq;

namespace GeographicLib.Tests
{
    [TestClass]
    public class GeodesicTests
    {
        [DataTestMethod]
        [DynamicData("Subsample", typeof(GeodesicTestData))]
        public void TestDirectFromPointOne(
            double lat1, double lon1, double azi1, double lat2, double lon2, double azi2, double s12, double a12, double m12, double S12)
        {
            var geodesic = new Geodesic(Ellipsoid.WGS84);
            var ra12 = geodesic.Direct(lat1, lon1, azi1, s12, out var rlat2, out var rlon2, out var razi2, out var rm12);

            Assert.AreEqual(lat2, rlat2, GeodesicTestData.Tolerance);
            Assert.AreEqual(lon2, rlon2, GeodesicTestData.Tolerance);
            Assert.AreEqual(azi2, razi2, GeodesicTestData.Tolerance);
            Assert.AreEqual(m12, rm12, GeodesicTestData.Tolerance);
            Assert.AreEqual(a12, ra12, GeodesicTestData.Tolerance);
        }

        [DataTestMethod]
        [DynamicData("Subsample", typeof(GeodesicTestData))]
        public void TestDirectFromPointTwo(
            double lat1, double lon1, double azi1, double lat2, double lon2, double azi2, double s12, double a12, double m12, double S12)
        {
            var geodesic = new Geodesic(Ellipsoid.WGS84);
            var ra12 = geodesic.Direct(lat2, lon2, azi2, -s12, out var rlat1, out var rlon1, out var razi1, out var rm12);

            Assert.AreEqual(lat1, rlat1, GeodesicTestData.Tolerance);
            Assert.AreEqual(lon1, rlon1, GeodesicTestData.Tolerance);
            Assert.AreEqual(azi1, razi1, GeodesicTestData.Tolerance);
            Assert.AreEqual(m12, -rm12, GeodesicTestData.Tolerance);
            Assert.AreEqual(a12, -ra12, GeodesicTestData.Tolerance);
        }

        [DataTestMethod]
        [DynamicData("Subsample", typeof(GeodesicTestData))]
        public void TestInverse(
            double lat1, double lon1, double azi1, double lat2, double lon2, double azi2, double s12, double a12, double m12, double S12)
        {
            var geodesic = new Geodesic(Ellipsoid.WGS84);
            var ra12 = geodesic.Inverse(lat1, lon1, lat2, lon2, out var rs12, out var razi1, out var razi2, out var rm12);

            // Accuracy decreased when running between vertices (azi1 = azi2 = 90°), verified with GeodSolve.
            var vertices = Math.Abs(90 - azi1) < 0.1 && Math.Abs(90 - azi2) < 0.1;

            Assert.AreEqual(azi1, razi1, vertices ? 0.001 : GeodesicTestData.Tolerance);
            Assert.AreEqual(azi2, razi2, vertices ? 0.001 : GeodesicTestData.Tolerance);

            Assert.AreEqual(s12, rs12, GeodesicTestData.Tolerance);
            Assert.AreEqual(m12, rm12, vertices ? 1e-5 : GeodesicTestData.Tolerance);
            Assert.AreEqual(a12, ra12, GeodesicTestData.Tolerance);
        }
    }
}
