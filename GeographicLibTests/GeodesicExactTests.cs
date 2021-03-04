using GeographicLib;
using GeographicLib.Tests;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib.Tests
{
    [TestClass]
    public class GeodesicExactTests
    {
        [DataTestMethod]
        [DynamicData("Subsample", typeof(GeodesicTestData))]
        public void TestDirectFromPointOne(
            double lat1, double lon1, double azi1, double lat2, double lon2, double azi2, double s12, double a12, double m12, double S12)
        {
            var geodesic = new GeodesicExact(Ellipsoid.WGS84);
            var ra12 = geodesic.Direct(lat1, lon1, azi1, s12, out var rlat2, out var rlon2, out var razi2, out var rm12);

            Assert.AreEqual(lat2, rlat2, GeodesicTestData.ToleranceExact);
            Assert.AreEqual(lon2, rlon2, GeodesicTestData.ToleranceExact);
            Assert.AreEqual(azi2, razi2, GeodesicTestData.ToleranceExact);
            Assert.AreEqual(m12, rm12, GeodesicTestData.ToleranceExact);
            Assert.AreEqual(a12, ra12, GeodesicTestData.ToleranceExact);
        }

        [DataTestMethod]
        [DynamicData("Subsample", typeof(GeodesicTestData))]
        public void TestDirectFromPointTwo(
            double lat1, double lon1, double azi1, double lat2, double lon2, double azi2, double s12, double a12, double m12, double S12)
        {
            var geodesic = new GeodesicExact(Ellipsoid.WGS84);
            var ra12 = geodesic.Direct(lat2, lon2, azi2, -s12, out var rlat1, out var rlon1, out var razi1, out var rm12);

            Assert.AreEqual(lat1, rlat1, GeodesicTestData.ToleranceExact);
            Assert.AreEqual(lon1, rlon1, GeodesicTestData.ToleranceExact);
            Assert.AreEqual(azi1, razi1, GeodesicTestData.ToleranceExact);
            Assert.AreEqual(m12, -rm12, GeodesicTestData.ToleranceExact);
            Assert.AreEqual(a12, -ra12, GeodesicTestData.ToleranceExact);
        }

        [DataTestMethod]
        [DynamicData("Subsample", typeof(GeodesicTestData))]
        public void TestInverse(
            double lat1, double lon1, double azi1, double lat2, double lon2, double azi2, double s12, double a12, double m12, double S12)
        {
            var geodesic = new GeodesicExact(Ellipsoid.WGS84);
            var ra12 = geodesic.Inverse(lat1, lon1, lat2, lon2, out var rs12, out var razi1, out var razi2, out var rm12);

            // Accuracy decreased when running between vertices (azi1 = azi2 = 90°), verified with GeodSolve.
            var vertices = Math.Abs(90 - azi1) < 0.1 && Math.Abs(90 - azi2) < 0.1;

            Assert.AreEqual(azi1, razi1, vertices ? 0.001 : GeodesicTestData.ToleranceExact);
            Assert.AreEqual(azi2, razi2, vertices ? 0.001 : GeodesicTestData.ToleranceExact);

            Assert.AreEqual(s12, rs12, GeodesicTestData.ToleranceExact);
            Assert.AreEqual(m12, rm12, vertices ? 1e-5 : GeodesicTestData.ToleranceExact);
            Assert.AreEqual(a12, ra12, GeodesicTestData.ToleranceExact);
        }
    }
}
