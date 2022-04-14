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

        [DataTestMethod]
        [DataRow(+0d, +180d, +180d)]
        [DataRow(-0d, -180d, -180d)]
        [DataRow(+180d, +180d, +0d)]
        [DataRow(-180d, -180d, -0d)]
        public void TestDirect_Azimuth_Equals_PlusMinus0_And_PlusMinus180(double azi1, double lon2, double azi2)
        {
            GeodesicExact.WGS84.GenDirect(0, 0, azi1, false, 15e6,
                GeodesicFlags.Longitude | GeodesicFlags.Azimuth | GeodesicFlags.LongUnroll,
                out _, out var lon2a, out var azi2a,
                out _, out _, out _, out _, out _);

            Assert.That.EqualsExactly(lon2, lon2a);
            Assert.That.EqualsExactly(azi2, azi2a);
        }

        [DataTestMethod]
        [DataRow(+180d, +90d)]
        [DataRow(-180d, -90d)]
        public void TestInversion_AntipodalPointsOnTheEquatorWithProlateEllipsoid(double lon2, double azi)
        {
            new GeodesicExact(6.4e6, -1 / 300d).Inverse(0, 0, 0, lon2, out var azi1, out var azi2);
            Assert.That.EqualsExactly(azi, azi1);
            Assert.That.EqualsExactly(azi, azi2);
        }

        [DataTestMethod]
        [DataRow(+0d, +0d, +180d, +0d, +180d)]
        [DataRow(-0d, -0d, +180d, +180d, +0d)]
        [DataRow(+0d, +0d, -180d, -0d, -180d)]
        [DataRow(-0d, -0d, -180d, -180d, -0d)]
        public void TestInverse_DirectionOfExactAntipodalEquatorialSolution(double lat1, double lat2, double lon2, double azi1, double azi2)
        {
            GeodesicExact.WGS84.Inverse(lat1, 0, lat2, lon2, out var azi1a, out var azi2a);
            Assert.That.EqualsExactly(azi1, azi1a);
            Assert.That.EqualsExactly(azi2, azi2a);
        }

        [DataTestMethod]
        [DataRow(+0d, +0d, 56d, 124d)]
        [DataRow(-0d, -0d, 124d, 56d)]
        public void TestInverse_DirectionOfNearlyAntipodalEquatorialSolution(double lat1, double lat2, double azi1, double azi2)
        {
            GeodesicExact.WGS84.Inverse(lat1, 0, lat2, 179.5, out var azi1a, out var azi2a);
            Assert.AreEqual(azi1, azi1a, 1);
            Assert.AreEqual(azi2, azi2a, 1);
        }

        [DataTestMethod]
        [DataRow(+0d, -0d, 180d)]
        [DataRow(-0d, +0d, 0d)]
        public void TestInverse_AzimuthOfGeodesicLineOnEquatorDeterminedBySignsOfLatitude(double lat1, double lat2, double azi)
        {
            GeodesicExact.WGS84.Inverse(lat1, 0, lat2, 0, out var azi1, out var azi2);
            Assert.That.EqualsExactly(azi, azi1);
            Assert.That.EqualsExactly(azi, azi2);
        }

        [DataTestMethod]
        [DynamicData("GeodTests", typeof(GeodesicTestData))]
        public void TestInverse_RandomGeodesicProblems(
            double lat1, double lon1, double azi1,
            double lat2, double lon2, double azi2,
            double s12, double a12, double m12,
            double M12, double M21, double S12)
        {
            var f = 2;
            var a12a = GeodesicExact.WGS84.GenInverse(lat1, lon1, lat2, lon2, GeodesicFlags.All,
                out var s12a, out var azi1a, out var azi2a,
                out var m12a, out var M12a, out var M21a, out var S12a);

            Assert.AreEqual(azi1, azi1a, 1e-13 * f);
            Assert.AreEqual(azi2, azi2a, 1e-13 * f);
            Assert.AreEqual(s12, s12a, 1e-8 * f);
            Assert.AreEqual(a12, a12a, 1e-13 * f);
            Assert.AreEqual(m12, m12a, 1e-8 * f);
            Assert.AreEqual(M12, M12a, 1e-15 * f);
            Assert.AreEqual(M21, M21a, 1e-15 * f);
            Assert.AreEqual(S12, S12a, 0.1 * f);
        }

        [DataTestMethod]
        [DynamicData("GeodTests", typeof(GeodesicTestData))]
        public void TestDirect_RandomGeodesicProblems(
            double lat1, double lon1, double azi1,
            double lat2, double lon2, double azi2,
            double s12, double a12, double m12,
            double M12, double M21, double S12)
        {
            var f = 2;
            var a12a = GeodesicExact.WGS84.GenDirect(lat1, lon1, azi1, false, s12, GeodesicFlags.All | GeodesicFlags.LongUnroll,
                out var lat2a, out var lon2a, out var azi2a, out var s12a,
                out var m12a, out var M12a, out var M21a, out var S12a);

            Assert.AreEqual(lat2, lat2a, 1e-13 * f);
            Assert.AreEqual(lon2, lon2a, 1e-13 * f);
            Assert.AreEqual(azi2, azi2a, 1e-13 * f);
            Assert.AreEqual(s12, s12a, 0 * f);
            Assert.AreEqual(a12, a12a, 1e-13 * f);
            Assert.AreEqual(m12, m12a, 1e-8 * f);
            Assert.AreEqual(M12, M12a, 1e-15 * f);
            Assert.AreEqual(M21, M21a, 1e-15 * f);
            Assert.AreEqual(S12, S12a, 0.1 * f);
        }

        [DataTestMethod]
        [DynamicData("GeodTests", typeof(GeodesicTestData))]
        public void TestArcDirect_RandomGeodesicProblems(
            double lat1, double lon1, double azi1,
            double lat2, double lon2, double azi2,
            double s12, double a12, double m12,
            double M12, double M21, double S12)
        {
            var f = 2;
            var a12a = GeodesicExact.WGS84.GenDirect(lat1, lon1, azi1, true, a12, GeodesicFlags.All | GeodesicFlags.LongUnroll,
                out var lat2a, out var lon2a, out var azi2a, out var s12a,
                out var m12a, out var M12a, out var M21a, out var S12a);

            Assert.AreEqual(lat2, lat2a, 1e-13 * f);
            Assert.AreEqual(lon2, lon2a, 1e-13 * f);
            Assert.AreEqual(azi2, azi2a, 1e-13 * f);
            Assert.AreEqual(s12, s12a, 1e-8 * f); // TODO: delta should be 0 here
            Assert.AreEqual(a12, a12a, 1e-13 * f);
            Assert.AreEqual(m12, m12a, 1e-8 * f);
            Assert.AreEqual(M12, M12a, 1e-15 * f);
            Assert.AreEqual(M21, M21a, 1e-15 * f);
            Assert.AreEqual(S12, S12a, 0.1 * f);
        }
    }
}
