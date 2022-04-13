﻿using Microsoft.VisualStudio.TestTools.UnitTesting;
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
        public static IEnumerable<object[]> TestCases =>
            new[]
            {
                new object[]{35.60777, -139.44815, 111.098748429560326,
     -11.17491, -69.95921, 129.289270889708762,
     8935244.5604818305, 80.50729714281974, 6273170.2055303837,
     0.16606318447386067, 0.16479116945612937, 12841384694976.432},
    new object[]{55.52454, 106.05087, 22.020059880982801,
     77.03196, 197.18234, 109.112041110671519,
     4105086.1713924406, 36.892740690445894, 3828869.3344387607,
     0.80076349608092607, 0.80101006984201008, 61674961290615.615},
    new object[]{-21.97856, 142.59065, -32.44456876433189,
     41.84138, 98.56635, -41.84359951440466,
     8394328.894657671, 75.62930491011522, 6161154.5773110616,
     0.24816339233950381, 0.24930251203627892, -6637997720646.717},
    new object[]{-66.99028, 112.2363, 173.73491240878403,
     -12.70631, 285.90344, 2.512956620913668,
     11150344.2312080241, 100.278634181155759, 6289939.5670446687,
     -0.17199490274700385, -0.17722569526345708, -121287239862139.744},
    new object[]{-17.42761, 173.34268, -159.033557661192928,
     -15.84784, 5.93557, -20.787484651536988,
     16076603.1631180673, 144.640108810286253, 3732902.1583877189,
     -0.81273638700070476, -0.81299800519154474, 97825992354058.708},
    new object[]{32.84994, 48.28919, 150.492927788121982,
     -56.28556, 202.29132, 48.113449399816759,
     16727068.9438164461, 150.565799985466607, 3147838.1910180939,
     -0.87334918086923126, -0.86505036767110637, -72445258525585.010},
    new object[]{6.96833, 52.74123, 92.581585386317712,
     -7.39675, 206.17291, 90.721692165923907,
     17102477.2496958388, 154.147366239113561, 2772035.6169917581,
     -0.89991282520302447, -0.89986892177110739, -1311796973197.995},
    new object[]{-50.56724, -16.30485, -105.439679907590164,
     -33.56571, -94.97412, -47.348547835650331,
     6455670.5118668696, 58.083719495371259, 5409150.7979815838,
     0.53053508035997263, 0.52988722644436602, 41071447902810.047},
    new object[]{-58.93002, -8.90775, 140.965397902500679,
     -8.91104, 133.13503, 19.255429433416599,
     11756066.0219864627, 105.755691241406877, 6151101.2270708536,
     -0.26548622269867183, -0.27068483874510741, -86143460552774.735},
    new object[]{-68.82867, -74.28391, 93.774347763114881,
     -50.63005, -8.36685, 34.65564085411343,
     3956936.926063544, 35.572254987389284, 3708890.9544062657,
     0.81443963736383502, 0.81420859815358342, -41845309450093.787},
    new object[]{-10.62672, -32.0898, -86.426713286747751,
     5.883, -134.31681, -80.473780971034875,
     11470869.3864563009, 103.387395634504061, 6184411.6622659713,
     -0.23138683500430237, -0.23155097622286792, 4198803992123.548},
    new object[]{-21.76221, 166.90563, 29.319421206936428,
     48.72884, 213.97627, 43.508671946410168,
     9098627.3986554915, 81.963476716121964, 6299240.9166992283,
     0.13965943368590333, 0.14152969707656796, 10024709850277.476},
    new object[]{-19.79938, -174.47484, 71.167275780171533,
     -11.99349, -154.35109, 65.589099775199228,
     2319004.8601169389, 20.896611684802389, 2267960.8703918325,
     0.93427001867125849, 0.93424887135032789, -3935477535005.785},
    new object[]{-11.95887, -116.94513, 92.712619830452549,
     4.57352, 7.16501, 78.64960934409585,
     13834722.5801401374, 124.688684161089762, 5228093.177931598,
     -0.56879356755666463, -0.56918731952397221, -9919582785894.853},
    new object[]{-87.85331, 85.66836, -65.120313040242748,
     66.48646, 16.09921, -4.888658719272296,
     17286615.3147144645, 155.58592449699137, 2635887.4729110181,
     -0.90697975771398578, -0.91095608883042767, 42667211366919.534},
    new object[]{1.74708, 128.32011, -101.584843631173858,
     -11.16617, 11.87109, -86.325793296437476,
     12942901.1241347408, 116.650512484301857, 5682744.8413270572,
     -0.44857868222697644, -0.44824490340007729, 10763055294345.653},
    new object[]{-25.72959, -144.90758, -153.647468693117198,
     -57.70581, -269.17879, -48.343983158876487,
     9413446.7452453107, 84.664533838404295, 6356176.6898881281,
     0.09492245755254703, 0.09737058264766572, 74515122850712.444},
    new object[]{-41.22777, 122.32875, 14.285113402275739,
     -7.57291, 130.37946, 10.805303085187369,
     3812686.035106021, 34.34330804743883, 3588703.8812128856,
     0.82605222593217889, 0.82572158200920196, -2456961531057.857},
    new object[]{11.01307, 138.25278, 79.43682622782374,
     6.62726, 247.05981, 103.708090215522657,
     11911190.819018408, 107.341669954114577, 6070904.722786735,
     -0.29767608923657404, -0.29785143390252321, 17121631423099.696},
    new object[]{-29.47124, 95.14681, -163.779130441688382,
     -27.46601, -69.15955, -15.909335945554969,
     13487015.8381145492, 121.294026715742277, 5481428.9945736388,
     -0.51527225545373252, -0.51556587964721788, 104679964020340.318 }
    };

        [DataTestMethod]
        [DynamicData(nameof(TestCases), typeof(GeodesicTests))]
        public void TestInverse_Incorporated(
            double lat1, double lon1, double azi1, double lat2, double lon2, double azi2,
            double s12, double a12, double m12, double M12, double M21, double S12)
        {
            var result = Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2, GeodesicFlags.All | GeodesicFlags.LongUnroll);

            Assert.AreEqual(azi1, result.Azimuth1, 1e-13);
            Assert.AreEqual(azi2, result.Azimuth2, 1e-13);
            Assert.AreEqual(s12, result.Distance, 1e-8);
            Assert.AreEqual(a12, result.ArcLength, 1e-13);
            Assert.AreEqual(m12, result.ReducedLength, 1e-8);
            Assert.AreEqual(M12, result.GeodesicScale12, 1e-15);
            Assert.AreEqual(M21, result.GeodesicScale21, 1e-15);
            Assert.AreEqual(S12, result.Area, 0.1);
        }

        [DataTestMethod]
        [DynamicData(nameof(TestCases), typeof(GeodesicTests))]
        public void TestDirect_Incorporated(
            double lat1, double lon1, double azi1, double lat2, double lon2, double azi2,
            double s12, double a12, double m12, double M12, double M21, double S12)
        {
            var result = Geodesic.WGS84.Direct(lat1, lon1, azi1, s12, GeodesicFlags.All | GeodesicFlags.LongUnroll);

            Assert.AreEqual(lat2, result.Latitude, 1e-13);
            Assert.AreEqual(lon2, result.Longitude, 1e-13);
            Assert.AreEqual(azi2, result.Azimuth2, 1e-13);
            Assert.AreEqual(a12, result.ArcLength, 1e-13);
            Assert.AreEqual(m12, result.ReducedLength, 1e-8);
            Assert.AreEqual(M12, result.GeodesicScale12, 1e-15);
            Assert.AreEqual(M21, result.GeodesicScale21, 1e-15);
            Assert.AreEqual(S12, result.Area, 0.1);
        }

        [DataTestMethod]
        [DynamicData(nameof(TestCases), typeof(GeodesicTests))]
        public void TestArcDirect_Incorporated(
            double lat1, double lon1, double azi1, double lat2, double lon2, double azi2,
            double s12, double a12, double m12, double M12, double M21, double S12)
        {
            var result = Geodesic.WGS84.ArcDirect(lat1, lon1, azi1, a12, GeodesicFlags.All | GeodesicFlags.LongUnroll);

            Assert.AreEqual(lat2, result.Latitude, 1e-13);
            Assert.AreEqual(lon2, result.Longitude, 1e-13);
            Assert.AreEqual(azi2, result.Azimuth2, 1e-13);
            Assert.AreEqual(a12, result.ArcLength, 1e-13);
            Assert.AreEqual(m12, result.ReducedLength, 1e-8);
            Assert.AreEqual(M12, result.GeodesicScale12, 1e-15);
            Assert.AreEqual(M21, result.GeodesicScale21, 1e-15);
            Assert.AreEqual(S12, result.Area, 0.1);
        }

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

        [DataTestMethod]
        [DataRow(+0d, -0d, 180d)]
        [DataRow(-0d, +0d, 0d)]
        public void TestInverse_AzimuthOfGeodesicLineOnEquatorDeterminedBySignsOfLatitude(double lat1, double lat2, double azi)
        {
            Geodesic.WGS84.Inverse(lat1, 0, lat2, 0, out var azi1, out var azi2);
            Assert.That.EqualsExactly(azi, azi1);
            Assert.That.EqualsExactly(azi, azi2);
        }

        [DataTestMethod]
        [DataRow(+0d, -0d, 180d)]
        [DataRow(-0d, +0d, 0d)]
        public void TestInverse_AzimuthOfGeodesicLineOnEquatorDeterminedBySignsOfLatitude_Exact(double lat1, double lat2, double azi)
        {
            GeodesicExact.WGS84.Inverse(lat1, 0, lat2, 0, out var azi1, out var azi2);
            Assert.That.EqualsExactly(azi, azi1);
            Assert.That.EqualsExactly(azi, azi2);
        }

        [DataTestMethod]
        [DataRow(+0d, +0d, 56d, 124d)]
        [DataRow(-0d, -0d, 124d, 56d)]
        public void TestInverse_DirectionOfAntipodalEquatorialSolution(double lat1, double lat2, double azi1, double azi2)
        {
            Geodesic.WGS84.Inverse(lat1, 0, lat2, 179.5, out var azi1a, out var azi2a);
            Assert.AreEqual(azi1, azi1a, 1);
            Assert.AreEqual(azi2, azi2a, 1);
        }

        [DataTestMethod]
        [DataRow(+0d, +0d, 56d, 124d)]
        [DataRow(-0d, -0d, 124d, 56d)]
        public void TestInverse_DirectionOfNearlyAntipodalEquatorialSolution_Exact(double lat1, double lat2, double azi1, double azi2)
        {
            GeodesicExact.WGS84.Inverse(lat1, 0, lat2, 179.5, out var azi1a, out var azi2a);
            Assert.AreEqual(azi1, azi1a, 1);
            Assert.AreEqual(azi2, azi2a, 1);
        }

        [DataTestMethod]
        [DataRow(+0d, +0d, +180d, +0d, +180d)]
        [DataRow(-0d, -0d, +180d, +180d, +0d)]
        [DataRow(+0d, +0d, -180d, -0d, -180d)]
        [DataRow(-0d, -0d, -180d, -180d, -0d)]
        public void TestInverse_DirectionOfExactAntipodalEquatorialSolution(double lat1, double lat2, double lon2, double azi1, double azi2)
        {
            Geodesic.WGS84.Inverse(lat1, 0, lat2, lon2, out var azi1a, out var azi2a);
            Assert.That.EqualsExactly(azi1, azi1a);
            Assert.That.EqualsExactly(azi2, azi2a);
        }

        [DataTestMethod]
        [DataRow(+0d, +0d, +180d, +0d, +180d)]
        [DataRow(-0d, -0d, +180d, +180d, +0d)]
        [DataRow(+0d, +0d, -180d, -0d, -180d)]
        [DataRow(-0d, -0d, -180d, -180d, -0d)]
        public void TestInverse_DirectionOfExactAntipodalEquatorialSolution_Exact(double lat1, double lat2, double lon2, double azi1, double azi2)
        {
            GeodesicExact.WGS84.Inverse(lat1, 0, lat2, lon2, out var azi1a, out var azi2a);
            Assert.That.EqualsExactly(azi1, azi1a);
            Assert.That.EqualsExactly(azi2, azi2a);
        }

        [DataTestMethod]
        [DataRow(+180d, +90d)]
        [DataRow(-180d, -90d)]
        public void TestInversion_AntipodalPointsOnTheEquatorWithProlateEllipsoid(double lon2, double azi)
        {
            new Geodesic(6.4e6, -1/300d).Inverse(0, 0, 0, lon2, out var azi1, out var azi2);
            Assert.That.EqualsExactly(azi, azi1);
            Assert.That.EqualsExactly(azi, azi2);
        }

        [DataTestMethod]
        [DataRow(+180d, +90d)]
        [DataRow(-180d, -90d)]
        public void TestInversion_AntipodalPointsOnTheEquatorWithProlateEllipsoid_Exact(double lon2, double azi)
        {
            new GeodesicExact(6.4e6, -1 / 300d).Inverse(0, 0, 0, lon2, out var azi1, out var azi2);
            Assert.That.EqualsExactly(azi, azi1);
            Assert.That.EqualsExactly(azi, azi2);
        }

        [DataTestMethod]
        [DataRow(+0d, +180d, +180d)]
        [DataRow(-0d, -180d, -180d)]
        [DataRow(+180d, +180d, +0d)]
        [DataRow(-180d, -180d, -0d)]
        public void TestDirect_Azimuth_Equals_PlusMinus0_And_PlusMinus180(double azi1, double lon2, double azi2)
        {
            Geodesic.WGS84.GenDirect(0, 0, azi1, false, 15e6, 
                GeodesicFlags.Longitude | GeodesicFlags.Azimuth | GeodesicFlags.LongUnroll,
                out _, out var lon2a, out var azi2a,
                out _, out _, out _, out _, out _);

            Assert.That.EqualsExactly(lon2, lon2a);
            Assert.That.EqualsExactly(azi2, azi2a);
        }

        [DataTestMethod]
        [DataRow(+0d, +180d, +180d)]
        [DataRow(-0d, -180d, -180d)]
        [DataRow(+180d, +180d, +0d)]
        [DataRow(-180d, -180d, -0d)]
        public void TestDirect_Azimuth_Equals_PlusMinus0_And_PlusMinus180_Exact(double azi1, double lon2, double azi2)
        {
            GeodesicExact.WGS84.GenDirect(0, 0, azi1, false, 15e6,
                GeodesicFlags.Longitude | GeodesicFlags.Azimuth | GeodesicFlags.LongUnroll,
                out _, out var lon2a, out var azi2a,
                out _, out _, out _, out _, out _);

            Assert.That.EqualsExactly(lon2, lon2a);
            Assert.That.EqualsExactly(azi2, azi2a);
        }
    }
}
