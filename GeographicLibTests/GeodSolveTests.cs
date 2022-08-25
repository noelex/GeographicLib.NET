using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using static System.Math;
#if NETCOREAPP2_1
using static GeographicLib.MathEx;
#endif

namespace GeographicLib.Tests
{
    [TestClass]
    public class GeodSolveTests
    {
        [TestMethod]
        public void GeodSolve0()
        {
            Geodesic.WGS84.Inverse(40.6, -73.8, 49.01666667, 2.55,
                out var s12, out var azi1, out var azi2);

            Assert.AreEqual(53.47022, azi1, 0.5e-5);
            Assert.AreEqual(111.59367, azi2, 0.5e-5);
            Assert.AreEqual(5853226, s12, 0.5);
        }

        [TestMethod]
        public void GeodSolve1()
        {
            Geodesic.WGS84.Direct(40.63972222, -73.77888889, 53.5, 5850e3,
                out var lat2, out var lon2, out var azi2);

            Assert.AreEqual(49.01467, lat2, 0.5e-5);
            Assert.AreEqual(2.56106, lon2, 0.5e-5);
            Assert.AreEqual(111.62947, azi2, 0.5);
        }

        [TestMethod]
        public void GeodSolve2_CheckFixForAntipodalProlateBug()
        {
            var geod = new Geodesic(6.4e6, -1 / 150.0);

            geod.Inverse(0.07476, 0, -0.07476, 180,
                out var s12, out var azi1, out var azi2);

            Assert.AreEqual(90.00078, azi1, 0.5e-5);
            Assert.AreEqual(90.00078, azi2, 0.5e-5);
            Assert.AreEqual(20106193, s12, 0.5);
        }

        [TestMethod]
        public void GeodSolve4_CheckFixForShortLineBug()
        {
            Geodesic.WGS84.Inverse(36.493349428792, 0, 36.49334942879201, 0.0000008, out var s12);

            Assert.AreEqual(0.072, s12, 0.5e-3);
        }

        [TestMethod]
        public void GeodSolve5_CheckFixForPoint2IsPoleBug()
        {
            Geodesic.WGS84.Direct(0.01777745589997, 30, 0, 10e6, out var lat2, out var lon2, out var azi2);

            Assert.AreEqual(90, lat2, 0.5e-5);
            if (lon2 < 0)
            {
                Assert.AreEqual(-150, lon2, 0.5e-5);
                Assert.AreEqual(180, Abs(azi2), 0.5e-5);
            }
            else
            {
                Assert.AreEqual(30, lon2, 0.5e-5);
                Assert.AreEqual(0, azi2, 0.5e-5);
            }
        }

        [TestMethod]
        public void GeodSolve6_CheckFixForVolatile_sbt12a_Bug()
        {
            Geodesic.WGS84.Inverse(88.202499451857, 0, -88.202499451857, 179.981022032992859592, out var s12);
            Assert.AreEqual(20003898.214, s12, 0.5e-3);

            Geodesic.WGS84.Inverse(89.262080389218, 0, -89.262080389218, 179.992207982775375662, out s12);
            Assert.AreEqual(20003925.854, s12, 0.5e-3);

            Geodesic.WGS84.Inverse(89.333123580033, 0, -89.333123580032997687, 179.99295812360148422, out s12);
            Assert.AreEqual(20003926.881, s12, 0.5e-3);
        }

        [TestMethod]
        public void GeodSolve9_CheckFixForVolatile_x_Bug()
        {
            Geodesic.WGS84.Inverse(56.320923501171, 0, -56.320923501171, 179.664747671772880215, out var s12);
            Assert.AreEqual(19993558.287, s12, 0.5e-3);
        }

        [TestMethod]
        public void GeodSolve10_CheckFixForAdjust_tol1_Bug()
        {
            Geodesic.WGS84.Inverse(52.784459512564, 0, -52.784459512563990912, 179.634407464943777557, out var s12);
            Assert.AreEqual(19991596.095, s12, 0.5e-3);
        }

        [TestMethod]
        public void GeodSolve11_CheckFixFor_bet2_Equals_Negative_bet1_BugBug()
        {
            Geodesic.WGS84.Inverse(48.522876735459, 0, -48.52287673545898293, 179.599720456223079643, out var s12);
            Assert.AreEqual(19989144.774, s12, 0.5e-3);
        }

        [TestMethod]
        public void GeodSolve12_CheckFixForInverseGeodesicsOnExtremeProlateOblateEllipsoids()
        {
            var geod = new Geodesic(89.8, -1.83);

            geod.Inverse(0, 0, -10, 160,
                out var s12, out var azi1, out var azi2);

            Assert.AreEqual(120.27, azi1, 1e-2);
            Assert.AreEqual(105.15, azi2, 1e-1);
            Assert.AreEqual(266.7, s12, 1e-1);
        }

        [TestMethod]
        public void GeodSolve14_CheckFixForInverseIgnoring_lon12_Equals_NaN()
        {
            Geodesic.WGS84.Inverse(0, 0, 1, double.NaN,
                out var s12, out var azi1, out var azi2);

            Assert.IsTrue(double.IsNaN(azi1));
            Assert.IsTrue(double.IsNaN(azi2));
            Assert.IsTrue(double.IsNaN(s12));
        }

        [TestMethod]
        public void GeodSolve15_CheckImplementationOfEAtahhEFor_e2_LessThan_Zero()
        {
            var geod = new Geodesic(6.4e6, -1 / 150.0);
            geod.GenDirect(1, 2, 3, false, 4, GeodesicFlags.Area,
                out _, out _, out _, out _, out _, out _, out _, out var S12);

            Assert.AreEqual(23700, S12, 0.5);
        }

        [TestMethod]
        public void GeodSolve17_CheckFixForLongUnrollBug()
        {
            Geodesic.WGS84.GenDirect(40, -75, -10, false, 2e7,
                GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth | GeodesicFlags.LongUnroll,
                out var lat2, out var lon2, out var azi2, out _, out _, out _, out _, out _);
            Assert.AreEqual(-39, lat2, 1);
            Assert.AreEqual(-254, lon2, 1);
            Assert.AreEqual(-170, azi2, 1);

            var line = Geodesic.WGS84.Line(40, -75, -10);
            line.GenPosition(
                false, 2e7,
                GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth | GeodesicFlags.LongUnroll,
                out lat2, out lon2, out azi2, out _, out _, out _, out _, out _);
            Assert.AreEqual(-39, lat2, 1);
            Assert.AreEqual(-254, lon2, 1);
            Assert.AreEqual(-170, azi2, 1);

            Geodesic.WGS84.Direct(40, -75, -10, 2e7, out lat2, out lon2, out azi2);
            Assert.AreEqual(-39, lat2, 1);
            Assert.AreEqual(105, lon2, 1);
            Assert.AreEqual(-170, azi2, 1);

            line.Position(2e7, out lat2, out lon2, out azi2);
            Assert.AreEqual(-39, lat2, 1);
            Assert.AreEqual(105, lon2, 1);
            Assert.AreEqual(-170, azi2, 1);
        }

        [TestMethod]
        public void GeodSolve26_Check0Div0ProblemWithAreaCalculationOnSphere()
        {
            var geod = new Geodesic(6.4e6, 0);
            geod.GenInverse(1, 2, 3, 4, GeodesicFlags.Area, out _, out _, out _, out _, out _, out _, out var S12);

            Assert.AreEqual(49911046115.0, S12, 0.5);
        }

        [TestMethod]
        public void GeodSolve28_CheckForBadPlacementOfAssignmentOf_r_a12_With_abs_f_GreaterThan0_01()
        {
            var geod = new Geodesic(6.4e6, 0.1);
            var a12 = geod.Direct(1, 2, 10, 5e6, out _, out _);

            Assert.AreEqual(48.55570690, a12, 0.5);
        }

        [TestMethod]
        public void GeodSolve33_CheckMax_Negative_Zero_Positive_Zero_Issues()
        {
            Geodesic.WGS84.Inverse(0, 0, 0, 179, out var s12, out var azi1, out var azi2);
            Assert.AreEqual(90.00000, azi1, 0.5e-5);
            Assert.AreEqual(90.00000, azi2, 0.5e-5);
            Assert.AreEqual(19926189, s12, 0.5);

            Geodesic.WGS84.Inverse(0, 0, 0, 179.5, out s12, out azi1, out azi2);
            Assert.AreEqual(55.96650, azi1, 0.5e-5);
            Assert.AreEqual(124.03350, azi2, 0.5e-5);
            Assert.AreEqual(19980862, s12, 0.5);

            Geodesic.WGS84.Inverse(0, 0, 0, 180, out s12, out azi1, out azi2);
            Assert.AreEqual(0.00000, azi1, 0.5e-5);
            Assert.AreEqual(180, Abs(azi2), 0.5e-5);
            Assert.AreEqual(20003931, s12, 0.5);

            Geodesic.WGS84.Inverse(0, 0, 1, 180, out s12, out azi1, out azi2);
            Assert.AreEqual(0.00000, azi1, 0.5e-5);
            Assert.AreEqual(180, Abs(azi2), 0.5e-5);
            Assert.AreEqual(19893357, s12, 0.5);

            var geod = new Geodesic(6.4e6, 0);
            geod.Inverse(0, 0, 0, 179, out s12, out azi1, out azi2);
            Assert.AreEqual(90.00000, azi1, 0.5e-5);
            Assert.AreEqual(90.00000, azi2, 0.5e-5);
            Assert.AreEqual(19994492, s12, 0.5);

            geod.Inverse(0, 0, 0, 180, out s12, out azi1, out azi2);
            Assert.AreEqual(0.00000, azi1, 0.5e-5);
            Assert.AreEqual(180, Abs(azi2), 0.5e-5);
            Assert.AreEqual(20106193, s12, 0.5);

            geod.Inverse(0, 0, 1, 180, out s12, out azi1, out azi2);
            Assert.AreEqual(0.00000, azi1, 0.5e-5);
            Assert.AreEqual(180, Abs(azi2), 0.5e-5);
            Assert.AreEqual(19994492, s12, 0.5);

            geod = new Geodesic(6.4e6, -1 / 300.0);
            geod.Inverse(0, 0, 0, 179, out s12, out azi1, out azi2);
            Assert.AreEqual(90.00000, azi1, 0.5e-5);
            Assert.AreEqual(90.00000, azi2, 0.5e-5);
            Assert.AreEqual(19994492, s12, 0.5);

            geod.Inverse(0, 0, 0, 180, out s12, out azi1, out azi2);
            Assert.AreEqual(90.00000, azi1, 0.5e-5);
            Assert.AreEqual(90.00000, azi2, 0.5e-5);
            Assert.AreEqual(20106193, s12, 0.5);

            geod.Inverse(0, 0, 0.5, 180, out s12, out azi1, out azi2);
            Assert.AreEqual(33.02493, azi1, 0.5e-5);
            Assert.AreEqual(146.97364, azi2, 0.5e-5);
            Assert.AreEqual(20082617, s12, 0.5);

            geod.Inverse(0, 0, 1, 180, out s12, out azi1, out azi2);
            Assert.AreEqual(0.00000, azi1, 0.5e-5);
            Assert.AreEqual(180, Abs(azi2), 0.5e-5);
            Assert.AreEqual(20027270, s12, 0.5);
        }

        [TestMethod]
        public void GeodSolve55_CheckFixForNaNPlusPointOnEquatorOrPoleNotReturningAllNaNsInInverse()
        {
            Geodesic.WGS84.Inverse(double.NaN, 0, 0, 90, out var s12, out var azi1, out var azi2);
            Assert.IsTrue(double.IsNaN(azi1));
            Assert.IsTrue(double.IsNaN(azi2));
            Assert.IsTrue(double.IsNaN(s12));

            Geodesic.WGS84.Inverse(double.NaN, 0, 90, 9, out s12, out azi1, out azi2);
            Assert.IsTrue(double.IsNaN(azi1));
            Assert.IsTrue(double.IsNaN(azi2));
            Assert.IsTrue(double.IsNaN(s12));
        }

        [TestMethod]
        public void GeodSolve59_CheckForPointsCloseWithLongitudesCloseTo180degApart()
        {
            Geodesic.WGS84.Inverse(5, 0.00000000000001, 10, 180, out var s12, out var azi1, out var azi2);
            Assert.AreEqual(0.000000000000035, azi1, 1.5e-14);
            Assert.AreEqual(179.99999999999996, azi2, 1.5e-14);
            Assert.AreEqual(18345191.174332713, s12, 5e-9);
        }

        [TestMethod]
        public void GeodSolve61_MakeSureSmallNegativeAzimuthsAreWestGoing()
        {
            Geodesic.WGS84.GenDirect(45, 0, -0.000000000000000003, false, 1e7,
                GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth | GeodesicFlags.LongUnroll,
                out var lat2, out var lon2, out var azi2, out _, out _, out _, out _, out _);
            Assert.AreEqual(45.30632, lat2, 0.5e-5);
            Assert.AreEqual(-180, lon2, 0.5e-5);
            Assert.AreEqual(180, Abs(azi2), 0.5e-5);

            var line = Geodesic.WGS84.InverseLine(45, 0, 80, -0.000000000000000003);
            line.GenPosition(false, 1e7,
                GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth | GeodesicFlags.LongUnroll,
                out lat2, out lon2, out azi2, out _, out _, out _, out _, out _);
            Assert.AreEqual(45.30632, lat2, 0.5e-5);
            Assert.AreEqual(-180, lon2, 0.5e-5);
            Assert.AreEqual(180, Abs(azi2), 0.5e-5);
        }

        [TestMethod]
        public void GeodSolve65_CheckForBugInEastGoingCheckInGeodesicLine()
        {
            var line = Geodesic.WGS84.InverseLine(30, -0.000000000000000001, -31, 180, GeodesicFlags.All);
            var a12 = line.GenPosition(false, 1e7, GeodesicFlags.All | GeodesicFlags.LongUnroll,
                out var lat2, out var lon2, out var azi2, out var s12, out var m12, out var M12, out var M21, out var S12);
            Assert.AreEqual(-60.23169, lat2, 0.5e-5);
            Assert.AreEqual(-0.00000, lon2, 0.5e-5);
            Assert.AreEqual(180.00000, Abs(azi2), 0.5e-5);
            Assert.AreEqual(10000000, s12, 0.5);
            Assert.AreEqual(90.06544, a12, 0.5e-5);
            Assert.AreEqual(6363636, m12, 0.5);
            Assert.AreEqual(-0.0012834, M12, 0.5e-7);
            Assert.AreEqual(0.0013749, M21, 0.5e-7);
            Assert.AreEqual(0, S12, 0.5);

            a12 = line.GenPosition(false, 2e7, GeodesicFlags.All | GeodesicFlags.LongUnroll,
                out lat2, out lon2, out azi2, out s12, out m12, out M12, out M21, out S12);
            Assert.AreEqual(-30.03547, lat2, 0.5e-5);
            Assert.AreEqual(-180.00000, lon2, 0.5e-5);
            Assert.AreEqual(-0.00000, azi2, 0.5e-5);
            Assert.AreEqual(20000000, s12, 0.5);
            Assert.AreEqual(179.96459, a12, 0.5e-5);
            Assert.AreEqual(54342, m12, 0.5);
            Assert.AreEqual(-1.0045592, M12, 0.5e-7);
            Assert.AreEqual(-0.9954339, M21, 0.5e-7);
            Assert.AreEqual(127516405431022.0, S12, 0.5);
        }

        [TestMethod]
        public void GeodSolve66_CheckForInverseLineIfLineIsSlightlyWestOfSAndThat_s12_IsCorrecltySet()
        {
            var line = Geodesic.WGS84.InverseLine(-5, -0.000000000000002, -10, 180);
            line.GenPosition(false, 2e7,
                GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth | GeodesicFlags.LongUnroll,
                out var lat2, out var lon2, out var azi2, out _, out _, out _, out _, out _);
            Assert.AreEqual(4.96445, lat2, 0.5e-5);
            Assert.AreEqual(-180.00000, lon2, 0.5e-5);
            Assert.AreEqual(-0.00000, azi2, 0.5e-5);

            line.GenPosition(false, 0.5 * line.Distance,
                GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth | GeodesicFlags.LongUnroll,
                out lat2, out lon2, out azi2, out _, out _, out _, out _, out _);
            Assert.AreEqual(-87.52461, lat2, 0.5e-5);
            Assert.AreEqual(-0.00000, lon2, 0.5e-5);
            Assert.AreEqual(-180.00000, azi2, 0.5e-5);
        }

        [TestMethod]
        public void GeodSolve71_CheckThatDirectLineSets_s13()
        {
            var line = Geodesic.WGS84.DirectLine(1, 2, 45, 1e7);
            line.GenPosition(false, 0.5 * line.Distance,
                GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth | GeodesicFlags.LongUnroll,
                out var lat2, out var lon2, out var azi2, out _, out _, out _, out _, out _);

            Assert.AreEqual(30.92625, lat2, 0.5e-5);
            Assert.AreEqual(37.54640, lon2, 0.5e-5);
            Assert.AreEqual(55.43104, azi2, 0.5e-5);
        }

        [TestMethod]
        public void GeodSolve73_CheckForBackwardsFromThePoleBug()
        {
            Geodesic.WGS84.Direct(90, 10, 180, -1e6, out var lat2, out var lon2, out var azi2);
            Assert.AreEqual(81.04623, lat2, 0.5e-5);
            Assert.AreEqual(-170, lon2, 0.5e-5);
            Assert.AreEqual(0, azi2, 0.5e-5);
            Assert.IsTrue(CopySign(1, azi2) > 0);
        }

        [TestMethod]
        public void GeodSolve74_CheckFixForInaccurateAreas()
        {
            var a12 = Geodesic.WGS84.Inverse(54.1589, 15.3872, 54.1591, 15.3877,
                out var s12, out var azi1, out var azi2, out var m12, out var M12, out var M21, out var S12);
            Assert.AreEqual(55.723110355, azi1, 5e-9);
            Assert.AreEqual(55.723515675, azi2, 5e-9);
            Assert.AreEqual(39.527686385, s12, 5e-9);
            Assert.AreEqual(0.000355495, a12, 5e-9);
            Assert.AreEqual(39.527686385, m12, 5e-9);
            Assert.AreEqual(0.999999995, M12, 5e-9);
            Assert.AreEqual(0.999999995, M21, 5e-9);
            Assert.AreEqual(286698586.30197, S12, 5e-4);
        }

        [TestMethod]
        public void GeodSolve76_TheDistanceFromWellingtonToSalamanca()
        {
            Geodesic.WGS84.Inverse(-(41 + 19 / 60.0), 174 + 49 / 60.0, 40 + 58 / 60.0, -(5 + 30 / 60.0),
                out var s12, out var azi1, out var azi2);

            Assert.AreEqual(160.39137649664, azi1, 0.5e-11);
            Assert.AreEqual(19.50042925176, azi2, 0.5e-11);
            Assert.AreEqual(19960543.857179, s12, 0.5e-6);
        }

        [TestMethod]
        public void GeodSolve78_AnExampleWhereTheNGSCalculatorFailsToConverge()
        {
            Geodesic.WGS84.Inverse(27.2, 0.0, -27.1, 179.5,
                out var s12, out var azi1, out var azi2);

            Assert.AreEqual(45.82468716758, azi1, 0.5e-11);
            Assert.AreEqual(134.22776532670, azi2, 0.5e-11);
            Assert.AreEqual(19974354.765767, s12, 0.5e-6);
        }

        [TestMethod]
        public void GeodSolve80_Coverage_ComputingScaleInSpecialCases_ZeroLengthGeodesic_UsingAnIncapableLine()
        {
            Geodesic.WGS84.GenInverse(0, 0, 0, 90, GeodesicFlags.GeodesicScale,
                out _, out _, out _, out _, out var M12, out var M21, out _);
            Assert.AreEqual(-0.00528427534, M12, 0.5e-10);
            Assert.AreEqual(-0.00528427534, M21, 0.5e-10);

            Geodesic.WGS84.GenInverse(0, 0, 1e-6, 1e-6, GeodesicFlags.GeodesicScale,
                out _, out _, out _, out _, out M12, out M21, out _);
            Assert.AreEqual(1, M12, 0.5e-10);
            Assert.AreEqual(1, M21, 0.5e-10);

            var a12 = Geodesic.WGS84.Inverse(20.001, 0, 20.001, 0,
                out var s12, out var azi1, out var azi2, out var m12, out M12, out M21, out var S12);
            Assert.AreEqual(0, a12, 1e-13);
            Assert.AreEqual(0, s12, 1e-8);
            Assert.AreEqual(180, azi1, 1e-13);
            Assert.AreEqual(180, azi2, 1e-13);
            Assert.AreEqual(0, m12, 1e-8);
            Assert.AreEqual(1, M12, 1e-15);
            Assert.AreEqual(1, M21, 1e-15);
            Assert.AreEqual(0, S12, 1e-10);

            a12 = Geodesic.WGS84.Inverse(90, 0, 90, 180,
                out s12, out azi1, out azi2, out m12, out M12, out M21, out S12);
            Assert.AreEqual(0, a12, 1e-13);
            Assert.AreEqual(0, s12, 1e-8);
            Assert.AreEqual(0, azi1, 1e-13);
            Assert.AreEqual(180, azi2, 1e-13);
            Assert.AreEqual(0, m12, 1e-8);
            Assert.AreEqual(1, M12, 1e-15);
            Assert.AreEqual(1, M21, 1e-15);
            Assert.AreEqual(127516405431022.0, S12, 0.5);

            var line = Geodesic.WGS84.Line(1, 2, 90, GeodesicFlags.Latitude);
            a12 = line.Position(1000, out _, out _);
            Assert.IsTrue(double.IsNaN(a12));
        }

        [TestMethod]
        public void GeodSolve84_CheckFixForRangeErrorsWith_fmod_sin_cos_inf()
        {
            Geodesic.WGS84.Direct(0, 0, 90, double.PositiveInfinity, out var lat2, out var lon2, out var azi2);
            Assert.IsTrue(double.IsNaN(lat2));
            Assert.IsTrue(double.IsNaN(lon2));
            Assert.IsTrue(double.IsNaN(azi2));

            Geodesic.WGS84.Direct(0, 0, 90, double.NaN, out lat2, out lon2, out azi2);
            Assert.IsTrue(double.IsNaN(lat2));
            Assert.IsTrue(double.IsNaN(lon2));
            Assert.IsTrue(double.IsNaN(azi2));

            Geodesic.WGS84.Direct(0, 0, double.PositiveInfinity, 1000, out lat2, out lon2, out azi2);
            Assert.IsTrue(double.IsNaN(lat2));
            Assert.IsTrue(double.IsNaN(lon2));
            Assert.IsTrue(double.IsNaN(azi2));

            Geodesic.WGS84.Direct(0, 0, double.NaN, 1000, out lat2, out lon2, out azi2);
            Assert.IsTrue(double.IsNaN(lat2));
            Assert.IsTrue(double.IsNaN(lon2));
            Assert.IsTrue(double.IsNaN(azi2));

            Geodesic.WGS84.Direct(0, double.PositiveInfinity,90, 1000, out lat2, out lon2, out azi2);
            Assert.AreEqual(0, lat2);
            Assert.IsTrue(double.IsNaN(lon2));
            Assert.AreEqual(90, azi2);

            Geodesic.WGS84.Direct(0, double.NaN, 90, 1000, out lat2, out lon2, out azi2);
            Assert.AreEqual(0, lat2);
            Assert.IsTrue(double.IsNaN(lon2));
            Assert.AreEqual(90, azi2);

            Geodesic.WGS84.Direct(double.PositiveInfinity, 0, 90, 1000, out lat2, out lon2, out azi2);
            Assert.IsTrue(double.IsNaN(lat2));
            Assert.IsTrue(double.IsNaN(lon2));
            Assert.IsTrue(double.IsNaN(azi2));

            Geodesic.WGS84.Direct(double.NaN, 0, 90, 1000, out lat2, out lon2, out azi2);
            Assert.IsTrue(double.IsNaN(lat2));
            Assert.IsTrue(double.IsNaN(lon2));
            Assert.IsTrue(double.IsNaN(azi2));
        }

        [TestMethod]
        public void GeodSolve92_CheckFixForInaccurateHypotWithPython()
        {
            var result = Geodesic.WGS84.Inverse(37.757540000000006, -122.47018, 37.75754, -122.470177);
            var resultExact = GeodesicExact.WGS84.Inverse(37.757540000000006, -122.47018, 37.75754, -122.470177);

            Assert.AreEqual(89.9999992, result.Azimuth1, 1e-7);
            Assert.AreEqual(90.000001, result.Azimuth2, 1e-6);
            Assert.AreEqual(0.264, result.Distance, 1e-3);

            Assert.AreEqual(89.9999992, resultExact.Azimuth1, 1e-7);
            Assert.AreEqual(90.000001, resultExact.Azimuth2, 1e-6);
            Assert.AreEqual(0.264, resultExact.Distance, 1e-3);
        }

        [TestMethod]
        public void GeodSolve94_CheckFixFor_lat2_eq_NaN_BeingTreatedAs_lat2_eq_Zero()
        {
            var result = Geodesic.WGS84.Inverse(0, 0, double.NaN, 90);

            Assert.IsTrue(double.IsNaN(result.Azimuth1));
            Assert.IsTrue(double.IsNaN(result.Azimuth2));
            Assert.IsTrue(double.IsNaN(result.Distance));
        }

        [TestMethod]
        public void GeodSolve94_CheckFixFor_lat2_eq_NaN_BeingTreatedAs_lat2_eq_Zero_Exact()
        {
            var result = GeodesicExact.WGS84.Inverse(0, 0, double.NaN, 90);

            Assert.IsTrue(double.IsNaN(result.Azimuth1));
            Assert.IsTrue(double.IsNaN(result.Azimuth2));
            Assert.IsTrue(double.IsNaN(result.Distance));
        }

        /// <summary>
        /// Failure with long doubles found with test case from Nowak + Nowak Da
        /// Costa (2022). Problem was using somg12 > 1 as a test that it needed
        /// to be set when roundoff could result in somg12 slightly bigger that 1.
        /// Found + fixed 2022-03-30.
        /// </summary>
        [TestMethod]
        public void GeodSolve96_97()
        {
            IGeodesic g = new Geodesic(6378137, 1 / 298.257222101);
            var result = g.Inverse(0, 0, 60.0832522871723, 89.8492185074635);
            Assert.AreEqual("42426932221845", result.Area.ToFixedString(0));

            g = new GeodesicExact(6378137, 1 / 298.257222101);
            result = g.Inverse(0, 0, 60.0832522871723, 89.8492185074635);
            Assert.AreEqual("42426932221845", result.Area.ToFixedString(0));
        }

        /// <summary>
        /// Area of on high eccentricity n = 0.94 ellipsoid.  This is the most
        /// oblate case that Nowak + Nowak Da Costa (2022) treats (convergence
        /// with kmax = 18504).  High precision result is 5910062452739.93899501;
        /// Accept xxx.930 thru xx.949
        /// </summary>
        [TestMethod]
        public void GeodSolve98()
        {
            var result = new GeodesicExact(6.4e6, 188 / 194d)
                .ArcDirect(0, 0, 45, 90, GeodesicFlags.Area | GeodesicFlags.LongUnroll);

            Assert.AreEqual(5910062452739.9395, result.Area, 0.0095);
        }
    }
}
