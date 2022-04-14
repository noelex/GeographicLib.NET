using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib.Tests
{
    [TestClass]
    public class PlanimeterTests
    {
        private readonly IPolygonArea 
            polygon = new PolygonArea(Geodesic.WGS84),
            polyline = new PolygonArea(Geodesic.WGS84, true),
            polygonR = new PolygonAreaRhumb(Rhumb.WGS84);

        private (int num, double perimeter, double area) Planimeter((double lat, double lon)[] points)
        {
            polygon.Clear();
            foreach (var (lat, lon) in points)
            {
                polygon.AddPoint(lat, lon);
            }

            return polygon.Compute(false, true);
        }

        private (int num, double perimeter, double area) PlanimeterRhumb((double lat, double lon)[] points)
        {
            polygonR.Clear();
            foreach (var (lat, lon) in points)
            {
                polygonR.AddPoint(lat, lon);
            }

            return polygonR.Compute(false, true);
        }

        [TestMethod]
        public void Planimeter0_CheckFixForPoleEncirclingBug()
        {
            var points = new[] { (89d, 0d), (89, 90), (89, 180), (89, 270) };
            var (_, perimeter, area) = Planimeter(points);
            Assert.AreEqual(631819.8745, perimeter, 1e-4);
            Assert.AreEqual(24952305678.0, area, 1);

            points = new[] { (-89d, 0d), (-89, 90), (-89, 180), (-89, 270) };
            (_, perimeter, area) = Planimeter(points);
            Assert.AreEqual(631819.8745, perimeter, 1e-4);
            Assert.AreEqual(-24952305678.0, area, 1);

            points = new[] { (0d, -1d), (-1, 0), (0, 1), (1, 0) };
            (_, perimeter, area) = Planimeter(points);
            Assert.AreEqual(627598.2731, perimeter, 1e-4);
            Assert.AreEqual(24619419146.0, area, 1);

            points = new[] { (90d, 0d), (0, 0), (0, 90) };
            (_, perimeter, area) = Planimeter(points);
            Assert.AreEqual(30022685, perimeter, 1);
            Assert.AreEqual(63758202715511.0, area, 1);
        }

        [TestMethod]
        public void Planimeter5_CheckFixForPlanimeterPoleCrossingBug()
        {
            var points = new[] { (89d, 0.1d), (89, 90.1), (89, -179.9) };
            var (_, perimeter, area) = Planimeter(points);
            Assert.AreEqual(539297, perimeter, 1);
            Assert.AreEqual(12476152838.5, area, 1);
        }

        [TestMethod]
        public void Planimeter6_CheckFixForPlanimeter_lon12_RoundingBug()
        {
            var points = new[] { (9d, -0.00000000000001d), (9, 180), (9, 0) };
            var (_, perimeter, area) = Planimeter(points);
            Assert.AreEqual(36026861, perimeter, 1);
            Assert.AreEqual(0, area, 1);

            points = new[] { (9d, 0.00000000000001d), (9, 0), (9, 180) };
            (_, perimeter, area) = Planimeter(points);
            Assert.AreEqual(36026861, perimeter, 1);
            Assert.AreEqual(0, area, 1);

            points = new[] { (9d, 0.00000000000001d), (9, 180), (9, 0) };
            (_, perimeter, area) = Planimeter(points);
            Assert.AreEqual(36026861, perimeter, 1);
            Assert.AreEqual(0, area, 1);

            points = new[] { (9d, -0.00000000000001d), (9, 0), (9, 180) };
            (_, perimeter, area) = Planimeter(points);
            Assert.AreEqual(36026861, perimeter, 1);
            Assert.AreEqual(0, area, 1);
        }

        [TestMethod]
        public void Planimeter11_AreaOfArcticCircle_Rhumb()
        {
            var points = new[] {
                (DMS.DecodeAngle("66:33:44"), 0d),
                (DMS.DecodeAngle("66:33:44"), 180),
                (DMS.DecodeAngle("66:33:44"), 360)
            };
            var (_, perimeter, area) = PlanimeterRhumb(points);
            Assert.AreEqual(15985058, perimeter, 1);
            Assert.AreEqual(21208418252300, area, 100);
        }

        [TestMethod]
        public void Planimeter11r_ReverseAreaOfArcticCircle_Rhumb()
        {
            var points = new[] {
                (DMS.DecodeAngle("66:33:44"), -0d),
                (DMS.DecodeAngle("66:33:44"), -180),
                (DMS.DecodeAngle("66:33:44"), -360)
            };
            var (_, perimeter, area) = PlanimeterRhumb(points);
            Assert.AreEqual(15985058, perimeter, 1);
            Assert.AreEqual(-21208418252300, area, 100);
        }

        [TestMethod]
        public void Planimeter12_AreaOfArcticCircle()
        {
            var points = new[] { 
                (DMS.DecodeAngle("66:33:44"), 0d),
                (DMS.DecodeAngle("66:33:44"), 180),
                (DMS.DecodeAngle("66:33:44"), 360)
            };
            var (_, perimeter, area) = Planimeter(points);
            Assert.AreEqual(10465729, perimeter, 1);
            Assert.AreEqual(0, area);
        }

        [TestMethod]
        public void Planimeter12r_ReverseAreaOfArcticCircle()
        {
            var points = new[] {
                (DMS.DecodeAngle("66:33:44"), -0d),
                (DMS.DecodeAngle("66:33:44"), -180),
                (DMS.DecodeAngle("66:33:44"), -360)
            };
            var (_, perimeter, area) = Planimeter(points);
            Assert.AreEqual(10465729, perimeter, 1);
            Assert.AreEqual(0, area);
        }

        [TestMethod]
        public void Planimeter13_CheckEncirclingPoleTwice()
        {
            var points = new[] { (89d, -360d), (89, -240), (89, -120), (89, 0), (89, 120), (89, 240) };
            var (_, perimeter, area) = Planimeter(points);
            Assert.AreEqual(1160741, perimeter, 1);
            Assert.AreEqual(32415230256.0, area, 1);
        }

        [TestMethod]
        public void Planimeter15_Coverage_Planimeter15ToPlanimeter18()
        {
            var lat = new[] { 2d, 1, 3 };
            var lon = new[] { 1d, 2, 3 };
            var r = 18454562325.45119;
            var a0 = 510065621724088.5093;
            polygon.Clear();
            polygon.AddPoint(lat[0], lon[0]);
            polygon.AddPoint(lat[1], lon[1]);

            var (_, _, area) = polygon.TestPoint(lat[2], lon[2], false, true);
            Assert.AreEqual(r, area, 0.5);
            (_, _, area) = polygon.TestPoint(lat[2], lon[2], false, false);
            Assert.AreEqual(r, area, 0.5);
            (_, _, area) = polygon.TestPoint(lat[2], lon[2], true, true);
            Assert.AreEqual(-r, area, 0.5);
            (_, _, area) = polygon.TestPoint(lat[2], lon[2], true, false);
            Assert.AreEqual(a0 - r, area, 0.5);

            Geodesic.WGS84.Inverse(lat[1], lon[1], lat[2], lon[2], out var s12, out var azi1, out _);

            (_, _, area) = polygon.TestEdge(azi1, s12, false, true);
            Assert.AreEqual(r, area, 0.5);
            (_, _, area) = polygon.TestEdge(azi1, s12, false, false);
            Assert.AreEqual(r, area, 0.5);
            (_, _, area) = polygon.TestEdge(azi1, s12, true, true);
            Assert.AreEqual(-r, area, 0.5);
            (_, _, area) = polygon.TestEdge(azi1, s12, true, false);
            Assert.AreEqual(a0 - r, area, 0.5);

            polygon.AddPoint(lat[2], lon[2]);

            (_, _, area) = polygon.Compute(false, true);
            Assert.AreEqual(r, area, 0.5);
            (_, _, area) = polygon.Compute(false, false);
            Assert.AreEqual(r, area, 0.5);
            (_, _, area) = polygon.Compute(true, true);
            Assert.AreEqual(-r, area, 0.5);
            (_, _, area) = polygon.Compute(true, false);
            Assert.AreEqual(a0 - r, area, 0.5);
        }

        [TestMethod]
        public void Planimeter19_Coverage_Planimeter19ToPlanimeter20()
        {
            polygon.Clear();
            var (_, perimeter, area) = polygon.Compute(false, true);
            Assert.AreEqual(0, area);
            Assert.AreEqual(0, perimeter);

            (_, perimeter, area) = polygon.TestPoint(1, 1, false, true);
            Assert.AreEqual(0, area);
            Assert.AreEqual(0, perimeter);

            (_, perimeter, area) = polygon.TestEdge(90, 1000, false, true);
            Assert.IsTrue(double.IsNaN(area));
            Assert.IsTrue(double.IsNaN(perimeter));

            polygon.AddPoint(1, 1);
            (_, perimeter, area) = polygon.Compute(false, true);
            Assert.AreEqual(0, area);
            Assert.AreEqual(0, perimeter);




            polyline.Clear();
            (_, perimeter, _) = polyline.Compute(false, true);
            Assert.AreEqual(0, perimeter);

            (_, perimeter, _) = polyline.TestPoint(1, 1, false, true);
            Assert.AreEqual(0, perimeter);

            (_, perimeter, _) = polyline.TestEdge(90, 1000, false, true);
            Assert.IsTrue(double.IsNaN(perimeter));

            polyline.AddPoint(1, 1);
            (_, perimeter, _) = polyline.Compute(false, true);
            Assert.AreEqual(0, perimeter);




            polygon.AddPoint(1, 1);
            (_, perimeter, _) = polyline.TestEdge(90, 1000, false, true);
            Assert.AreEqual(1000, perimeter, 1e-10);

            (_, perimeter, _) = polyline.TestPoint(2, 2, false, true);
            Assert.AreEqual(156876.149, perimeter, 0.5e-3);
        }

        [TestMethod]
        public void Planimeter21_Coverage_MultipleCirclingsOfPole_InvocationsViaTestPointAndTestEdge()
        {
            double lat = 45,
                 azi = 39.2144607176828184218,
                 s = 8420705.40957178156285,
                 r = 39433884866571.4277,   // Area for one circuit
                 a0 = 510065621724088.5093;  // Ellipsoid area

            polygon.Clear();
            polygon.AddPoints(new[] {
                new GeoCoords(lat,60),
                new GeoCoords(lat,180),
                new GeoCoords(lat,-60),
                new GeoCoords(lat,60),
                new GeoCoords(lat,180),
                new GeoCoords(lat,-60),
            });

            for(var i = 3; i <= 4; i++)
            {
                polygon.AddPoint(lat, 60);
                polygon.AddPoint(lat, 180);

                var (_, _, area) = polygon.TestPoint(lat, -60, false, true);
                Assert.AreEqual(area, i * r, 0.5);
                (_, _, area) = polygon.TestPoint(lat, -60, false, false);
                Assert.AreEqual(area, i * r, 0.5);
                (_, _, area) = polygon.TestPoint(lat, -60, true, true);
                Assert.AreEqual(area, -i * r, 0.5);
                (_, _, area) = polygon.TestPoint(lat, -60, true, false);
                Assert.AreEqual(area, -i * r + a0, 0.5);

                (_, _, area) = polygon.TestEdge(azi, s, false, true);
                Assert.AreEqual(area, i * r, 0.5);
                (_, _, area) = polygon.TestEdge(azi, s, false, false);
                Assert.AreEqual(area, i * r, 0.5);
                (_, _, area) = polygon.TestEdge(azi, s, true, true);
                Assert.AreEqual(area, -i * r, 0.5);
                (_, _, area) = polygon.TestEdge(azi, s, true, false);
                Assert.AreEqual(area, -i * r + a0, 0.5);

                polygon.AddPoint(lat, -60);

                (_, _, area) = polygon.Compute(false, true);
                Assert.AreEqual(area, i * r, 0.5);
                (_, _, area) = polygon.Compute(false, false);
                Assert.AreEqual(area, i * r, 0.5);
                (_, _, area) = polygon.Compute(true, true);
                Assert.AreEqual(area, -i * r, 0.5);
                (_, _, area) = polygon.Compute(true, false);
                Assert.AreEqual(area, -i * r + a0, 0.5);
            }
        }

        [TestMethod]
        public void Planimeter29_CheckFixToTransitDirectVSTransitZeroHandlingInconsistency()
        {
            polygon.Clear();
            polygon.AddPoint(0, 0);
            polygon.AddEdge(90, 1000);
            polygon.AddEdge(0, 1000);
            polygon.AddEdge(-90, 1000);

            var (_, _, area) = polygon.Compute(false, true);
            Assert.AreEqual(1000000.0, area, 0.01);
        }
    }
}
