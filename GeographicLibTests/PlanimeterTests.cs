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
        private PolygonArea polygon = new PolygonArea(Geodesic.WGS84);
        private PolygonArea polyline = new PolygonArea(Geodesic.WGS84, false);

        private (int num, double perimeter, double area) Planimeter((double lat, double lon)[] points)
        {
            polygon.Clear();
            foreach(var (lat,lon) in points)
            {
                polygon.AddPoint(lat, lon);
            }

            return polygon.Compute(false, true);
        }

        [TestMethod]
        public void TestPlanimeter0_CheckFixForPoleEncirclingBug()
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
        public void TestPlanimeter5_CheckFixForPlanimeterPoleCrossingBug()
        {
            var points = new[] { (89d, 0.1d), (89, 90.1), (89, -179.9) };
            var (_, perimeter, area) = Planimeter(points);
            Assert.AreEqual(539297, perimeter, 1);
            Assert.AreEqual(12476152838.5, area, 1);
        }

        [TestMethod]
        public void TestPlanimeter6_CheckFixForPlanimeter_lon12_RoundingBug()
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
        public void TestPlanimeter12_AreaOfArcticCircle()
        {
            var points = new[] { (66.562222222, 0d), (66.562222222, 180) };
            var (_, perimeter, area) = Planimeter(points);
            Assert.AreEqual(10465729, perimeter, 1);
            Assert.AreEqual(0, area, 1);
        }

        [TestMethod]
        public void TestPlanimeter13_CheckEncirclingPoleTwice()
        {
            var points = new[] { (89d, -360d), (89, -240), (89, -120), (89, 0), (89, 120), (89, 240) };
            var (_, perimeter, area) = Planimeter(points);
            Assert.AreEqual(1160741, perimeter, 1);
            Assert.AreEqual(32415230256.0, area, 1);
        }

        [TestMethod]
        public void TestPlanimeter15_Coverage_Planimeter15_Planimeter18()
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
            Assert.AreEqual(a0-r, area, 0.5);

            Geodesic.WGS84.Inverse(lat[1], lon[1], lat[2], lon[2], out var s12, out var azi1, out _);

            (_, _, area) = polygon.TestEdge(azi1, s12, false, true);
            Assert.AreEqual(r, area, 0.5);
            (_, _, area) = polygon.TestEdge(azi1, s12, false, false);
            Assert.AreEqual(r, area, 0.5);
            (_, _, area) = polygon.TestEdge(azi1, s12, true, true);
            Assert.AreEqual(-r, area, 0.5);
            (_, _, area) = polygon.TestEdge(azi1, s12, true, false);
            Assert.AreEqual(a0-r, area, 0.5);

            polygon.AddPoint(lat[2], lon[2]);

            (_, _, area) = polygon.Compute(false, true);
            Assert.AreEqual(r, area, 0.5);
            (_, _, area) = polygon.Compute(false, false);
            Assert.AreEqual(r, area, 0.5);
            (_, _, area) = polygon.Compute(true, true);
            Assert.AreEqual(-r, area, 0.5);
            (_, _, area) = polygon.Compute(true, false);
            Assert.AreEqual(a0-r, area, 0.5);
        }
    }
}
