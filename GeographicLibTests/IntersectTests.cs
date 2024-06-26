﻿using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;

namespace GeographicLib.Tests
{

    [TestClass]
    public class IntersectTests
    {
        public static IEnumerable<object[]> EquatorialSegmentData
        {
            get
            {
                for (int exact = 0; exact < 2; ++exact)
                {
                    for (int fi = -1; fi <= 1; ++fi)
                    {
                        var setup = new object[] { fi / 10.0, exact != 0 };
                        yield return setup.Concat(new object[] { 0, 40, -20, -10, -5, 15 }).ToArray();
                        yield return setup.Concat(new object[] { 0, 40, -10, 0, 0, 10 }).ToArray();
                        yield return setup.Concat(new object[] { 0, 40, -8, 2, 1, 9 }).ToArray();
                        yield return setup.Concat(new object[] { 0, 40, -2, 8, 4, 6 }).ToArray();
                        yield return setup.Concat(new object[] { 0, 40, 0, 10, 5, 5 }).ToArray();
                        yield return setup.Concat(new object[] { 0, 40, 2, 12, 7, 5 }).ToArray();
                        yield return setup.Concat(new object[] { 0, 40, 15, 25, 20, 5 }).ToArray();
                        yield return setup.Concat(new object[] { 0, 40, 30, 40, 35, 5 }).ToArray();
                        yield return setup.Concat(new object[] { 0, 40, 32, 42, 36, 4 }).ToArray();
                        yield return setup.Concat(new object[] { 0, 40, 38, 48, 39, 1 }).ToArray();
                        yield return setup.Concat(new object[] { 0, 40, 40, 50, 40, 0 }).ToArray();
                        yield return setup.Concat(new object[] { 0, 40, 50, 60, 45, -5 }).ToArray();
                        yield return setup.Concat(new object[] { 40, 0, -20, -10, 40 + 5, 15 }).ToArray();
                        yield return setup.Concat(new object[] { 40, 0, -10, -0, 40, 10 }).ToArray();
                        yield return setup.Concat(new object[] { 40, 0, -8, 2, 40 - 1, 9 }).ToArray();
                        yield return setup.Concat(new object[] { 40, 0, -2, 8, 40 - 4, 6 }).ToArray();
                        yield return setup.Concat(new object[] { 40, 0, 0, 10, 40 - 5, 5 }).ToArray();
                        yield return setup.Concat(new object[] { 40, 0, 2, 12, 40 - 7, 5 }).ToArray();
                        yield return setup.Concat(new object[] { 40, 0, 15, 25, 40 - 20, 5 }).ToArray();
                        yield return setup.Concat(new object[] { 40, 0, 30, 40, 40 - 35, 5 }).ToArray();
                        yield return setup.Concat(new object[] { 40, 0, 32, 42, 40 - 36, 4 }).ToArray();
                        yield return setup.Concat(new object[] { 40, 0, 38, 48, 40 - 39, 1 }).ToArray();
                        yield return setup.Concat(new object[] { 40, 0, 40, 50, 40 - 40, 0 }).ToArray();
                        yield return setup.Concat(new object[] { 40, 0, 50, 60, 40 - 45, -5 }).ToArray();
                    }
                }
            }
        }

        [TestMethod]
        public void TestClosest()
        {
            Geodesic geod = Geodesic.WGS84;
            Intersect inter = new Intersect(geod);
            IGeodesicLine
                lineX = geod.Line(0, 0, 45),
                lineY = geod.Line(45, 10, 135);

            Point point = inter.Closest(lineX, lineY);
            lineX.Position(point.X, out double latx, out double lonx);
            lineY.Position(point.Y, out double laty, out double lony);

            Assert.AreEqual(latx, laty, 1e-12);
            Assert.AreEqual(lonx, lony, 1e-12);
        }

        [DataTestMethod]
        [DynamicData(nameof(EquatorialSegmentData), typeof(IntersectTests))]
        public void EquatorialSegments(
            double f, bool exact,
            double lonx1, double lonx2, double lony1, double lony2, double px, double py)
        {
            const double a = 180 / Math.PI, eps = 1 / 1000000.0;

            var geod = new Geodesic(a, f, exact);
            var inter = new Intersect(geod);

            var p1 = inter.Segment(0, lonx1, 0, lonx2,
                                   0, lony1, 0, lony2, out _);
            var p2 = inter.Segment(0, lony1, 0, lony2,
                                   0, lonx1, 0, lonx2, out _);

            Assert.AreEqual(p1.X, px, eps);
            Assert.AreEqual(p1.Y, py, eps);
            Assert.AreEqual(p2.X, py, eps);
            Assert.AreEqual(p2.Y, px, eps);
        }

        private (double latX, double lonX, double aziX,
                 double latY, double lonY, double aziY) Parse(string input)
        {
            var args = input.Split(' ');
            var (latX, lonX) = DMS.Decode(args[0], args[1]);
            var aziX = DMS.DecodeAzimuth(args[2]);
            var (latY, lonY) = DMS.Decode(args[3], args[4]);
            var aziY = DMS.DecodeAzimuth(args[5]);
            return (latX, lonX, aziX, latY, lonY, aziY);
        }

        [TestMethod]
        public void TestNext()
        {
            var coord = new GeoCoords("50N 4W");
            var aziX = DMS.DecodeAzimuth("147.7W");
            var aziY = DMS.DecodeAzimuth("45");

            var inter = new Intersect(Geodesic.WGS84);
            var p = inter.Next(coord.Latitude, coord.Longitude, aziX, aziY, out var c);

            Assert.AreEqual(-39901512.037, p.X, 1);
            Assert.AreEqual(-117095.317, p.Y, 1);
            Assert.AreEqual(0, c);
        }

        [DataTestMethod]
        [DataRow(0, -40212385.966, -0.0, 0)]
        [DataRow(-0.3, -21976861.924, 21975015.260, 0)]
        [DataRow(0.3, -23840207.577, -11920024.164, 0)]
        public void HiglyProlateOblateSphere(double f, double x, double y, double c)
        {
            var coord = new GeoCoords("50N 4W");
            var aziX = DMS.DecodeAzimuth("147.7W");
            var aziY = DMS.DecodeAzimuth("45");

            var inter = new Intersect(new Geodesic(6.4e6, f));
            var p = inter.Next(coord.Latitude, coord.Longitude, aziX, aziY, out var _c);

            Assert.AreEqual(x, p.X, 1);
            Assert.AreEqual(y, p.Y, 1);
            Assert.AreEqual(c, _c);
        }

        [DataTestMethod]
        [DataRow(-0.5)]
        [DataRow(0.5)]
        [ExpectedException(typeof(GeographicException))]
        public void TooEccentric(double f)
        {
            new Intersect(new Geodesic(6.4e6, f));
            Assert.Fail();
        }

        [TestMethod]
        public void Intersect1()
        {
            var (latX, lonX, aziX, latY, lonY, aziY) = Parse("50N 4W 147.7W 0 0 90");
            var inter = new Intersect(Geodesic.WGS84);
            var p = inter.Closest(latX, lonX, aziX, latY, lonY, aziY, out var c);
            Assert.AreEqual(6058049, p.X, 1);
            Assert.AreEqual(-3311253, p.Y, 1);
            Assert.AreEqual(0, c);
        }

        [TestMethod]
        public void Intersect2()
        {
            var (latX, lonX, aziX, latY, lonY, aziY) = Parse("50N 4W 147.7W 0 180 0");
            var inter = new Intersect(Geodesic.WGS84);
            var p = inter.All(latX, lonX, aziX, latY, lonY, aziY, 2.6e7, out var c);

            Assert.AreEqual(2, p.Length);

            Assert.AreEqual(-494582, p[0].X, 1);
            Assert.AreEqual(14052230, p[0].Y, 1);
            Assert.AreEqual(0, c[0]);
            Assert.AreEqual(14546812, Intersect.Dist(p[0]), 1);

            Assert.AreEqual(19529110, p[1].X, 1);
            Assert.AreEqual(-5932344, p[1].Y, 1);
            Assert.AreEqual(0, c[1]);
            Assert.AreEqual(25461454, Intersect.Dist(p[1]), 1);
        }
    }
}