using Microsoft.VisualStudio.TestTools.UnitTesting;
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