using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Collections.Generic;
using System.Linq;

namespace GeographicLib.Tests
{
    [TestClass]
    public class PolygonAreaTests
    {
        public static IEnumerable<object[]> PointData =>
            new[]
            {
                new object[]{ new[] { "18n 500000 4400000", "18n 600000 4400000", "18n 600000 4500000", "18n 500000 4500000" }, 4, 400139.5329585969,10007388597.19133  },
                new object[]{ new[] { "52 0", "41 -74", "-23 -43", "-26 28" },4, 29506941.1551780142, 65690027591345.67188 },
            };

        public static IEnumerable<object[]> EdgeData =>
            new[]
            {
                new object[]{ (0.0, 0.0), new[] { (90.0,1000.0), (0.0,1000.0), (-90.0, 1000.0) }, 4, 3999.9999876263041, 999999.99793803820 },
            };

        [DataTestMethod]
        [DynamicData("PointData", typeof(PolygonAreaTests))]
        public void TestComputeWithPoints(string[] coords, int n, double perimeter, double area)
        {
            var p = new PolygonArea(Geodesic.WGS84);

            p.AddPoints(coords.Select(x => new GeoCoords(x)));
            var (c, peri, a) = p.Compute(false, false);

            Assert.AreEqual(n, c);
            Assert.AreEqual(perimeter, peri, 1e-8);
            Assert.AreEqual(area, a, 1e-8);
        }

        [DataTestMethod]
        [DynamicData("EdgeData", typeof(PolygonAreaTests))]
        public void TestComputeWithEdges((double lat, double lon) startp, (double azi, double s)[] edges, int n, double perimeter, double area)
        {
            var p = new PolygonArea(Geodesic.WGS84);
            p.AddPoint(startp.lat, startp.lon);

            foreach (var (azimuth, distance) in edges)
            {
                p.AddEdge(azimuth, distance);
            }
            var (c, peri, a) = p.Compute(false, false);

            Assert.AreEqual(n, c);
            Assert.AreEqual(perimeter, peri, 1e-8);
            Assert.AreEqual(area, a, 1e-8);
        }

        [TestMethod]
        public void GH29_PolygonAreaExactGivesIncorrectResult()
        {
            var poly = new PolygonAreaExact(GeodesicExact.WGS84);
            poly.AddPoint(46.227634, 7.9019728);
            poly.AddPoint(46.21863761, 7.9019728);
            poly.AddPoint(46.21863761, 7.914933282);
            poly.AddPoint(46.227634, 7.914933282);

            var (_, _, area) = poly.Compute(false, false);

            Assert.AreEqual(1_000_000, area, 100);
        }
    }
}
