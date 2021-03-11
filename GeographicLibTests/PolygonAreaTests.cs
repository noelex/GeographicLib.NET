using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib.Tests
{
    [TestClass]
    public class PolygonAreaTests
    {
        [TestMethod]
        public void TestCompute()
        {
            var data = new[] { "18n 500000 4400000", "18n 600000 4400000", "18n 600000 4500000", "18n 500000 4500000" };

            var p = new PolygonArea(Geodesic.WGS84);

            p.AddPoints(data.Select(x => new GeoCoords(x)));
            var (c, peri, area) = p.Compute(false, false);

            Assert.AreEqual(4, c);
            Assert.AreEqual(400139.5329585969, peri, 1e-8);
            Assert.AreEqual(10007388597.19133, area, 1e-8);

        }
    }
}
