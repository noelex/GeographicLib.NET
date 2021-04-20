using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using GeographicLib.Projections;

namespace GeographicLib.Tests
{
    [TestClass]
    public class ConicProjTests
    {
        [TestMethod]
        public void ConicProj0_CheckFixForAlbersEqualArea_ReverseBug()
        {
            var aea = new AlbersEqualArea(Ellipsoid.WGS84, DMS.DecodeAzimuth("40d58"), DMS.DecodeAzimuth("39d56"), 1);
            var (lat, lon) = aea.Reverse(DMS.DecodeAzimuth("77d45W"), 220e3, -52e3, out var gamma, out var k);

            Assert.AreEqual(39.95, lat, 0.01);
            Assert.AreEqual(-75.17, lon, 0.01);
            Assert.AreEqual(1.67, gamma, 0.01);
            Assert.AreEqual(0.99, k, 0.01);
        }

        [TestMethod]
        public void ConicProj1_CheckFixForAlbersEqualAreaProlateBug()
        {
            var aea = new AlbersEqualArea(6.4e6, -0.5, 0, 0, 1);
            var (lat, _) = aea.Reverse(0, 0, 8605508);

            Assert.AreEqual(85.00, lat, 0.01);
        }

        [TestMethod]
        public void ConicProj2_CheckFixForLambertConformalConic_ForwardBug()
        {
            var lcc = new LambertConformalConic(Ellipsoid.WGS84, -30, -30, 1);
            var (x, y) = lcc.Forward(0, -30, 0, out var gamma, out var k);

            Assert.AreEqual(0.0, Math.Abs(x), 0.01);
            Assert.AreEqual(0.0, Math.Abs(y), 0.01);
            Assert.AreEqual(0.0, Math.Abs(gamma), 0.01);
            Assert.AreEqual(1.0, k, 0.01);
        }

        [TestMethod]
        public void ConicProj3_CheckFixForLambertConformalConic_ReverseOverflowBugs()
        {
            var lcc = new LambertConformalConic(Ellipsoid.WGS84, 0, 0, 1);
            var (lat, lon) = lcc.Reverse(0, 1113195, -1e10);

            Assert.AreEqual(-90.0, lat, 0.01);
            Assert.AreEqual(10.00, lon, 0.001);
        }

        [TestMethod]
        public void ConicProj4()
        {
            var lcc = new LambertConformalConic(Ellipsoid.WGS84, 0, 0, 1);
            var (lat, lon) = lcc.Reverse(0, 1113195, double.PositiveInfinity);

            Assert.AreEqual(90.0, lat, 0.01);
            Assert.AreEqual(10.00, lon, 0.001);
        }

        [TestMethod]
        public void ConicProj5()
        {
            var lcc = new LambertConformalConic(Ellipsoid.WGS84, 45, 45, 1);
            var (lat, lon) = lcc.Reverse(0, 0, -1e100);

            Assert.AreEqual(-90.0, lat, 0.01);
            Assert.AreEqual(0.00, Math.Abs(lon), 0.001);
        }

        [TestMethod]
        public void ConicProj6()
        {
            var lcc = new LambertConformalConic(Ellipsoid.WGS84, 45, 45, 1);
            var (lat, lon) = lcc.Reverse(0, 0, double.NegativeInfinity);

            Assert.AreEqual(-90.0, lat, 0.01);
            Assert.AreEqual(0.00, Math.Abs(lon), 0.001);
        }

        [TestMethod]
        public void ConicProj7()
        {
            var lcc = new LambertConformalConic(Ellipsoid.WGS84, 90, 90, 1);
            var (lat, lon) = lcc.Reverse(0, 0, -1e150);

            Assert.AreEqual(-90.0, lat, 0.01);
            Assert.AreEqual(0.00, Math.Abs(lon), 0.001);
        }

        [TestMethod]
        public void ConicProj8()
        {
            var lcc = new LambertConformalConic(Ellipsoid.WGS84, 90, 90, 1);
            var (lat, lon) = lcc.Reverse(0, 0, double.NegativeInfinity);

            Assert.AreEqual(-90.0, lat, 0.01);
            Assert.AreEqual(0.00, Math.Abs(lon), 0.001);
        }

        [TestMethod]
        public void ConicProj9_CheckFixToInfiniteLoopInAlbersEqualAreaWith_e2_LessThan_MinusOne()
        {
            var lcc = new AlbersEqualArea(6.4e6,-0.5, -10, 40, 1);
            var (x, y) = lcc.Forward(0, 85, 10);

            Assert.AreEqual(609861, x, 1);
            Assert.AreEqual(7566522, y, 1);
        }
    }
}
