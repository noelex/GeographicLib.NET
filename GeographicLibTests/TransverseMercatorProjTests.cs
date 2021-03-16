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
    public class TransverseMercatorProjTests
    {
        [TestMethod]
        public void TransverseMercatorProj0_TestFixToBadMeridianConvergence()
        {
            var tme = new TransverseMercatorExact(Ellipsoid.WGS84, 1);
            var (x, y) = tme.Forward(0, 90, 75, out var gamma, out var k);
            Assert.AreEqual(0.0, x, 0.01);
            Assert.AreEqual(10001965.7293, y, 0.0001);
            Assert.AreEqual(75.0, gamma, 0.01);
            Assert.AreEqual(1.0, k, 0.01);

            var tm = new TransverseMercator(Ellipsoid.WGS84, 1);
            (x, y) = tm.Forward(0, 90, 75, out gamma, out k);
            Assert.AreEqual(0.0, x, 0.01);
            Assert.AreEqual(10001965.7293, y, 0.0001);
            Assert.AreEqual(75.0, gamma, 0.01);
            Assert.AreEqual(1.0, k, 0.01);
        }

        [TestMethod]
        public void TransverseMercatorProj2_TestFixToBadScaleAtPoleWithTransverseMercatorExact()
        {
            var tme = new TransverseMercatorExact(Ellipsoid.WGS84, 1);
            var (lat, lon) = tme.Reverse(0, 0, 10001965.7293127228, out var gamma, out var k);
            Assert.AreEqual(90.0, lat, 0.01);
            Assert.AreEqual(0, lon, 0.01);
            Assert.AreEqual(0, gamma, 0.01);
            Assert.AreEqual(1, k, 1e-5);

            var tm = new TransverseMercator(Ellipsoid.WGS84, 1);
            (lat, lon) = tm.Reverse(0, 0, 10001965.7293127228, out gamma, out k);
            Assert.AreEqual(90.0, lat, 0.01);
            Assert.AreEqual(180, Math.Abs(lon), 0.01);
            Assert.AreEqual(180, Math.Abs(gamma), 0.01);
            Assert.AreEqual(1, k, 1e-5);
        }

        [TestMethod]
        public void TransverseMercatorProj4_CheckUseOfComplexArithmeticToDoClenshawSum()
        {
            var tme = new TransverseMercatorExact(6.4e6, 1d/150, Constants.UTM_k0);
            var (x, y) = tme.Forward(0, 20,30, out var gamma, out var k);
            Assert.AreEqual(3266035.453860, x, 1e-6);
            Assert.AreEqual(2518371.552676, y, 1e-6);
            Assert.AreEqual(11.207356502141, gamma, 1e-12);
            Assert.AreEqual(1.134138960741, k, 1e-12);

            var tm = new TransverseMercator(6.4e6, 1d / 150, Constants.UTM_k0);
            (x, y) = tm.Forward(0, 20, 30, out gamma, out k);
            Assert.AreEqual(3266035.453860, x, 1e-6);
            Assert.AreEqual(2518371.552676, y, 1e-6);
            Assert.AreEqual(11.207356502141, gamma, 1e-12);
            Assert.AreEqual(1.134138960741, k, 1e-12);
        }

        [TestMethod]
        public void TransverseMercatorProj6()
        {
            var tme = new TransverseMercatorExact(6.4e6, 1d / 150, Constants.UTM_k0);
            var (lat, lon) = tme.Reverse(0, 3.3e6, 2.5e6, out var gamma, out var k);
            Assert.AreEqual(19.80370996793, lat, 1e-11);
            Assert.AreEqual(30.24919702282, lon, 1e-11);
            Assert.AreEqual(11.214378172893, gamma, 1e-12);
            Assert.AreEqual(1.137025775759, k, 1e-12);

            var tm = new TransverseMercator(6.4e6, 1d / 150, Constants.UTM_k0);
            (lat, lon) = tm.Reverse(0, 3.3e6, 2.5e6, out gamma, out k);
            Assert.AreEqual(19.80370996793, lat, 1e-11);
            Assert.AreEqual(30.24919702282, lon, 1e-11);
            Assert.AreEqual(11.214378172893, gamma, 1e-12);
            Assert.AreEqual(1.137025775759, k, 1e-12);
        }
    }
}