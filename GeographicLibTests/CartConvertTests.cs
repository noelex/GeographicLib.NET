using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib.Tests
{
    [TestClass]
    public class CartConvertTests
    {
        [TestMethod]
        public void CartConvert0()
        {
            var cart = new Geocentric(6.4e6, 1.0 / 100);
            var (lat, lon, h) = cart.Reverse(10e3, 0, 1e3);

            Assert.AreEqual(85.57, lat, 0.01);
            Assert.AreEqual(0.00, lon, 0.001);
            Assert.AreEqual(-6334614, h, 1);
        }

        [TestMethod]
        public void CartConvert1()
        {
            var cart = new Geocentric(6.4e6, -1.0 / 100);
            var (lat, lon, h) = cart.Reverse(1e3, 0, 10e3);

            Assert.AreEqual(4.42, lat, 0.01);
            Assert.AreEqual(0.00, lon, 0.001);
            Assert.AreEqual(-6398614, h, 1);
        }
    }
}
