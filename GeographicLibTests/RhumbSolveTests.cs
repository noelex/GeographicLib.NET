using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib.Tests
{
    [TestClass]
    public class RhumbSolveTests
    {
        /// <summary>
        /// RhumbSovle0 - RhumbSolve4
        /// </summary>
        [TestMethod]
        public void RhumbSolve0_TestFixToBadHandlingOfPole()
        {
            var rhe= new Rhumb(Ellipsoid.WGS84);
            var rh = new Rhumb(Ellipsoid.WGS84, false);

            rhe.Inverse(0, 0, 90, 0, out var s12, out var azi12, out var S12);
            Assert.AreEqual(0, azi12, 0.01);
            Assert.AreEqual(10001965.729, s12, 0.001);
            Assert.AreEqual(0, S12, 0.1);

            rh.Inverse(0, 0, 90, 0, out s12, out azi12, out S12);
            Assert.AreEqual(0, azi12, 0.01);
            Assert.AreEqual(10001965.729, s12, 0.001);
            Assert.AreEqual(0, S12, 0.1);

            rhe.Inverse(0, 0, 90, 30, out s12, out azi12, out S12);
            Assert.AreEqual(0.41222947, azi12, 1e-8);
            Assert.AreEqual(10002224.609, s12, 0.001);
            Assert.AreEqual(21050634712282, S12, 0.1);

            rh.Inverse(0, 0, 90, 30, out s12, out azi12, out S12);
            Assert.AreEqual(0.41222947, azi12, 1e-8);
            Assert.AreEqual(10002224.609, s12, 0.001);
            Assert.AreEqual(21050634712282, S12, 0.1);
        }
    }
}
