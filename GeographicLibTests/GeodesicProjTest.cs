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
    public class GeodesicProjTest
    {
        [TestMethod]
        public void GeodesicProj0_TestFixToCassiniSoldner_ForwardBug()
        {
            var cass = new CassiniSoldner();
            var (x, _) = cass.Forward(90, 80, out var azi, out _);
            Assert.AreEqual(0.0, x, 0.1);
            Assert.AreEqual(170.0, azi, 0.1);
        }
    }
}
