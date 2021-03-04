using GeographicLib;
using GeographicLib.Geocodes;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib.Tests
{
    [TestClass]
    public class MGRSTests
    {
        [TestMethod]
        public void UTMUPSSantityCheck()
        {
            MGRS.Check();
        }

        [DataTestMethod]
        [DataRow(33.3, 44.4, "38SMB4414084706")]
        [DataRow(0, 0, "31NAA6602100000")]
        [DataRow(90, 180, "ZAH0000000000")]
        [DataRow(90, 0, "ZAH0000000000")]
        [DataRow(-90, 180, "BAN0000000000")]
        [DataRow(-90, 0, "BAN0000000000")]
        [DataRow(0, -180, "01NAA6602100000")]
        [DataRow(12.3456, 78.9012, "44PKU7177565663")]
        public void TestForward(double lat, double lon, string expeted)
        {
            var (zone, northp, x, y) = UTMUPS.Forward(lat, lon);
            var mgrs = MGRS.Forward(zone, northp, x, y, 5);

            Assert.AreEqual(expeted, mgrs);
        }

        [DataTestMethod]
        [DataRow(33.3, 44.4, "38SMB4414084706")]
        [DataRow(0, 0, "31NAA6602100000")]
        [DataRow(90, 180, "ZAH0000000000")]
        [DataRow(90, 0, "ZAH0000000000")]
        [DataRow(-90, 180, "BAN0000000000")]
        [DataRow(-90, 0, "BAN0000000000")]
        [DataRow(0, -180, "01NAA6602100000")]
        [DataRow(12.3456, 78.9012, "44PKU7177565663")]
        public void TestReverse(double lat, double lon, string input)
        {
            var (zone, northp, x, y, _) = MGRS.Reverse(input);
            var (olat, olon) = UTMUPS.Reverse(zone, northp, x, y);

            Assert.AreEqual(olat, lat, 1e-5);

            // No need to check longitude when near pole.
            if (Math.Abs(lat) != 90)
            {
                Assert.AreEqual(olon, lon, 1e-5);
            }
        }
    }
}
