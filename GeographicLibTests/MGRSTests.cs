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

        [DataTestMethod]
        [DataRow(+0d, 0d, true, "31NEA0000")]
        [DataRow(-0d, 10e6, false, "31MEV0099")]
        public void TestForward_LatitudeSign(double lat, double y, bool northp, string mgrs)
        {
            var (zone, northpa, x, ya) = UTMUPS.Forward(lat, 3);
            Assert.That.EqualsExactly(y, ya);
            Assert.AreEqual(northp, northpa);

            var mgrsa = MGRS.Forward(zone, northp, x, y, 2);
            Assert.AreEqual(mgrs, mgrsa);

            mgrsa = MGRS.Forward(zone, northp, x, y, +0, 2);
            Assert.AreEqual(mgrs, mgrsa);

            mgrsa = MGRS.Forward(zone, northp, x, y, -0, 2);
            Assert.AreEqual(mgrs, mgrsa);
        }

        [DataTestMethod]
        [DataRow("38SMB4488", "38S", "MB", "44", "88")]
        [DataRow("3SMB4488", "3S", "MB", "44", "88")]
        [DataRow("SMB4488", "S", "MB", "44", "88")]
        [DataRow("SMB48", "S", "MB", "4", "8")]
        public void TestDecode(string mgrs, string gridzone, string block, string easting, string northing)
        {
            var (z, b, e, n) = MGRS.Decode(mgrs.AsSpan());
            Assert.AreEqual(gridzone, z);
            Assert.AreEqual(block, b);
            Assert.AreEqual(easting, e);
            Assert.AreEqual(northing, n);
        }

        [DataTestMethod]
        [DataRow("INV", "INV", "", "", "")]
        [DataRow("iNv", "iNv", "", "", "")]
        [DataRow("inv1234567", "inv", "", "", "")]
        public void TestDecode_InvHandling(string mgrs, string gridzone, string block, string easting, string northing)
        {
            var (z, b, e, n) = MGRS.Decode(mgrs.AsSpan());
            Assert.AreEqual(gridzone, z);
            Assert.AreEqual(block, b);
            Assert.AreEqual(easting, e);
            Assert.AreEqual(northing, n);
        }

        [DataTestMethod]
        [ExpectedException(typeof(GeographicException))]
        [DataRow("AABCDDD")]
        [DataRow("0011223")]
        [DataRow("0AB1223")]
        [DataRow("0ABC223")]
        [DataRow("0ABCD23")]
        [DataRow("00OOD23")]
        [DataRow("00IID23")]
        public void TestDecode_SanityCheck(string mgrs)
        {
            MGRS.Decode(mgrs.AsSpan());
        }
    }
}
