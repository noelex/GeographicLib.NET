using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib.Tests
{
    [TestClass]
    public class DMSTests
    {
        [DataTestMethod]
        [DataRow(33.3, 44.4, "33d18'00.0\"N 044d24'00.0\"E")]
        [DataRow(0.0, 123.456, "00d00'00.0\"N 123d27'21.6\"E")]
        public void TestEncode(double lat, double lon, string dms)
        {
            var latStr = DMS.Encode(lat, 5, HemisphereIndicator.Latitude);
            var lonStr = DMS.Encode(lon, 5, HemisphereIndicator.Longitude);
            Assert.AreEqual(dms, $"{latStr} {lonStr}");
        }

        [DataTestMethod]
        [DataRow(33.3, 44.4, "33d18'00.0\"N 044d24'00.0\"E")]
        [DataRow(0.0, 123.456, "00d00'00.0\"N 123d27'21.6\"E")]
        public void TestDecode(double lat, double lon, string dms)
        {
            var components = dms.Split(' ');
            var (latVal, lonVal) = DMS.Decode(components[0], components[1]);
            Assert.AreEqual(lat, latVal, 0.001);
            Assert.AreEqual(lon, lonVal, 0.001);
        }
    }
}
