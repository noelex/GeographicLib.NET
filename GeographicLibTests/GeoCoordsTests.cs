using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib.Tests
{
    [TestClass]
    public class GeoCoordsTests
    {
        [DataTestMethod]
        [DataRow("38SMB", "38n 450000 3650000", true)]
        [DataRow("38SMB4484", "38n 444500 3684500", true)]
        [DataRow("38SMB44148470", "38n 444145 3684705", true)]
        [DataRow("38SMB", "38n 400000 3600000", false)]
        [DataRow("38SMB4484", "38n 444000 3684000", false)]
        [DataRow("38SMB44148470", "38n 444140 3684700", false)]
        public void MGRSToUTMUPSTest(string mgrs, string utmups, bool centerp)
        {
            var coord = new GeoCoords(mgrs, centerp);
            Assert.AreEqual(utmups, coord.ToUTMUPSString());
        }

        [DataTestMethod]
        [DataRow("38SMB", "38n 450000 3650000", true, -5)]
        [DataRow("38SMB4484", "38n 444500 3684500", true,-3)]
        [DataRow("38SMB44148470", "38n 444145 3684705", true,-1)]
        [DataRow("38SMB", "38n 400000 3600000", false, -5)]
        [DataRow("38SMB4484", "38n 444000 3684000", false, -3)]
        [DataRow("38SMB44148470", "38n 444140 3684700", false, -1)]
        public void UTMUPSToMGRSTest(string mgrs, string utmups, bool centerp,int prec)
        {
            var coord = new GeoCoords(utmups,centerp);
            Assert.AreEqual(mgrs, coord.ToMGRSString(prec));
        }
    }
}
