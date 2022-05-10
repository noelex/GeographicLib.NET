using GeographicLib.Geocodes;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
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
        [DataRow("38SMB4484", "38n 444500 3684500", true, -3)]
        [DataRow("38SMB44148470", "38n 444145 3684705", true, -1)]
        [DataRow("38SMB", "38n 400000 3600000", false, -5)]
        [DataRow("38SMB4484", "38n 444000 3684000", false, -3)]
        [DataRow("38SMB44148470", "38n 444140 3684700", false, -1)]
        public void UTMUPSToMGRSTest(string mgrs, string utmups, bool centerp, int prec)
        {
            var coord = new GeoCoords(utmups, centerp);
            Assert.AreEqual(mgrs, coord.ToMGRSString(prec));
        }

        [TestMethod]
        public void GeoConvert0()
        {
            var coord = new GeoCoords("33.3 44.4");
            Assert.AreEqual("38SMB4484", coord.ToMGRSString(-3));
        }

        [TestMethod]
        public void GeoConvert1()
        {
            var coord = new GeoCoords("38smb");
            Assert.That.MatchesRegex("32d59'14\\.1\"N 044d27'53\\.4\"E", coord.ToDMSString());
        }

        [TestMethod]
        public void GeoConvert2()
        {
            var coord = new GeoCoords("30d30'30\" 30.50833");
            Assert.AreEqual("30.508 30.508", coord.ToGeoString(-2));
        }

        [TestMethod]
        [ExpectedException(typeof(GeographicException))]
        public void GeoConvert3_4()
        {
            new GeoCoords("garbage");
        }

        [TestMethod]
        [ExpectedException(typeof(GeographicException))]
        public void GeoConvert5_CheckFixForDMSDecodeBug()
        {
            new GeoCoords("5d. 0");
        }

        /// <summary>
        /// Check fix for DMS::Decode double rounding bug fixed on 2012-11-15
        /// This test is known to fail for VC 11 and 12 bug reported 2013-01-10
        /// OK to skip this test for these compilers because this is a question
        /// of accuracy of the least significant bit.  The bug is fixed in VC 14.
        /// 
        /// N.B. 179.99999999999998578 = 180 - 0.50032 * 0.5^45 which (as a
        /// double) rounds to 180 - 0.5^45 = 179.9999999999999716
        /// </summary>
        [TestMethod]
        public void GeoConvert6_CheckFixForDMSDevoceDoubleRoundingBug()
        {
            var coord = new GeoCoords("0 179.99999999999998578");

            // For .NET implementation:
            // This test fails on .NET Core 2.1, as 179.99999999999998578 is rounded up to 180.
#if NETCOREAPP2_1
            Assert.That.MatchesRegex("180", coord.ToGeoString(9));
#else
            Assert.That.MatchesRegex("179\\.9999999999999[7-9]", coord.ToGeoString(9));
#endif
        }

        [TestMethod]
        public void GeoConvert7_UTMUPSSantityCheck()
        {
            MGRS.Check();
        }

        [TestMethod]
        public void GeoConvert8_CheckFixToPolarStereographicInitializationBlunder()
        {
            var coord = new GeoCoords("86 0");
            Assert.That.MatchesRegex("n 2000000\\.0* 1555731\\.570643", coord.ToUTMUPSString(6));
        }

        /// <summary>
        /// Check that integer(minutes) >= 60 and decimal(minutes) > 60 fail.
        /// Latter used to succeed; fixed 2015-06-11.
        /// </summary>
        /// <param name="s"></param>
        [DataTestMethod]
        [DataRow("5d70.0 10")]
        [DataRow("5d60 10")]
        [ExpectedException(typeof(GeographicException))]
        public void GeoConvert9_10(string s)
        {
            new GeoCoords(s);
        }

        /// <summary>
        /// Check that integer(minutes) < 60 and decimal(minutes) <= 60 succeed.
        /// Latter used to fail with 60.; fixed 2015-06-11.
        /// </summary>
        /// <param name="s"></param>
        [DataTestMethod]
        [DataRow("5d59 10")]
        [DataRow("5d60. 10")]
        [DataRow("5d60.0 10")]
        public void GeoConvert11_12_13(string s)
        {
            new GeoCoords(s);
            Assert.IsTrue(true);
        }

        [DataTestMethod]
        [DataRow("5.25 5.75", -4, '\0', "05.2N 005.8E")]
        [DataRow("5.03125 5.09375", -1, ':', "05:01:52N 005:05:38E")]
        public void GeoConvert14_15_CheckDMSEncodeDoesRoundTiesToEven(string s, int prec, char sep, string r)
        {
            var coord = new GeoCoords(s);
            Assert.AreEqual(r, coord.ToDMSString(prec, dmssep: sep));
        }

        [TestMethod]
        public void GeoConvert16_CheckMGRSForwardImprovedRoundingFix()
        {
            var coord = new GeoCoords("38n 444140.6 3684706.3");
            Assert.AreEqual("38SMB4414060084706300", coord.ToMGRSString(3));
        }

        [TestMethod]
        public void GeoConvert17_18_CheckMGRSForwardDigitConsistencyFix()
        {
            var coord = new GeoCoords("38n 500000 63.811");
            Assert.AreEqual("38NNF0000000000063811", coord.ToMGRSString(3));
            Assert.AreEqual("38NNF000000000000638110", coord.ToMGRSString(4));
        }

        [DataTestMethod]
        [DataRow("s 2746000 1515000", -6, "B")]
        [DataRow("s 2746000 1515000", -5, "BKH")]
        [DataRow("s 2746000 1515000", -4, "BKH41")]
        public void GeoConvert19_20_21_Check_prec_eq_minus_6_ForUPS(string s, int prec, string r)
        {
            var coord = new GeoCoords(s);
            Assert.AreEqual(r, coord.ToMGRSString(prec));
        }

        [TestMethod]
        public void GeoConvert22_23_CheckDMSEncodeRoundTiesToEventForWholeDegrees()
        {
            var coord = new GeoCoords("5.5 6.5");
            Assert.AreEqual("6 6", coord.ToGeoString(-5));
            Assert.AreEqual("06N 006E", coord.ToDMSString(-5));
        }
    }
}
