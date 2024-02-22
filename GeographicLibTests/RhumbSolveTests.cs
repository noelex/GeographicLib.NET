using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Data;

namespace GeographicLib.Tests
{
    [TestClass]
    public class RhumbSolveTests
    {
        /// <summary>
        /// RhumbSovle0 - RhumbSolve1
        /// </summary>
        [DataTestMethod]
        [DataRow(false, 0, 0, 90, 0, 10001965.729, 0, 0, DisplayName = "RhumbSovle0")]
        [DataRow(true, 0, 0, 90, 0, 10001965.729, 0, 0, DisplayName = "RhumbSovle1")]
        [DataRow(false, 0, 0, 90, 30, 10001965.729, 0, 21252734238504, DisplayName = "RhumbSovle2")]
        [DataRow(true, 0, 0, 90, 30, 10001965.729, 0, 21252734238504, DisplayName = "RhumbSolve3")]
        public void PathsBetweenEquatorAndPole(
            bool exact,
            double lat1, double lon1, double lat2, double lon2,
            double s12, double azi12, double S12)
        {
            var r = new Rhumb(Ellipsoid.WGS84, exact);

            r.Inverse(lat1, lon1, lat2, lon2, out var _s12, out var _azi12, out var _S12);
            Assert.AreEqual(s12, _s12, 0.001);
            Assert.AreEqual(azi12, _azi12, 0.01);
            Assert.AreEqual(S12, _S12, 1);
        }


        /// <summary>
        /// RhumbSovle4 - RhumbSolve7
        /// </summary>
        [DataTestMethod]
        [DataRow(false, 6.4e6, 0, 20, 0, 37.1, 8399000, 79.97173787686, 90.00285375451, 52220156106664, DisplayName = "RhumbSolve4")]
        [DataRow(true, 6.4e6, 0, 20, 0, 37.1, 8399000, 79.97173787686, 90.00285375451, 52220156106664, DisplayName = "RhumbSolve5")]
        [DataRow(true, 6.4e6, 0.5, 20, 0, 51.5, 8083000, 80.02647493551, 90.02807944939, 29256421226004, DisplayName = "RhumbSolve6")]
        [DataRow(true, 6.4e6, -1, 20, 0, 30.0, 8229000, 77.59323119420, 90.39005992280, 100655035400805, DisplayName = "RhumbSolve7")]
        public void HiglyProlateOblateSphere(
            bool exact,
            double a, double f,
            double lat1, double lon1, double azi12, double s12,
            double lat2, double lon2, double S12)
        {
            var r = new Rhumb(a, f, exact);
            r.Direct(lat1, lon1, azi12, s12, out var _lat2, out var _lon2, out var _S12);
            Assert.AreEqual(lat2, _lat2, 1e-11);
            Assert.AreEqual(lon2, _lon2, 1e-11);
            Assert.AreEqual(S12, _S12, 1);
        }

        /// <summary>
        /// RhumbSovle8 - RhumbSolve11
        /// </summary>
        [DataTestMethod]
        [DataRow(0, 0, 0, 10001965, 89.99999347, 0.0, 0, false, DisplayName = "RhumbSovle8")]
        [DataRow(0, 0, 60, 20003930, 89.99999347, 1654.69811, 1123576495779545, true, DisplayName = "RhumbSovle9")]
        [DataRow(0, 0, 0, 10001966, 89.99999758, double.NaN, double.NaN, false, DisplayName = "RhumbSovle10")]
        [DataRow(0, 0, 60, 20003932, 89.99999758, double.NaN, double.NaN, true, DisplayName = "RhumbSovle11")]
        public void RhumbLinesWithDestinationsCloseToPole(
            double lat1, double lon1, double azi12, double s12,
            double lat2, double lon2, double S12,
            bool unroll)
        {
            var r = Rhumb.WGS84;
            r.GenDirect(
                lat1, lon1, azi12, s12,
                !unroll ? GeodesicFlags.All : GeodesicFlags.All | GeodesicFlags.LongUnroll,
                out var _lat2, out var _lon2, out var _S12);
            Assert.AreEqual(lat2, _lat2, 1e-8);

            if (double.IsNaN(lon2))
            {
                Assert.IsTrue(double.IsNaN(_lon2));
            }
            else
            {
                Assert.AreEqual(lon2, _lon2, 1e-5);
            }

            if (double.IsNaN(S12))
            {
                Assert.IsTrue(double.IsNaN(_S12));
            }
            else
            {
                Assert.AreEqual(S12, _S12, 1);
            }
        }
    }
}
