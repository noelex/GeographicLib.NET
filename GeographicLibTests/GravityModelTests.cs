using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib.Tests
{
    [TestClass]
    public class GravityModelTests
    {
        [TestMethod]
        public void DefaultGravityPath_Test()
        {
            var expected = OperatingSystem.IsWindows() ? @"C:\ProgramData\GeographicLib\gravity" : @"/usr/local/share/GeographicLib/gravity";
            Assert.AreEqual(expected, GravityModel.DefaultGravityPath);
        }

        [TestMethod]
        public void DefaultGravityName_Test()
        {
            Assert.AreEqual(@"egm96", GravityModel.DefaultGravityName);
        }

        [TestMethod]
        public void Test_LoadModel()
        {
            var model = new GravityModel("egm96",
                Path.Combine(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location), "gravity"));
            {
                Assert.AreEqual("Earth Gravity Model 1996", model.Description);
                Assert.AreEqual(DateTime.Parse("1997-01-01"), model.DateTime);
                Assert.AreEqual("egm96", model.GravityModelName);
                Assert.AreEqual(3986004.418e8, model.MassConstant);
                Assert.AreEqual(7292115e-11, model.AngularVelocity);
            }
        }

        [DataTestMethod]
        [DataRow("27:59:17N 86:55:32E", 8820, -0.00021, 0.00084, -9.76653)]
        [DataRow("0 0", 0, -0.00002, 0.00001, - 9.78037)]
        [DataRow("90 0", 0, -0.00007, -0.00006, -9.83208)]
        [DataRow("-90 0", 0, 0.00001, 0.00009, - 9.83204)]
        public void Test_Evaluate(string coord, double height, double GX, double GY, double GZ)
        {
            var c = new GeoCoords(coord);
            var model = new GravityModel("egm96",
                Path.Combine(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location), "gravity"));
            {
                model.Gravity(c.Latitude, c.Longitude, height, out var gx, out var gy, out var gz);
                Assert.AreEqual(GX, gx, 1e-5);
                Assert.AreEqual(GY, gy, 1e-5);
                Assert.AreEqual(GZ, gz, 1e-5);
            }
        }
    }
}
