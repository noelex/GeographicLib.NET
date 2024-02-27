using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.IO;
using System.Reflection;

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

                Assert.IsNotNull(model.GravityFile);
                Assert.IsNotNull(model.GravityModelDirectory);
            }
        }

        [TestMethod]
        public void Test_LoadModelFromStream()
        {
            var metadataFile = Path.Combine(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location), "gravity", "egm96.egm");
            var coeffFile = metadataFile + ".cof";

            var metadataStream = File.OpenRead(metadataFile);
            var coeffStream = File.OpenRead(coeffFile);
            var model = new GravityModel(metadataStream, coeffStream);

            try
            {
                _ = metadataStream.Length;
                Assert.Fail("Expects ObjectDisposedException to be thrown.");
            }
            catch (ObjectDisposedException)
            {

            }

            try
            {
                _ = coeffStream.Length;
                Assert.Fail("Expects ObjectDisposedException to be thrown.");
            }
            catch (ObjectDisposedException)
            {

            }

            Assert.AreEqual("Earth Gravity Model 1996", model.Description);
            Assert.AreEqual(DateTime.Parse("1997-01-01"), model.DateTime);
            Assert.AreEqual("egm96", model.GravityModelName);
            Assert.AreEqual(3986004.418e8, model.MassConstant);
            Assert.AreEqual(7292115e-11, model.AngularVelocity);

            Assert.IsNull(model.GravityFile);
            Assert.IsNull(model.GravityModelDirectory);
        }

        [TestMethod]
        public void Test_LoadModelFromStream_LeaveOpen()
        {
            var metadataFile = Path.Combine(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location), "gravity", "egm96.egm");
            var coeffFile = metadataFile + ".cof";

            var metadataStream = File.OpenRead(metadataFile);
            var coeffStream = File.OpenRead(coeffFile);
            _ = new GravityModel(metadataStream, coeffStream, leaveOpen: true);

            _ = metadataStream.Length;
            _ = coeffStream.Length;
        }

        [TestMethod]
        public void Test_LoadModelFromByteArray()
        {
            var metadataFile = Path.Combine(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location), "gravity", "egm96.egm");
            var coeffFile = metadataFile + ".cof";

            var metadataBytes = File.ReadAllBytes(metadataFile);
            var coeffBytes = File.ReadAllBytes(coeffFile);
            var model = new GravityModel(metadataBytes, coeffBytes);

            Assert.AreEqual("Earth Gravity Model 1996", model.Description);
            Assert.AreEqual(DateTime.Parse("1997-01-01"), model.DateTime);
            Assert.AreEqual("egm96", model.GravityModelName);
            Assert.AreEqual(3986004.418e8, model.MassConstant);
            Assert.AreEqual(7292115e-11, model.AngularVelocity);

            Assert.IsNull(model.GravityFile);
            Assert.IsNull(model.GravityModelDirectory);
        }

        [TestMethod]
        public void Test_GravityCircle()
        {
            var model = new GravityModel("egm96",
                Path.Combine(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location), "gravity"));
            {
                var circle = model.Circle(0, 0, GravityFlags.All);

                var (V, GX, GY, GZ) = circle.V(0);
                Assert.IsFalse(double.IsNaN(V));
                Assert.IsFalse(double.IsNaN(GX));
                Assert.IsFalse(double.IsNaN(GY));
                Assert.IsFalse(double.IsNaN(GZ));

                var (T, dX, dY, dZ) = circle.Td(0);
                Assert.IsFalse(double.IsNaN(T));
                Assert.IsFalse(double.IsNaN(dX));
                Assert.IsFalse(double.IsNaN(dY));
                Assert.IsFalse(double.IsNaN(dZ));

                var h = circle.GeoidHeight(0);
                Assert.IsFalse(double.IsNaN(h));

                var (Dg01, xi, eta) = circle.SphericalAnomaly(0);
                Assert.IsFalse(double.IsNaN(Dg01));
                Assert.IsFalse(double.IsNaN(xi));
                Assert.IsFalse(double.IsNaN(eta));

                circle = model.Circle(0, 0, GravityFlags.None);

                (V, GX, GY, GZ) = circle.V(0);
                Assert.IsTrue(double.IsNaN(V));
                Assert.IsTrue(double.IsNaN(GX));
                Assert.IsTrue(double.IsNaN(GY));
                Assert.IsTrue(double.IsNaN(GZ));

                (T, dX, dY, dZ) = circle.Td(0);
                Assert.IsTrue(double.IsNaN(T));
                Assert.IsTrue(double.IsNaN(dX));
                Assert.IsTrue(double.IsNaN(dY));
                Assert.IsTrue(double.IsNaN(dZ));

                h = circle.GeoidHeight(0);
                Assert.IsTrue(double.IsNaN(h));

                (Dg01, xi, eta) = circle.SphericalAnomaly(0);
                Assert.IsTrue(double.IsNaN(Dg01));
                Assert.IsTrue(double.IsNaN(xi));
                Assert.IsTrue(double.IsNaN(eta));
            }
        }

        [DataTestMethod]
        [DataRow("27:59:17N 86:55:32E", 8820, -0.00021, 0.00084, -9.76653)]
        [DataRow("0 0", 0, -0.00002, 0.00001, -9.78037)]
        [DataRow("90 0", 0, -0.00007, -0.00006, -9.83208)]
        [DataRow("-90 0", 0, 0.00001, 0.00009, -9.83204)]
        public void Test_Evaluate(string coord, double height, double GX, double GY, double GZ)
        {
            var c = new GeoCoords(coord);
            var model = new GravityModel("egm96",
                Path.Combine(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location), "gravity"));
            {
                var (_, gx, gy, gz) = model.Gravity(c.Latitude, c.Longitude, height);
                Assert.AreEqual(GX, gx, 1e-5);
                Assert.AreEqual(GY, gy, 1e-5);
                Assert.AreEqual(GZ, gz, 1e-5);
            }
        }
    }
}
