using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using GeographicLib;
using System.Reflection;

namespace GeographicLib.Tests
{
    [TestClass]
    public class GeoidTests
    {
        class UnseekableMemoryStream : MemoryStream
        {
            public override bool CanSeek => false;
        }

        [TestMethod]
        public void DefaultGeoidPath_Test()
        {
            var expected = OperatingSystem.IsWindows() ? @"C:\ProgramData\GeographicLib\geoids" : @"/usr/local/share/GeographicLib/geoids";
            Assert.AreEqual(expected, Geoid.DefaultGeoidPath);
        }

        [TestMethod]
        public void DefaultGeoidName_Test()
        {
            Assert.AreEqual(@"egm96-5", Geoid.DefaultGeoidName);
        }

        private void CheckLoadedGeoid(Geoid geoid)
        {
            Assert.AreEqual("WGS84 EGM84, 30-minute grid", geoid.Description);
            Assert.AreEqual(DateTime.Parse("2009-08-29 18:45:02"), geoid.DateTime);
            Assert.AreEqual(0.274, geoid.MaxError);
            Assert.AreEqual(0.014, geoid.RMSError);
            Assert.AreEqual(-108, geoid.Offset);
            Assert.AreEqual(0.003, geoid.Scale);
        }

        [TestMethod]
        public void Test_Load()
        {
            using (var geoid = new Geoid("egm84-30",
                Path.Combine(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location), "geoids")))
            {
                CheckLoadedGeoid(geoid);
            }
        }

        [TestMethod]
        public void Test_LoadFromByteArray()
        {
            var buffer = File.ReadAllBytes(Path.Combine(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location), "geoids", "egm84-30.pgm"));
            using (var geoid = new Geoid(buffer))
            {
                Assert.AreEqual("WGS84 EGM84, 30-minute grid", geoid.Description);
                Assert.AreEqual(DateTime.Parse("2009-08-29 18:45:02"), geoid.DateTime);
                Assert.AreEqual(0.274, geoid.MaxError);
                Assert.AreEqual(0.014, geoid.RMSError);
                Assert.AreEqual(-108, geoid.Offset);
                Assert.AreEqual(0.003, geoid.Scale);
            }
        }

        [TestMethod]
        [ExpectedException(typeof(ObjectDisposedException))]
        public void Test_LoadFromStream()
        {
            var fs = File.OpenRead(Path.Combine(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location), "geoids", "egm84-30.pgm"));
            using (var geoid = new Geoid(fs))
            {
                CheckLoadedGeoid(geoid);
            }

            // This should throw as the stream is disposed by the geoid.
            fs.Seek(0, SeekOrigin.Begin);
        }

        [TestMethod]
        public void Test_LoadFromStreamAndLeaveOpen()
        {
            using (var fs = File.OpenRead(Path.Combine(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location), "geoids", "egm84-30.pgm")))
            {
                using (var geoid = new Geoid(fs, leaveOpen: true))
                {
                    CheckLoadedGeoid(geoid);
                }

                Assert.IsFalse(fs.SafeFileHandle.IsClosed);
            }
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Test_LoadFromUnseekableStream()
        {
            new Geoid(new UnseekableMemoryStream());
        }

        [DataTestMethod]
        [DataRow("531595 4468135 28n", 48.7602)]
        [DataRow("16:46:33N 3:00:34W", 31.3031)]
        public void Test_EvaluateFromCache(string c, double height)
        {
            using (var geoid = new Geoid("egm84-30",
                Path.Combine(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location), "geoids")))
            {
                geoid.CacheAll();

                var coord = new GeoCoords(c);
                var h = geoid.Evaluate(coord);

                Assert.AreEqual(height, h, 1e-4);
            }
        }

        [DataTestMethod]
        [DataRow("531595 4468135 28n", 48.7602)]
        [DataRow("16:46:33N 3:00:34W", 31.3031)]
        public void Test_EvaluateWithoutCache(string c, double height)
        {
            using (var geoid = new Geoid("egm84-30",
                Path.Combine(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location), "geoids")))
            {
                var coord = new GeoCoords(c);
                var h = geoid.Evaluate(coord);

                Assert.AreEqual(height, h, 1e-4);
            }
        }
    }
}
