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
    public class MagneticModelTests
    {
        [TestMethod]
        public void DefaultMagneticPath_Test()
        {
            var expected = OperatingSystem.IsWindows() ? @"C:\ProgramData\GeographicLib\magnetic" : @"/usr/local/share/GeographicLib/magnetic";
            Assert.AreEqual(expected, MagneticModel.DefaultMagneticPath);
        }

        [TestMethod]
        public void DefaultMagneticName_Test()
        {
            Assert.AreEqual(@"wmm2020", MagneticModel.DefaultMagneticName);
        }

        [TestMethod]
        public void Test_LoadModel()
        {
            var model = new MagneticModel("wmm2020",
                Path.Combine(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location), "magnetic"));

            Assert.AreEqual("World Magnetic Model 2020", model.Description);
            Assert.AreEqual(DateTime.Parse("2019-12-10"), model.DateTime);
            Assert.AreEqual("wmm2020", model.MagneticModelName);
            Assert.AreEqual(2020, model.MinTime);
            Assert.AreEqual(2025, model.MaxTime);
            Assert.AreEqual(-1000, model.MinHeight);
            Assert.AreEqual(850000, model.MaxHeight);
        }

        [TestMethod]
        public void Test_LoadModel_EMM2017()
        {
            var model = new MagneticModel("emm2017",
                Path.Combine(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location), "magnetic"));

            Assert.AreEqual("Enhanced Magnetic Model 2017", model.Description);
            Assert.AreEqual(DateTime.Parse("2017-07-05"), model.DateTime);
            Assert.AreEqual("emm2017", model.MagneticModelName);
            Assert.AreEqual(2000, model.MinTime);
            Assert.AreEqual(2022, model.MaxTime);
            Assert.AreEqual(-20000, model.MinHeight);
            Assert.AreEqual(10000000, model.MaxHeight);
        }

        [DataTestMethod]
        [DataRow(2020, "16:46:33N 3:00:34W", 300,
            -1.59908987711, 12.00271322006, 33973.4874018165, 33960.2567138655, -948.0559950443, 7222.9691748258, 34732.8249634531,
            0.13240764812, -0.02108995839, 21.7574504999, 23.9398866911, 77.8732771740, -8.4447635876, 19.5256275320)]
        [DataRow(2020, "0 0", 0,
            -4.66834226021, -30.11087207520, 27628.0573899198, 27536.4015392405, -2248.5874255670, -16022.4296255844, 31937.8741660666,
            0.16276610551, -0.04413635078, -11.2706069741, -4.8454196603, 79.1428222375, -21.9041837809, 1.2390701594)]
        [DataRow(2020, "27:59:17N 86:55:32E", 8820,
            0.24105840589, 44.11617090945, 35030.5802595924, 35030.2702208329, 147.3824249174, 33966.1707417194, 48793.8757241047,
            -0.02065698417, 0.11143777174, 0.0081318874, 0.0612679529, -12.6295127346, 132.1962061103, 92.0296595883)]
        [DataRow(2020, "90 0", 100,
            3.65450069107, 88.15458612498, 1827.6747526472, 1823.9582665902, 116.4956788964, 56725.3976949931, 56754.8336149142,
            2.04651258030, 0.02476377883, -23.7553654223, -27.8680973065, 63.6346764477, 24.4366267718, 23.6589591548)]
        [DataRow(2020, "-90 0", 200,
            -30.73972682597, -72.12011164244, 16781.2992204543, 14423.4943407499, -8577.5762619038, -52018.3243612949, 54658.1930993269,
            -0.14428589981, 0.03527455492, 13.8318354690, -9.7121657134, -43.3921571332, 66.7277229518, -59.2582371089)]
        public void Test_Evaluate(double t, string c, double h,
            double D_, double I_, double H_, double By_, double Bx_, double Bz_, double F_,
            double Dt_, double It_, double Ht_, double Byt_, double Bxt_, double Bzt_, double Ft_)
        {
            var tolerance = 1e-10;
            var coord = new GeoCoords(c);

            var model = new MagneticModel("wmm2020",
                Path.Combine(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location), "magnetic"));

            var (Bx, By, Bz) = model.Evaluate(t, coord.Latitude, coord.Longitude, h, out var Bxt, out var Byt, out var Bzt);

            var (H, F, D, I) = MagneticModel.FieldComponents(Bx, By, Bz, Bxt, Byt, Bzt,
                                          out var Ht, out var Ft, out var Dt, out var It);

            Assert.AreEqual(D_, D, tolerance);
            Assert.AreEqual(I_, I, tolerance);
            Assert.AreEqual(H_, H, tolerance);
            Assert.AreEqual(By_, By, tolerance);
            Assert.AreEqual(Bx_, Bx, tolerance);
            Assert.AreEqual(Bz_, -Bz, tolerance);
            Assert.AreEqual(F_, F, tolerance);

            Assert.AreEqual(Dt_, Dt, tolerance);
            Assert.AreEqual(It_, It, tolerance);
            Assert.AreEqual(Ht_, Ht, tolerance);
            Assert.AreEqual(Byt_, Byt, tolerance);
            Assert.AreEqual(Bxt_, Bxt, tolerance);
            Assert.AreEqual(Bzt_, -Bzt, tolerance);
            Assert.AreEqual(Ft_, Ft, tolerance);
        }
    }
}
