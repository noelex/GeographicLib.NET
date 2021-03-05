using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib.Tests
{
    [TestClass]
    public class SetupAssemblyInitializer
    {
        [AssemblyInitialize]
        public static void AssemblyInit(TestContext context)
        {
            // Initalization code goes here
            if(context.Properties.Contains("cmath") && (string)context.Properties["cmath"] == "native")
            {
                MathEx.UseManagedCMath = false;
            }
            else
            {
                MathEx.UseManagedCMath = true;
            }

            context.WriteLine($"Executing tests with {(MathEx.UseManagedCMath ? "managed" : "native")} C mathematical functions.");
        }
    }
}
