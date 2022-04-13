using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib.Tests
{
    internal static class TestExtensions
    {
        public static void EqualsExactly(this Assert assert, double expected, double actual)
        {
            var equals = (double.IsNaN(expected) && double.IsNaN(actual)) ||
                 (expected == actual && MathEx.SignBit(expected) == MathEx.SignBit(actual));

            if (!equals)
            {
                var expectedStr = $"{(MathEx.SignBit(expected) ? "" : "+")}{expected}";
                var actualStr= $"{(MathEx.SignBit(actual) ? "" : "+")}{actual}";
                throw new AssertFailedException($"Expected: {expectedStr}, Actual: {actualStr}.");
            }
        }
        
    }
}

