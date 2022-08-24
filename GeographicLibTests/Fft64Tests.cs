using GeographicLib.Dsp;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Numerics;

namespace GeographicLib.Tests
{
    [TestClass]
    public class Fft64Tests
    {
        private static double Tolerance = 1e-7;
        private static ComplexComparer complexComparer = new ComplexComparer();

        public static DftTheoryTestDataProvider TestData { get; } = new DftTheoryTestDataProvider();

        [DataTestMethod]
        [DynamicData("TestData", typeof(Fft64Tests))]
        public void ArrayDftTheory(Complex[] signalExpected, Complex[] spectrumExpected)
        {
            int windowSize = signalExpected.Length;
            var calculator = new Fft64(windowSize);

            Complex[] spectrumActual = new Complex[windowSize];
            calculator.DFT(signalExpected, spectrumActual);
            Assert.That.AreEqual(spectrumExpected, spectrumActual, complexComparer);

            Complex[] signalActual = new Complex[windowSize];
            calculator.IDFT(spectrumActual, signalActual);
            Assert.That.AreEqual(signalExpected, signalActual, complexComparer);
        }

        [DataTestMethod]
        [DynamicData("TestData", typeof(Fft64Tests))]
        public void SpanDftTheory(Complex[] signalExpected, Complex[] spectrumExpected)
        {
            int windowSize = signalExpected.Length;
            var calculator = new Fft64(windowSize);

            Span<double> signalExpectedReal = stackalloc double[windowSize];
            Span<double> signalExpectedImaginary = stackalloc double[windowSize];
            ComplexSpan signalExpectedSpan = new ComplexSpan(signalExpectedReal, signalExpectedImaginary);
            CopyTo(signalExpected,signalExpectedSpan);
            Span<double> spectrumExpectedReal = stackalloc double[windowSize];
            Span<double> spectrumExpectedImaginary = stackalloc double[windowSize];
            ComplexSpan spectrumExpectedSpan = new ComplexSpan(spectrumExpectedReal, spectrumExpectedImaginary);
            CopyTo(spectrumExpected, spectrumExpectedSpan);

            Span<double> signalActualReal = stackalloc double[windowSize];
            Span<double> signalActualImaginary = stackalloc double[windowSize];
            ComplexSpan signalActualSpan = new ComplexSpan(signalActualReal, signalActualImaginary);
            Span<double> spectrumActualReal = stackalloc double[windowSize];
            Span<double> spectrumActualImaginary = stackalloc double[windowSize];
            ComplexSpan spectrumActualSpan = new ComplexSpan(spectrumActualReal, spectrumActualImaginary);

            calculator.DFT(signalExpectedSpan, spectrumActualSpan);
            Assert.That.AreEqual(spectrumExpected, ToArray(spectrumActualSpan), complexComparer);

            calculator.IDFT(spectrumActualSpan, signalActualSpan);
            Assert.That.AreEqual(signalExpected, ToArray(signalActualSpan), complexComparer);
        }

        private static void CopyTo(Complex[] array, ComplexSpan span)
        {
            Span<double> realSpan = span.Real;
            Span<double> imaginarySpan = span.Imaginary;
            for (int i = 0; i < realSpan.Length; i++)
            {
                realSpan[i] = array[i].Real;
                imaginarySpan[i] = array[i].Imaginary;
            }
        }

        private static Complex[] ToArray(ComplexSpan span)
        {
            Complex[] array = new Complex[span.Real.Length];
            CopyTo(span, array);
            return array;
        }

        private static void CopyTo(ComplexSpan span, Complex[] array)
        {
            Span<double> realSpan = span.Real;
            Span<double> imaginarySpan = span.Imaginary;
            for (int i = 0; i < realSpan.Length; i++)
            {
                array[i] = new Complex(realSpan[i], imaginarySpan[i]);
            }
        }

        private class ComplexComparer : IEqualityComparer<Complex>
        {
            public bool Equals(Complex x, Complex y)
            {
                return (x - y).Magnitude / (x + y).Magnitude < Tolerance;
            }

            public int GetHashCode(Complex obj)
            {
                throw new NotImplementedException();
            }
        }

        public class DftTheoryTestDataProvider : IEnumerable<object[]>
        {
            public IEnumerator<object[]> GetEnumerator()
            {
                return GenerateTestData().GetEnumerator();
            }

            IEnumerator IEnumerable.GetEnumerator()
            {
                return GetEnumerator();
            }

            private IEnumerable<object[]> GenerateTestData()
            {
                yield return SingleDataPoint;
                yield return DualDataPoint;
                yield return QuadDataPoint;
                yield return new object[] { FftTestData4096.Signal, FftTestData4096.Spectrum };
            }

            private object[] SingleDataPoint
            {
                get
                {
                    Complex[] data = new Complex[]
                    {
                        new Complex(1, 1)
                    };
                    return new object[] { data, data };
                }
            }

            private object[] DualDataPoint
            {
                get
                {
                    Complex[] signal = new Complex[]
                    {
                        new Complex(1.7, 3.2),
                        new Complex(-2.3, 1.6)
                    };

                    Complex[] spectrum = new Complex[]
                    {
                        new Complex(-0.6, 4.8),
                        new Complex(4, 1.6)
                    };

                    return new object[] { signal, spectrum };
                }
            }

            private object[] QuadDataPoint
            {
                get
                {
                    Complex[] signal = new Complex[]
                    {
                        new Complex(0.3, -2.7),
                        new Complex(-2.2, -1.8),
                        new Complex(3.1, -2.5),
                        new Complex(3.1, 0.7)
                    };

                    Complex[] spectrum = new Complex[]
                    {
                        new Complex(4.3, -6.3),
                        new Complex(-5.3, 5.1),
                        new Complex(2.5, -4.1),
                        new Complex(-0.3, -5.5)
                    };

                    return new object[] { signal, spectrum };
                }
            }
        }
    }
}
