using System;
using System.Numerics;

namespace GeographicLib.Dsp
{
    internal partial class Fft64
    {
        private readonly int windowSize;
        private readonly int depth;
        private readonly double[] dftCoefficientReal;
        private readonly double[] dftCoefficientImaginary;
        private readonly double[] idftCoefficientReal;
        private readonly double[] idftCoefficientImaginary;

        public Fft64(int windowSize)
        {
            if (windowSize <= 0)
            {
                throw new ArgumentException("n must be a positive integer.");
            }
            if ((windowSize & (windowSize - 1)) != 0)
            {
                throw new ArgumentException("n must be a power of 2.");
            }
            this.windowSize = windowSize;
            depth = Log2(windowSize);

            // Initialize coefficient
            Complex angleInterval = -2 * Complex.ImaginaryOne * Math.PI / windowSize;
            dftCoefficientReal = new double[windowSize];
            dftCoefficientImaginary = new double[windowSize];
            for (int i = 0; i < windowSize; i++)
            {
                Complex coefficient = Complex.Exp(i * angleInterval);
                dftCoefficientReal[i] = coefficient.Real;
                dftCoefficientImaginary[i] = coefficient.Imaginary;
            }

            idftCoefficientReal = new double[windowSize];
            idftCoefficientImaginary = new double[windowSize];
            for (int i = 0; i < windowSize; i++)
            {
                Complex coefficient = 1 / Complex.Exp(i * angleInterval);
                idftCoefficientReal[i] = coefficient.Real;
                idftCoefficientImaginary[i] = coefficient.Imaginary;
            }
        }

        public int WindowSize { get => windowSize; }

        private static int Log2(int n)
        {
            int log = 0;
            while (n > 1)
            {
                n = n >> 1;
                log++;
            }
            return log;
        }
    }
}