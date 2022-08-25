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

        public void DFT(Complex[] signal, Complex[] spectrum)
        {
            Span<double> signalReal = stackalloc double[windowSize];
            Span<double> signalImaginary = stackalloc double[windowSize];
            Span<double> spectrumReal = stackalloc double[windowSize];
            Span<double> spectrumImaginary = stackalloc double[windowSize];

            for (int i = 0; i < windowSize; i++)
            {
                signalReal[i] = signal[i].Real;
                signalImaginary[i] = signal[i].Imaginary;
            }

            DFT(new ComplexSpan(signalReal, signalImaginary), new ComplexSpan(spectrumReal, spectrumImaginary));

            for (int i = 0; i < windowSize; i++)
            {
                spectrum[i] = new Complex(spectrumReal[i], spectrumImaginary[i]);
            }
        }

        public void DFT(ReadOnlySpan<double> signal, ComplexSpan spectrum)
        {
            Span<double>
                reals = stackalloc double[signal.Length / 2],
                imags = stackalloc double[signal.Length / 2];

            var cpx = new ComplexSpan(reals, imags);
            cpx.CopyFromScalar(signal);

            FFT(cpx, spectrum, new ComplexSpan(dftCoefficientReal, dftCoefficientImaginary));
        }

        public void IDFT(Complex[] spectrum, Complex[] signal)
        {
            Span<double> signalReal = stackalloc double[windowSize];
            Span<double> signalImaginary = stackalloc double[windowSize];
            Span<double> spectrumReal = stackalloc double[windowSize];
            Span<double> spectrumImaginary = stackalloc double[windowSize];

            for (int i = 0; i < windowSize; i++)
            {
                spectrumReal[i] = spectrum[i].Real;
                spectrumImaginary[i] = spectrum[i].Imaginary;
            }

            IDFT(new ComplexSpan(spectrumReal, spectrumImaginary), new ComplexSpan(signalReal, signalImaginary));

            for (int i = 0; i < windowSize; i++)
            {
                signal[i] = new Complex(signalReal[i], signalImaginary[i]);
            }
        }

        public void DFT(ReadOnlyComplexSpan signal, ComplexSpan spectrum)
        {
            FFT(signal, spectrum, new ComplexSpan(dftCoefficientReal, dftCoefficientImaginary));
        }

        public void IDFT(ReadOnlyComplexSpan spectrum, ComplexSpan signal)
        {
            FFT(spectrum, signal, new ComplexSpan(idftCoefficientReal, idftCoefficientImaginary));
            double coefficientAdjustment = 1D / windowSize;
            Span<double> signalReal = signal.Real;
            Span<double> signalImaginary = signal.Imaginary;
            for (int i = 0; i < windowSize; i++)
            {
                signalReal[i] *= coefficientAdjustment;
                signalImaginary[i] *= coefficientAdjustment;
            }
        }

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