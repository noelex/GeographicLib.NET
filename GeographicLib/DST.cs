using GeographicLib.Dsp;
using System;
using System.Numerics;
using static System.Math;
using static System.Numerics.Complex;

namespace GeographicLib
{
    public class DST
    {
        private readonly Fft64 _fft;
        private readonly int _N;

        /// <summary>
        /// Constructor specifying the number of points to use.
        /// </summary>
        /// <param name="N">the number of points to use.</param>
        public DST(int N = 0)
        {
            _N = N;
            _fft = new Fft64(2 * N);
        }

        /// <summary>
        /// Return the number of points.
        /// </summary>
        public int N => _N;

        /// <summary>
        /// Determine first <i>N</i> terms in the Fourier series.
        /// </summary>
        /// <param name="f"></param>
        /// <param name="F"></param>
        public void Transform(Func<double, double> f, Span<double> F)
        {
            Span<double> data = stackalloc double[4 * _N];
            var d = PI / (2 * _N);
            for (int i = 1; i <= _N; ++i)
                data[i] = f(i * d);
            fft_transform(data, F, false);
        }

        /// <summary>
        /// Refine the Fourier series by doubling the number of points sampled
        /// </summary>
        /// <param name="f"></param>
        /// <param name="F"></param>
        public void Refine(Func<double, double> f, Span<double> F)
        {
            Span<double> data = stackalloc double[4 * _N];
            var d = PI / (4 * _N);
            for (int i = 0; i < _N; ++i)
                data[i] = f((2 * i + 1) * d);
            fft_transform2(data, F);
        }

        /// <summary>
        /// Evaluate the Fourier sum given the sine and cosine of the angle.
        /// </summary>
        /// <param name="sinx"></param>
        /// <param name="cosx"></param>
        /// <param name="F"></param>
        /// <param name="N"></param>
        /// <returns></returns>
        public double Eval(double sinx, double cosx, ReadOnlySpan<double> F, int N)
        {
            // Evaluate
            // y = sum(F[i] * sin((2*i+1) * x), i, 0, N-1)
            // using Clenshaw summation.
            // Approx operation count = (N + 5) mult and (2 * N + 2) add
            double
              ar = 2 * (cosx - sinx) * (cosx + sinx), // 2 * cos(2 * x)
              y0 = (N & 1) != 0 ? F[--N] : 0, y1 = 0;          // accumulators for sum
                                                               // Now N is even
            while (N > 0)
            {
                // Unroll loop x 2, so accumulators return to their original role
                y1 = ar * y0 - y1 + F[--N];
                y0 = ar * y1 - y0 + F[--N];
            }
            return sinx * (y0 + y1);    // sin(x) * (y0 + y1)
        }

        /// <summary>
        /// Evaluate the integral of Fourier sum given the sine and cosine of the
        /// angle.
        /// </summary>
        /// <param name="sinx"></param>
        /// <param name="cosx"></param>
        /// <param name="F"></param>
        /// <param name="N"></param>
        /// <returns></returns>
        public double Integral(double sinx, double cosx, ReadOnlySpan<double> F, int N)
        {
            // Evaluate
            // y = -sum(F[i]/(2*i+1) * cos((2*i+1) * x), i, 0, N-1)
            // using Clenshaw summation.
            // Approx operation count = (N + 5) mult and (2 * N + 2) add
            int l = N;
            double
              ar = 2 * (cosx - sinx) * (cosx + sinx), // 2 * cos(2 * x)
              y0 = (N & 1) != 0 ? F[--N] / (2 * (--l) + 1) : 0, y1 = 0; // accumulators for sum
                                                                        // Now N is even
            while (N > 0)
            {
                // Unroll loop x 2, so accumulators return to their original role
                y1 = ar * y0 - y1 + F[--N] / (2 * (--l) + 1);
                y0 = ar * y1 - y0 + F[--N] / (2 * (--l) + 1);
            }
            return cosx * (y1 - y0);    // cos(x) * (y1 - y0)
        }

        private void fft_transform(Span<double> data, Span<double> F, bool centerp)
        {
            // Implement DST-III (centerp = false) or DST-IV (centerp = true).

            // Elements (0,N], resp. [0,N), of data should be set on input for centerp
            // = false, resp. true.  F must have a size of at least N and on output
            // elements [0,N) of F contain the transform.
            if (_N == 0) return;
            if (centerp)
            {
                for (int i = 0; i < _N; ++i)
                {
                    data[_N + i] = data[_N - 1 - i];
                    data[2 * _N + i] = -data[i];
                    data[3 * _N + i] = -data[_N - 1 - i];
                }
            }
            else
            {
                data[0] = 0;            // set [0]
                for (int i = 1; i < _N; ++i) data[_N + i] = data[_N - i]; // set [N+1,2*N-1]
                for (int i = 0; i < 2 * _N; ++i) data[2 * _N + i] = -data[i]; // [2*N, 4*N-1]
            }

            Span<double>
                reals = stackalloc double[2 * _N],
                imags = stackalloc double[2 * _N];

            var ctemp = new ComplexSpan(reals, imags);
            _fft.DFT(data, ctemp);

            if (centerp)
            {
                var d = -PI / (4 * _N);
                for (int i = 0, j = 1; i < _N; ++i, j += 2)
                    ctemp[j] *= Exp(new Complex(0, j * d));
            }
            for (int i = 0, j = 1; i < _N; ++i, j += 2)
            {
                F[i] = -ctemp[j].Imaginary / (2 * _N);
            }
        }

        private void fft_transform2(Span<double> data, Span<double> F)
        {
            // Elements [0,N), of data should be set to the N grid center values and F
            // should have size of at least 2*N.  On input elements [0,N) of F contain
            // the size N transform; on output elements [0,2*N) of F contain the size
            // 2*N transform.
            fft_transform(data, F.Slice(_N), true);
            // Copy DST-IV order N tx to [0,N) elements of data
            for (int i = 0; i < _N; ++i) data[i] = F[i + _N];
            for (int i = _N; i < 2 * _N; ++i)
                // (DST-IV order N - DST-III order N) / 2
                F[i] = (data[2 * _N - 1 - i] - F[2 * _N - 1 - i]) / 2;
            for (int i = 0; i < _N; ++i)
                // (DST-IV order N + DST-III order N) / 2
                F[i] = (data[i] + F[i]) / 2;
        }

        private void TransformScalar()
        {

        }
    }
}