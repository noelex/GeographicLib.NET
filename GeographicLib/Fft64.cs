/*
 * This file is translated from KISS FFT - https://github.com/mborgerding/kissfft
 * 
 * SPDX-FileCopyrightText: 2003-2010, Mark Borgerding
 * SPDX-License-Identifier: BSD-3-Clause
 */

using System;
using System.Collections.Generic;
using System.Numerics;
using System.Runtime.InteropServices;
using static System.Math;
using static System.Numerics.Complex;

namespace GeographicLib
{
    internal class Fft64
    {
        private readonly int _nfft;
        private readonly bool _inverse;

        private readonly ReadOnlyMemory<Complex> _twiddles;
        private readonly List<int> _stageRadix = new List<int>();
        private readonly List<int> _stageRemainder = new List<int>();

        private Complex[] _scratchbuf = Array.Empty<Complex>();

        public Fft64(int nfft, bool inverse)
        {
            _nfft = nfft;
            _inverse = inverse;

            if (_nfft == 0) return;
            // fill twiddle factors
            var twiddles = new Complex[_nfft];
            {
                double s = _inverse ? 1 : -1;
                double d = Acos(-1) / (2 * _nfft);
                int i = 0, N = _nfft; // signed ints needed for subtractions
                                      // enforce trigonometric symmetries by evaluating sin and cos
                                      // with arguments in the range [-pi/4, pi/4]
                for (; 8 * i < N; ++i)
                    twiddles[i] = new Complex(+Cos((4 * i) * d),
                                         +s * Sin((4 * i) * d)); // pi/4*[0, 1)
                for (; 8 * i < 3 * N; ++i)
                    twiddles[i] = new Complex(-Sin((4 * i - N) * d),
                                         +s * Cos((4 * i - N) * d)); // pi/4*[1, 3)
                for (; 8 * i < 5 * N; ++i)
                    twiddles[i] = new Complex(-Cos((4 * i - 2 * N) * d),
                                         -s * Sin((4 * i - 2 * N) * d)); // pi/4*[3, 5)
                for (; 8 * i < 7 * N; ++i)
                    twiddles[i] = new Complex(+Sin((4 * i - 3 * N) * d),
                                         -s * Cos((4 * i - 3 * N) * d)); // pi/4*[5, 7)
                for (; i < N; ++i)
                    twiddles[i] = new Complex(+Cos((4 * i - 4 * N) * d),
                                         +s * Sin((4 * i - 4 * N) * d)); // pi/4*[5, 8)

                _twiddles = twiddles;
            }
            //factorize
            //start factoring out 4's, then 2's, then 3,5,7,9,...
            var n = _nfft;
            var p = 4;
            do
            {
                while ((n % p) != 0)
                {
                    switch (p)
                    {
                        case 4: p = 2; break;
                        case 2: p = 3; break;
                        default: p += 2; break;
                    }
                    if (p * p > n)
                        p = n;// no more factors
                }
                n /= p;
                _stageRadix.Add(p);
                _stageRemainder.Add(n);
            } while (n > 1);
        }

        public void Transform(ReadOnlySpan<Complex> fft_in, Span<Complex> fft_out, int stage = 0, int fstride = 1, int in_stride = 1)
        {
            if (_nfft == 0) return;
            var p = _stageRadix[stage];
            var m = _stageRemainder[stage];

            int fft_out_i = 0,
                fft_in_i = 0,
                Fout_beg = fft_out_i,
                Fout_end = fft_out_i + p * m;

            if (m == 1)
            {
                do
                {
                    fft_out[fft_out_i] = fft_in[fft_in_i];
                    fft_in_i += fstride * in_stride;
                } while (++fft_out_i != Fout_end);
            }
            else
            {
                do
                {
                    // recursive call:
                    // DFT of size m*p performed by doing
                    // p instances of smaller DFTs of size m,
                    // each one takes a decimated version of the input
                    Transform(fft_in.Slice(fft_in_i), fft_out.Slice(fft_out_i), stage + 1, fstride * p, in_stride);
                    fft_in_i += fstride * in_stride;
                } while ((fft_out_i += m) != Fout_end);
            }

            fft_out_i = Fout_beg;

            // recombine the p smaller DFTs
            switch (p)
            {
                case 2: kf_bfly2(fft_out, fstride, m); break;
                case 3: kf_bfly3(fft_out, fstride, m); break;
                case 4: kf_bfly4(fft_out, fstride, m); break;
                case 5: kf_bfly5(fft_out, fstride, m); break;
                default: kf_bfly_generic(fft_out, fstride, m, p); break;
            }
        }

        public void Transform(ReadOnlySpan<double> src, Span<Complex> dst)
        {
            int N = _nfft;
            if (N == 0)
                return;

            // perform complex FFT
            Transform(MemoryMarshal.Cast<double, Complex>(src), dst);

            // post processing for k = 0 and k = N
            dst[0] = new Complex(dst[0].Real + dst[0].Imaginary,
                               dst[0].Real - dst[0].Imaginary);

            // post processing for all the other k = 1, 2, ..., N-1
            double pi = Acos(-1);
            double half_phi_inc = (_inverse ? pi : -pi) / N;
            var twiddle_mul = Exp(new Complex(0, half_phi_inc));
            for (var k = 1; 2 * k < N; ++k)
            {
                var w = 0.5 * new Complex(
                     dst[k].Real + dst[N - k].Real,
                     dst[k].Imaginary - dst[N - k].Imaginary);
                var z = 0.5 * new Complex(
                     dst[k].Imaginary + dst[N - k].Imaginary,
                    -dst[k].Real + dst[N - k].Real);
                var twiddle =
                    k % 2 == 0 ?
                    _twiddles.Span[k / 2] :
                    _twiddles.Span[k / 2] * twiddle_mul;
                dst[k] = w + twiddle * z;
                dst[N - k] = Conjugate(w - twiddle * z);
            }
            if (N % 2 == 0)
                dst[N / 2] = Conjugate(dst[N / 2]);
        }

        void kf_bfly2(Span<Complex> Fout, int fstride, int m)
        {
            for (var k = 0; k < m; ++k)
            {
                var t = Fout[m + k] * _twiddles.Span[k * fstride];
                Fout[m + k] = Fout[k] - t;
                Fout[k] += t;
            }
        }

        void kf_bfly3(Span<Complex> Fout, int fstride, int m)
        {
            int k = m;
            int m2 = 2 * m;
            int tw1, tw2, Fout_i = 0;

            Span<Complex> scratch = stackalloc Complex[5];
            var epi3 = _twiddles.Span[fstride * m];

            tw1 = tw2 = 0;

            do
            {
                scratch[1] = Fout[m] * _twiddles.Span[tw1];
                scratch[2] = Fout[m2] * _twiddles.Span[tw2];

                scratch[3] = scratch[1] + scratch[2];
                scratch[0] = scratch[1] - scratch[2];
                tw1 += fstride;
                tw2 += fstride * 2;

                Fout[Fout_i + m] = Fout[0] - scratch[3] * 0.5;
                scratch[0] *= epi3.Imaginary;

                Fout[Fout_i + 0] += scratch[3];

                Fout[Fout_i + m2] = new Complex(Fout[m].Real + scratch[0].Imaginary, Fout[m].Imaginary - scratch[0].Real);

                Fout[Fout_i + m] += new Complex(-scratch[0].Imaginary, scratch[0].Real);
                ++Fout_i;
            } while ((--k) != 0);
        }

        void kf_bfly4(Span<Complex> Fout, int fstride, int m)
        {
            Span<Complex> scratch = stackalloc Complex[7];
            double negative_if_inverse = _inverse ? -1 : +1;
            for (var k = 0; k < m; ++k)
            {
                scratch[0] = Fout[k + m] * _twiddles.Span[k * fstride];
                scratch[1] = Fout[k + 2 * m] * _twiddles.Span[k * fstride * 2];
                scratch[2] = Fout[k + 3 * m] * _twiddles.Span[k * fstride * 3];
                scratch[5] = Fout[k] - scratch[1];

                Fout[k] += scratch[1];
                scratch[3] = scratch[0] + scratch[2];
                scratch[4] = scratch[0] - scratch[2];
                scratch[4] = new Complex(scratch[4].Imaginary * negative_if_inverse,
                                      -scratch[4].Real * negative_if_inverse);

                Fout[k + 2 * m] = Fout[k] - scratch[3];
                Fout[k] += scratch[3];
                Fout[k + m] = scratch[5] + scratch[4];
                Fout[k + 3 * m] = scratch[5] - scratch[4];
            }
        }

        void kf_bfly5(Span<Complex> Fout, int fstride, int m)
        {
            int Fout0, Fout1, Fout2, Fout3, Fout4;
            Span<Complex> scratch = stackalloc Complex[13];
            var ya = _twiddles.Span[fstride * m];
            var yb = _twiddles.Span[fstride * 2 * m];

            Fout0 = 0;
            Fout1 = Fout0 + m;
            Fout2 = Fout0 + 2 * m;
            Fout3 = Fout0 + 3 * m;
            Fout4 = Fout0 + 4 * m;

            for (var u = 0; u < m; ++u)
            {
                scratch[0] = Fout[Fout0];

                scratch[1] = Fout[Fout1] * _twiddles.Span[u * fstride];
                scratch[2] = Fout[Fout2] * _twiddles.Span[2 * u * fstride];
                scratch[3] = Fout[Fout3] * _twiddles.Span[3 * u * fstride];
                scratch[4] = Fout[Fout4] * _twiddles.Span[4 * u * fstride];

                scratch[7] = scratch[1] + scratch[4];
                scratch[10] = scratch[1] - scratch[4];
                scratch[8] = scratch[2] + scratch[3];
                scratch[9] = scratch[2] - scratch[3];

                Fout[Fout0] += scratch[7];
                Fout[Fout0] += scratch[8];

                scratch[5] = scratch[0] + new Complex(
                        scratch[7].Real * ya.Real + scratch[8].Real * yb.Real,
                        scratch[7].Imaginary * ya.Real + scratch[8].Imaginary * yb.Real
                        );

                scratch[6] = new Complex(
                         scratch[10].Imaginary * ya.Imaginary + scratch[9].Imaginary * yb.Imaginary,
                                -scratch[10].Real * ya.Imaginary - scratch[9].Real * yb.Imaginary

                        );

                Fout[Fout1] = scratch[5] - scratch[6];
                Fout[Fout4] = scratch[5] + scratch[6];

                scratch[11] = scratch[0] +
                    new Complex(
                            scratch[7].Real * yb.Real + scratch[8].Real * ya.Real,
                            scratch[7].Imaginary * yb.Real + scratch[8].Imaginary * ya.Real
                            );

                scratch[12] = new Complex(
                                -scratch[10].Imaginary * yb.Imaginary + scratch[9].Imaginary * ya.Imaginary,
                         scratch[10].Real * yb.Imaginary - scratch[9].Real * ya.Imaginary

                        );

                Fout[Fout2] = scratch[11] + scratch[12];
                Fout[Fout3] = scratch[11] - scratch[12];

                ++Fout0;
                ++Fout1;
                ++Fout2;
                ++Fout3;
                ++Fout4;
            }
        }

        /* perform the butterfly for one stage of a mixed radix FFT */
        void kf_bfly_generic(
                Span<Complex> Fout,
                int fstride,
                int m,
                int p
                )
        {
            if (p > _scratchbuf.Length) Array.Resize(ref _scratchbuf, p);

            for (var u = 0; u < m; ++u)
            {
                int k = u;
                for (var q1 = 0; q1 < p; ++q1)
                {
                    _scratchbuf[q1] = Fout[k];
                    k += m;
                }

                k = u;
                for (var q1 = 0; q1 < p; ++q1)
                {
                    int twidx = 0;
                    Fout[k] = _scratchbuf[0];
                    for (var q = 1; q < p; ++q)
                    {
                        twidx += fstride * k;
                        if (twidx >= _nfft)
                            twidx -= _nfft;
                        Fout[k] += _scratchbuf[q] * _twiddles.Span[twidx];
                    }
                    k += m;
                }
            }
        }
    }
}