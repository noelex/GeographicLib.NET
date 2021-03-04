using System;
using System.Collections.Generic;
using System.Text;

using static System.Math;
using static GeographicLib.MathEx;

namespace GeographicLib.SphericalHarmonics
{
    /// <summary>
    /// Spherical harmonic sums for a circle.
    /// </summary>
    /// <remarks>
    /// The class is a companion to <see cref="SphericalEngine"/>. 
    /// If the results of a spherical harmonic sum are needed for several points 
    /// on a circle of constant latitude <i>lat</i> and height <i>h</i>, then <see cref="SphericalEngine.Circle"/> can compute
    /// the inner sum, which is independent of longitude lon, and produce a <see cref="CircularEngine"/> object.
    /// <see cref="Evaluate(double)"/> can then be used to perform the outer sum for particular vales of <i>lon</i>.
    /// This can lead to substantial improvements in computational speed for high degree sum (approximately by a factor of <i>N</i> / 2
    /// where <i>N</i> is the maximum degree).
    /// <para>
    /// <see cref="CircularEngine"/> is tightly linked to the internals of <see cref="SphericalEngine"/>. 
    /// For that reason, the constructor for this class is private.
    /// Use <see cref="SphericalHarmonic.Circle"/>, <see cref="SphericalHarmonic1.Circle"/>, and <see cref="SphericalHarmonic2.Circle"/>
    /// to create instances of this class.
    /// </para>
    /// <para>
    /// <see cref="CircularEngine"/> stores the coefficients needed to allow the summation over order to be performed in 2 or 6 vectors of length <i>M</i> + 1
    /// (depending on whether gradients are to be calculated).
    /// </para>
    /// </remarks>
    public class CircularEngine
    {
        private readonly int _M;
        private readonly bool _gradp;
        private readonly Normalization _norm;
        private readonly double _a, _r, _u, _t;
        private readonly Memory<double> _wc, _ws, _wrc, _wrs, _wtc, _wts;
        private readonly double _q, _uq, _uq2;

        internal CircularEngine(int M, bool gradp, Normalization norm, double a, double r, double u, double t)
        {
            _M = M;
            _gradp = gradp;
            _norm = norm;
            _a = a;
            _r = r;
            _u = u;
            _t = t;

            _wc = new double[_M + 1];
            _ws = new double[_M + 1];
            _wrc = gradp ? new double[_M + 1] : null;
            _wrs = gradp ? new double[_M + 1] : null;
            _wtc = gradp ? new double[_M + 1] : null;
            _wts = gradp ? new double[_M + 1] : null;

            _q = _a / _r;
            _uq = _u * _q;
            _uq2 = Sq(_uq);
        }

        /// <summary>
        /// A default constructor. Calling <see cref="Evaluate(double)"/> and its overloads on the resulting object returns zero. 
        /// The resulting object can be assigned to the result of <see cref="SphericalHarmonic.Circle"/>.
        /// </summary>
        internal CircularEngine()
        {
            _M = -1;
            _gradp = true;
            _u = 0;
            _t = 1;
        }

        /// <summary>
        /// Evaluate the sum for a particular longitude given in terms of its sine and cosine.
        /// </summary>
        /// <param name="sinlon">the sine of the longitude.</param>
        /// <param name="coslon">the cosine of the longitude.</param>
        /// <returns><i>V</i>, the value of the sum.</returns>
        /// <remarks>
        /// The arguments must satisfy <paramref name="sinlon"/>^2 + <paramref name="coslon"/>^2 = 1.
        /// </remarks>
        public double Evaluate(double sinlon, double coslon) => Value(false, sinlon, coslon, out _, out _, out _);

        /// <summary>
        /// Evaluate the sum for a particular longitude.
        /// </summary>
        /// <param name="lon">the longitude (degrees).</param>
        /// <returns><i>V</i>, the value of the sum.</returns>
        public double Evaluate(double lon)
        {
            SinCosd(lon, out var sinlon, out var coslon);
            return Evaluate(sinlon, coslon);
        }

        /// <summary>
        /// Evaluate the sum and its gradient for a particular longitude given in terms of its sine and cosine.
        /// </summary>
        /// <param name="sinlon">the sine of the longitude.</param>
        /// <param name="coslon">the cosine of the longitude.</param>
        /// <param name="gradx"><i>x</i> component of the gradient.</param>
        /// <param name="grady"><i>y</i> component of the gradient.</param>
        /// <param name="gradz"><i>z</i> component of the gradient.</param>
        /// <returns><i>V</i>, the value of the sum.</returns>
        /// <remarks>
        /// The gradients will only be computed if the CircularEngine object was created with this capability
        /// (e.g., via <i>gradp</i> = <see langword="true"/> in <see cref="SphericalHarmonic.Circle"/>).
        /// If not, <paramref name="gradx"/>, etc., will be set to <see cref="double.NaN"/>.
        /// The arguments must satisfy <paramref name="sinlon"/>^2 + <paramref name="coslon"/>^2 = 1.
        /// </remarks>
        public double Evaluate(double sinlon, double coslon,
                 out double gradx, out double grady, out double gradz)
            => Value(true, sinlon, coslon, out gradx, out grady, out gradz);

        internal void SetCoeff(int m, double wc, double ws)
        {
            _wc.Span[m] = wc;
            _ws.Span[m] = ws;
        }

        internal void SetCoeff(int m, double wc, double ws,
                      double wrc, double wrs, double wtc, double wts)
        {
            _wc.Span[m] = wc;
            _ws.Span[m] = ws;
            if (_gradp)
            {
                _wrc.Span[m] = wrc; _wrs.Span[m] = wrs;
                _wtc.Span[m] = wtc; _wts.Span[m] = wts;
            }
        }

        private double Value(bool gradp, double sl, double cl,
                 out double gradx, out double grady, out double gradz)
        {
            gradx = grady = gradz = double.NaN;
            gradp = _gradp && gradp;
            var root = SphericalEngine.SqrtTable;

            // Initialize outer sum
            double vc = 0, vc2 = 0, vs = 0, vs2 = 0;       // v [N + 1], v [N + 2]
                                                           // vr, vt, vl and similar w variable accumulate the sums for the
                                                           // derivatives wrt r, theta, and lambda, respectively.
            double vrc = 0, vrc2 = 0, vrs = 0, vrs2 = 0;   // vr[N + 1], vr[N + 2]
            double vtc = 0, vtc2 = 0, vts = 0, vts2 = 0;   // vt[N + 1], vt[N + 2]
            double vlc = 0, vlc2 = 0, vls = 0, vls2 = 0;   // vl[N + 1], vl[N + 2]
            for (int m = _M; m >= 0; --m)
            {   // m = M .. 0
                // Now Sc[m] = wc, Ss[m] = ws
                // Sc'[m] = wtc, Ss'[m] = wtc
                if (m > 0)
                {
                    double v, A, B;           // alpha[m], beta[m + 1]
                    switch (_norm)
                    {
                        case Normalization.Full:
                            v = root[2] * root[2 * m + 3] / root[m + 1];
                            A = cl * v * _uq;
                            B = -v * root[2 * m + 5] / (root[8] * root[m + 2]) * _uq2;
                            break;
                        case Normalization.Schmidt:
                            v = root[2] * root[2 * m + 1] / root[m + 1];
                            A = cl * v * _uq;
                            B = -v * root[2 * m + 3] / (root[8] * root[m + 2]) * _uq2;
                            break;
                        default:
                            A = B = 0; break;
                    }
                    v = A * vc + B * vc2 + _wc.Span[m]; vc2 = vc; vc = v;
                    v = A * vs + B * vs2 + _ws.Span[m]; vs2 = vs; vs = v;
                    if (gradp)
                    {
                        v = A * vrc + B * vrc2 + _wrc.Span[m]; vrc2 = vrc; vrc = v;
                        v = A * vrs + B * vrs2 + _wrs.Span[m]; vrs2 = vrs; vrs = v;
                        v = A * vtc + B * vtc2 + _wtc.Span[m]; vtc2 = vtc; vtc = v;
                        v = A * vts + B * vts2 + _wts.Span[m]; vts2 = vts; vts = v;
                        v = A * vlc + B * vlc2 + m * _ws.Span[m]; vlc2 = vlc; vlc = v;
                        v = A * vls + B * vls2 - m * _wc.Span[m]; vls2 = vls; vls = v;
                    }
                }
                else
                {
                    double A, B, qs;
                    switch (_norm)
                    {
                        case Normalization.Full:
                            A = root[3] * _uq;       // F[1]/(q*cl) or F[1]/(q*sl)
                            B = -root[15] / 2 * _uq2; // beta[1]/q
                            break;
                        case Normalization.Schmidt:
                            A = _uq;
                            B = -root[3] / 2 * _uq2;
                            break;
                        default:
                            A = B = 0; break;
                    }
                    qs = _q / SphericalEngine.Scale;
                    vc = qs * (_wc.Span[m] + A * (cl * vc + sl * vs) + B * vc2);
                    if (gradp)
                    {
                        qs /= _r;
                        // The components of the gradient in circular coordinates are
                        // r: dV/dr
                        // theta: 1/r * dV/dtheta
                        // lambda: 1/(r*u) * dV/dlambda
                        vrc = -qs * (_wrc.Span[m] + A * (cl * vrc + sl * vrs) + B * vrc2);
                        vtc = qs * (_wtc.Span[m] + A * (cl * vtc + sl * vts) + B * vtc2);
                        vlc = qs / _u * (A * (cl * vlc + sl * vls) + B * vlc2);
                    }
                }
            }

            if (gradp)
            {
                // Rotate into cartesian (geocentric) coordinates
                gradx = cl * (_u * vrc + _t * vtc) - sl * vlc;
                grady = sl * (_u * vrc + _t * vtc) + cl * vlc;
                gradz = _t * vrc - _u * vtc;
            }
            return vc;
        }
    }
}
