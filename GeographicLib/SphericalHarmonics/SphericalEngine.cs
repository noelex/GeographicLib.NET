using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.IO;
using System.Text;

using static System.Math;
using static GeographicLib.MathEx;

namespace GeographicLib.SphericalHarmonics
{
    /// <summary>
    /// The evaluation engine for <see cref="SphericalHarmonic"/>.
    /// </summary>
    /// <remarks>
    /// This serves as the backend to <see cref="SphericalHarmonic"/>, <see cref="SphericalHarmonic1"/>, and <see cref="SphericalHarmonic2"/>.
    /// Typically end-users will not have to access this class directly.
    /// </remarks>
    public static class SphericalEngine
    {
        private static readonly List<double> sqrttable_ = new List<double>();

        /// <summary>
        /// Gets a value representing a table of the square roots of integers.
        /// </summary>
        internal static ReadOnlyCollection<double> SqrtTable { get; } = new ReadOnlyCollection<double>(sqrttable_);

        /// <summary>
        /// Gets a value representing an internal scaling of the coefficients to avoid overflow in intermediate calculations.
        /// </summary>
        internal static double Scale { get; }
            = Pow(FLT_RADIX,
                     -3 * (DBL_MAX_EXP < (1 << 14) ?
                           DBL_MAX_EXP : (1 << 14))
                     / 5);

        /// <summary>
        /// Move latitudes near the pole off the axis by this amount.
        /// </summary>
        internal static double Eps { get; } = DBL_EPSILON * Sqrt(DBL_EPSILON);

        /// <summary>
        /// Check that the static table of square roots is big enough and enlarge it if necessary.
        /// </summary>
        /// <param name="N">
        /// the maximum degree to be used in <see cref="SphericalEngine"/>.
        /// </param>
        /// <remarks>
        /// Typically, there's no need for an end-user to call this routine, because the constructors for <see cref="Coeff"/> do so.
        /// However, since this updates a static table, there's a possible race condition in a multi-threaded environment.
        /// Because this routine does nothing if the table is already large enough, one way to avoid race conditions is to call this
        /// routine at program start up (when it's still single threaded), supplying the largest degree that your program will use. E.g.,
        /// <code>SphericalEngine.RootTable(2190);</code>
        /// suffices to accommodate extant magnetic and gravity models.
        /// </remarks>
        public static void RootTable(int N)
        {
            // Need square roots up to max(2 * N + 5, 15).
            int L = Max(2 * N + 5, 15) + 1, oldL = sqrttable_.Count;
            if (oldL >= L)
                return;

            for (int l = oldL; l < L; ++l)
                sqrttable_.Add(Sqrt(l));
        }

        /// <summary>
        /// Clear the static table of square roots and release the memory. 
        /// Call this only when you are sure you no longer will be using <see cref="SphericalEngine"/>. 
        /// Your program will crash if you call <see cref="SphericalEngine"/> after calling this routine.
        /// </summary>
        /// <remarks>
        /// <para>WARNING: It's safest not to call this routine at all. (The space used by the table is modest.)</para>
        /// </remarks>
        public static void ClearRootTable() => sqrttable_.Clear();

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c">an array of coeff objects.</param>
        /// <param name="f">array of coefficient multipliers.  <c>f[0]</c> should be <c>1</c>.</param>
        /// <param name="x">the <i>x</i> component of the cartesian position.</param>
        /// <param name="y">the <i>y</i> component of the cartesian position.</param>
        /// <param name="z">the <i>z</i> component of the cartesian position.</param>
        /// <param name="a">the normalizing radius.</param>
        /// <param name="gradx">the <i>x</i> component of the gradient.</param>
        /// <param name="grady">the <i>y</i> component of the gradient.</param>
        /// <param name="gradz">the <i>z</i> component of the gradient.</param>
        /// <param name="gradp">should the gradient be calculated.</param>
        /// <param name="norm">the normalization for the associated Legendre polynomials.</param>
        /// <returns>the spherical harmonic sum.</returns>
        /// <remarks>
        /// See the SphericalHarmonic class for the definition of the sum. 
        /// The coefficients used by this function are, for example, <c>c[0].Cv + f[1] * c[1].Cv + ... + f[L−1] * c[L−1].Cv</c>.
        /// (Note that <c>f[0]</c> is <i>not</i> used.) The upper limits on the sum are determined by <c>c[0].Nmx</c> and <c>c[0].Mmx</c>;
        /// these limits apply to <i>all</i> the components of the coefficients.
        /// <para>
        /// Clenshaw summation is used which permits the evaluation of the sum without the need to allocate temporary arrays. 
        /// Thus this function never throws an exception.
        /// </para>
        /// </remarks>
        public static double Value(bool gradp, Normalization norm,
                                   ReadOnlySpan<Coeff> c, ReadOnlySpan<double> f,
                                   double x, double y, double z, double a,
                                   out double gradx, out double grady, out double gradz)
        {
            var L = c.Length;
            Debug.Assert(Enum.IsDefined(norm.GetType(), norm), "Unknown normalization");

            gradx = grady = gradz = double.NaN;

            int N = c[0].Nmx, M = c[0].Mmx;

            double
              p = Hypot(x, y),
              cl = p != 0 ? x / p : 1,  // cos(lambda); at pole, pick lambda = 0
              sl = p != 0 ? y / p : 0,  // sin(lambda)
              r = Hypot(z, p),
              t = r != 0 ? z / r : 0,   // cos(theta); at origin, pick theta = pi/2
              u = r != 0 ? Max(p / r, Eps) : 1, // sin(theta); but avoid the pole
              q = a / r;
            double
              q2 = Sq(q),
              uq = u * q,
              uq2 = Sq(uq),
              tu = t / u;
            // Initialize outer sum
            double vc = 0, vc2 = 0, vs = 0, vs2 = 0;       // v [N + 1], v [N + 2]
                                                           // vr, vt, vl and similar w variable accumulate the sums for the
                                                           // derivatives wrt r, theta, and lambda, respectively.
            double vrc = 0, vrc2 = 0, vrs = 0, vrs2 = 0;   // vr[N + 1], vr[N + 2]
            double vtc = 0, vtc2 = 0, vts = 0, vts2 = 0;   // vt[N + 1], vt[N + 2]
            double vlc = 0, vlc2 = 0, vls = 0, vls2 = 0;   // vl[N + 1], vl[N + 2]

            Span<int> k = stackalloc int[L];
            var root = SqrtTable;

            for (int m = M; m >= 0; --m)
            {   // m = M .. 0
                // Initialize inner sum
                double
                  wc = 0, wc2 = 0, ws = 0, ws2 = 0, // w [N - m + 1], w [N - m + 2]
                  wrc = 0, wrc2 = 0, wrs = 0, wrs2 = 0, // wr[N - m + 1], wr[N - m + 2]
                  wtc = 0, wtc2 = 0, wts = 0, wts2 = 0; // wt[N - m + 1], wt[N - m + 2]
                for (int l = 0; l < L; ++l)
                    k[l] = c[l].IndexOf(N, m) + 1;
                for (int n = N; n >= m; --n)
                {             // n = N .. m; l = N - m .. 0
                    double w, A = 0, Ax = 0, B = 0, R;    // alpha[l], beta[l + 1]
                    switch (norm)
                    {
                        case Normalization.Full:
                            w = root[2 * n + 1] / (root[n - m + 1] * root[n + m + 1]);
                            Ax = q * w * root[2 * n + 3];
                            A = t * Ax;
                            B = -q2 * root[2 * n + 5] /
                              (w * root[n - m + 2] * root[n + m + 2]);
                            break;
                        case Normalization.Schmidt:
                            w = root[n - m + 1] * root[n + m + 1];
                            Ax = q * (2 * n + 1) / w;
                            A = t * Ax;
                            B = -q2 * w / (root[n - m + 2] * root[n + m + 2]);
                            break;
                        default: break;       // To suppress warning message from Visual Studio
                    }
                    R = c[0].Cv(--k[0]);
                    for (int l = 1; l < L; ++l)
                        R += c[l].Cv(--k[l], n, m, f[l]);
                    R *= Scale;
                    w = A * wc + B * wc2 + R; wc2 = wc; wc = w;
                    if (gradp)
                    {
                        w = A * wrc + B * wrc2 + (n + 1) * R; wrc2 = wrc; wrc = w;
                        w = A * wtc + B * wtc2 - u * Ax * wc2; wtc2 = wtc; wtc = w;
                    }
                    if (m != 0)
                    {
                        R = c[0].Sv(k[0]);
                        for (int l = 1; l < L; ++l)
                            R += c[l].Sv(k[l], n, m, f[l]);
                        R *= Scale;
                        w = A * ws + B * ws2 + R; ws2 = ws; ws = w;
                        if (gradp)
                        {
                            w = A * wrs + B * wrs2 + (n + 1) * R; wrs2 = wrs; wrs = w;
                            w = A * wts + B * wts2 - u * Ax * ws2; wts2 = wts; wts = w;
                        }
                    }
                }
                // Now Sc[m] = wc, Ss[m] = ws
                // Sc'[m] = wtc, Ss'[m] = wtc
                if (m != 0)
                {
                    double v = 0, A = 0, B = 0;           // alpha[m], beta[m + 1]
                    switch (norm)
                    {
                        case Normalization.Full:
                            v = root[2] * root[2 * m + 3] / root[m + 1];
                            A = cl * v * uq;
                            B = -v * root[2 * m + 5] / (root[8] * root[m + 2]) * uq2;
                            break;
                        case Normalization.Schmidt:
                            v = root[2] * root[2 * m + 1] / root[m + 1];
                            A = cl * v * uq;
                            B = -v * root[2 * m + 3] / (root[8] * root[m + 2]) * uq2;
                            break;
                        default: break;       // To suppress warning message from Visual Studio
                    }
                    v = A * vc + B * vc2 + wc; vc2 = vc; vc = v;
                    v = A * vs + B * vs2 + ws; vs2 = vs; vs = v;
                    if (gradp)
                    {
                        // Include the terms Sc[m] * P'[m,m](t) and Ss[m] * P'[m,m](t)
                        wtc += m * tu * wc; wts += m * tu * ws;
                        v = A * vrc + B * vrc2 + wrc; vrc2 = vrc; vrc = v;
                        v = A * vrs + B * vrs2 + wrs; vrs2 = vrs; vrs = v;
                        v = A * vtc + B * vtc2 + wtc; vtc2 = vtc; vtc = v;
                        v = A * vts + B * vts2 + wts; vts2 = vts; vts = v;
                        v = A * vlc + B * vlc2 + m * ws; vlc2 = vlc; vlc = v;
                        v = A * vls + B * vls2 - m * wc; vls2 = vls; vls = v;
                    }
                }
                else
                {
                    double A = 0, B = 0, qs;
                    switch (norm)
                    {
                        case Normalization.Full:
                            A = root[3] * uq;       // F[1]/(q*cl) or F[1]/(q*sl)
                            B = -root[15] / 2 * uq2; // beta[1]/q
                            break;
                        case Normalization.Schmidt:
                            A = uq;
                            B = -root[3] / 2 * uq2;
                            break;
                        default: break;       // To suppress warning message from Visual Studio
                    }
                    qs = q / Scale;
                    vc = qs * (wc + A * (cl * vc + sl * vs) + B * vc2);
                    if (gradp)
                    {
                        qs /= r;
                        // The components of the gradient in spherical coordinates are
                        // r: dV/dr
                        // theta: 1/r * dV/dtheta
                        // lambda: 1/(r*u) * dV/dlambda
                        vrc = -qs * (wrc + A * (cl * vrc + sl * vrs) + B * vrc2);
                        vtc = qs * (wtc + A * (cl * vtc + sl * vts) + B * vtc2);
                        vlc = qs / u * (A * (cl * vlc + sl * vls) + B * vlc2);
                    }
                }
            }

            if (gradp)
            {
                // Rotate into cartesian (geocentric) coordinates
                gradx = cl * (u * vrc + t * vtc) - sl * vlc;
                grady = sl * (u * vrc + t * vtc) + cl * vlc;
                gradz = t * vrc - u * vtc;
            }
            return vc;
        }

        /// <summary>
        /// Create a <see cref="CircularEngine"/> object.
        /// </summary>
        /// <param name="gradp">should the gradient be calculated.</param>
        /// <param name="norm">the normalization for the associated Legendre polynomials.</param>
        /// <param name="c">an array of coeff objects.</param>
        /// <param name="f">array of coefficient multipliers.  <c>f[0]</c> should be <c>1</c>.</param>
        /// <param name="p">the radius of the circle = sqrt(<i>x</i>^2 + <i>y</i>^2).</param>
        /// <param name="z">the height of the circle.</param>
        /// <param name="a">the normalizing radius.</param>
        /// <returns>A <see cref="CircularEngine"/> instance.</returns>
        /// <remarks>
        /// If you need to evaluate the spherical harmonic sum for several points with constant 
        /// <i>f</i>, <i>p</i> = sqrt(<i>x</i>^2 + <i>y</i>^2), <i>z</i>, and <i>a</i>, it is more efficient to construct call
        /// <see cref="Circle(bool, Normalization, ReadOnlySpan{Coeff}, ReadOnlySpan{double}, double, double, double)"/>
        /// to give a <see cref="CircularEngine"/> object and then call <see cref="CircularEngine.Evaluate(double,double)"/> 
        /// with arguments <i>x</i>/<i>p</i> and <i>y</i>/<i>p</i>.
        /// </remarks>
        public static CircularEngine Circle(bool gradp, Normalization norm,
                                   ReadOnlySpan<Coeff> c, ReadOnlySpan<double> f,
                                   double p, double z, double a)
        {
            var L = c.Length;
            Debug.Assert(Enum.IsDefined(norm.GetType(), norm), "Unknown normalization");

            int N = c[0].Nmx, M = c[0].Mmx;

            double
              r = Hypot(z, p),
              t = r != 0 ? z / r : 0,   // cos(theta); at origin, pick theta = pi/2
              u = r != 0 ? Max(p / r, Eps) : 1, // sin(theta); but avoid the pole
              q = a / r;
            double
              q2 = Sq(q),
              tu = t / u;
            var circ = new CircularEngine(M, gradp, norm, a, r, u, t);
            Span<int> k = stackalloc int[L];
            var root = SqrtTable;

            for (int m = M; m >= 0; --m)
            {   // m = M .. 0
                // Initialize inner sum
                double
                  wc = 0, wc2 = 0, ws = 0, ws2 = 0, // w [N - m + 1], w [N - m + 2]
                  wrc = 0, wrc2 = 0, wrs = 0, wrs2 = 0, // wr[N - m + 1], wr[N - m + 2]
                  wtc = 0, wtc2 = 0, wts = 0, wts2 = 0; // wt[N - m + 1], wt[N - m + 2]
                for (int l = 0; l < L; ++l)
                    k[l] = c[l].IndexOf(N, m) + 1;
                for (int n = N; n >= m; --n)
                {             // n = N .. m; l = N - m .. 0
                    double w, A = 0, Ax = 0, B = 0, R;    // alpha[l], beta[l + 1]
                    switch (norm)
                    {
                        case Normalization.Full:
                            w = root[2 * n + 1] / (root[n - m + 1] * root[n + m + 1]);
                            Ax = q * w * root[2 * n + 3];
                            A = t * Ax;
                            B = -q2 * root[2 * n + 5] /
                              (w * root[n - m + 2] * root[n + m + 2]);
                            break;
                        case Normalization.Schmidt:
                            w = root[n - m + 1] * root[n + m + 1];
                            Ax = q * (2 * n + 1) / w;
                            A = t * Ax;
                            B = -q2 * w / (root[n - m + 2] * root[n + m + 2]);
                            break;
                        default: break;       // To suppress warning message from Visual Studio
                    }
                    R = c[0].Cv(--k[0]);
                    for (int l = 1; l < L; ++l)
                        R += c[l].Cv(--k[l], n, m, f[l]);
                    R *= Scale;
                    w = A * wc + B * wc2 + R; wc2 = wc; wc = w;
                    if (gradp)
                    {
                        w = A * wrc + B * wrc2 + (n + 1) * R; wrc2 = wrc; wrc = w;
                        w = A * wtc + B * wtc2 - u * Ax * wc2; wtc2 = wtc; wtc = w;
                    }
                    if (m != 0)
                    {
                        R = c[0].Sv(k[0]);
                        for (int l = 1; l < L; ++l)
                            R += c[l].Sv(k[l], n, m, f[l]);
                        R *= Scale;
                        w = A * ws + B * ws2 + R; ws2 = ws; ws = w;
                        if (gradp)
                        {
                            w = A * wrs + B * wrs2 + (n + 1) * R; wrs2 = wrs; wrs = w;
                            w = A * wts + B * wts2 - u * Ax * ws2; wts2 = wts; wts = w;
                        }
                    }
                }
                if (!gradp)
                    circ.SetCoeff(m, wc, ws);
                else
                {
                    // Include the terms Sc[m] * P'[m,m](t) and  Ss[m] * P'[m,m](t)
                    wtc += m * tu * wc; wts += m * tu * ws;
                    circ.SetCoeff(m, wc, ws, wrc, wrs, wtc, wts);
                }
            }

            return circ;
        }

        /// <summary>
        /// Package up coefficients for <see cref="SphericalEngine"/>.
        /// </summary>
        public readonly struct Coeff
        {
            private readonly int _Nx, _nmx, _mmx;
            private readonly ReadOnlyMemory<double> _Cnm, _Snm;

            /// <summary>
            /// The constructor for full coefficient vectors.
            /// </summary>
            /// <param name="C">a vector of coefficients for the cosine terms.</param>
            /// <param name="S">a vector of coefficients for the sine terms.</param>
            /// <param name="N">the degree giving storage layout for <paramref name="C"/> and <paramref name="S"/>.</param>
            public Coeff(ReadOnlyMemory<double> C, ReadOnlyMemory<double> S, int N)
            {
                _Nx = N;
                _nmx = N;
                _mmx = N;
                _Cnm = C;
                _Snm = S;

                if (!(_Nx >= -1))
                    throw new GeographicException("Bad indices for coeff");
                if (!(IndexOf(_nmx, _mmx) < C.Length &&
                      IndexOf(_nmx, _mmx) < S.Length + (_Nx + 1)))
                    throw new GeographicException("Arrays too small in coeff");

                RootTable(_nmx);
            }

            /// <summary>
            /// The general constructor.
            /// </summary>
            /// <param name="C">a vector of coefficients for the cosine terms.</param>
            /// <param name="S">a vector of coefficients for the sine terms.</param>
            /// <param name="N">the degree giving storage layout for <paramref name="C"/> and <paramref name="S"/>.</param>
            /// <param name="nmx">the maximum degree to be used.</param>
            /// <param name="mmx">the maximum order to be used.</param>
            public Coeff(ReadOnlyMemory<double> C, ReadOnlyMemory<double> S, int N, int nmx, int mmx)
            {
                _Nx = N;
                _nmx = nmx;
                _mmx = mmx;
                _Cnm = C;
                _Snm = S;

                if (!((_Nx >= _nmx && _nmx >= _mmx && _mmx >= 0) ||
                  // If mmx = -1 then the sums are empty so require nmx = -1 also.
                  (_nmx == -1 && _mmx == -1)))
                    throw new GeographicException("Bad indices for coeff");
                if (!(IndexOf(_nmx, _mmx) < C.Length &&
                      IndexOf(_nmx, _mmx) < S.Length + (_Nx + 1)))
                    throw new GeographicException("Arrays too small in coeff");

                RootTable(_nmx);
            }

            /// <summary>
            /// Get one-dimensional index into <i>C</i> and <i>S</i>.
            /// </summary>
            /// <param name="n">the degree.</param>
            /// <param name="m">the order.</param>
            /// <returns>the one-dimensional index.</returns>
            public int IndexOf(int n, int m) => m * _Nx - m * (m - 1) / 2 + n;

            /// <summary>
            /// Gets a value representing the degree giving storage layout for <i>C</i> and <i>S</i>.
            /// </summary>
            public int N => _Nx;

            /// <summary>
            /// Gets a value representing the maximum degree to be used.
            /// </summary>
            public int Nmx => _nmx;

            /// <summary>
            /// Gets a value representing the maximum order to be used.
            /// </summary>
            public int Mmx => _mmx;

            /// <summary>
            /// Gets an element of <i>C</i>.
            /// </summary>
            /// <param name="k">the one-dimensional index.</param>
            /// <returns>the value of the <i>C</i> coefficient.</returns>
            public double Cv(int k) => _Cnm.Span[k];

            /// <summary>
            /// Gets an element of <i>S</i>.
            /// </summary>
            /// <param name="k">the one-dimensional index.</param>
            /// <returns>the value of the <i>S</i> coefficient.</returns>
            public double Sv(int k) => _Snm.Span[k - (_Nx + 1)];

            /// <summary>
            /// Gets an element of <i>C</i> with checking.
            /// </summary>
            /// <param name="k">the one-dimensional index.</param>
            /// <param name="n">the requested degree.</param>
            /// <param name="m">the requested order.</param>
            /// <param name="f">a multiplier.</param>
            /// <returns>
            /// the value of the <i>C</i> coefficient multiplied by <paramref name="f"/> in <paramref name="n"/> and <paramref name="m"/>
            /// are in range else <c>0</c>.
            /// </returns>
            public double Cv(int k, int n, int m, double f) => m > _mmx || n > _nmx ? 0 : _Cnm.Span[k] * f;

            /// <summary>
            /// Gets an element of <i>S</i> with checking.
            /// </summary>
            /// <param name="k">the one-dimensional index.</param>
            /// <param name="n">the requested degree.</param>
            /// <param name="m">the requested order.</param>
            /// <param name="f">a multiplier.</param>
            /// <returns>
            /// the value of the <i>S</i> coefficient multiplied by <paramref name="f"/> in <paramref name="n"/> and <paramref name="m"/>
            /// are in range else <c>0</c>.
            /// </returns>
            public double Sv(int k, int n, int m, double f) => m > _mmx || n > _nmx ? 0 : _Snm.Span[k - (_Nx + 1)] * f;

            /// <summary>
            /// Gets the size of the coefficient vector for the cosine terms.
            /// </summary>
            /// <param name="N">the maximum degree.</param>
            /// <param name="M">the maximum order.</param>
            /// <returns>the size of the vector of cosine terms as stored in column major order.</returns>
            public static int Csize(int N, int M) => (M + 1) * (2 * N - M + 2) / 2;

            /// <summary>
            /// Gets size of the coefficient vector for the sine terms.
            /// </summary>
            /// <param name="N">the maximum degree.</param>
            /// <param name="M">the maximum order.</param>
            /// <returns>the size of the vector of sine terms as stored in column major order.</returns>
            public static int Ssize(int N, int M) => Csize(N, M) - (N + 1);

            /// <summary>
            /// Load coefficients from a binary stream.
            /// </summary>
            /// <param name="stream">the input stream.</param>
            /// <param name="N">The maximum degree of the coefficients.</param>
            /// <param name="M">The maximum order of the coefficients.</param>
            /// <param name="truncate">
            /// if <see langword="false"/> (the default) then <paramref name="N"/> and <paramref name="M"/> are determined by the values in the binary stream; 
            /// otherwise, the input values of <paramref name="N"/> and <paramref name="M"/> are used to truncate the coefficients read from the stream at the given degree and order.
            /// </param>
            /// <remarks>
            /// <paramref name="N"/> and <paramref name="M"/> are read as 4-byte ints.
            /// are resized to accommodate all the coefficients (with the <i>m</i> = 0 coefficients for <i>S</i> excluded)
            /// and the data for these coefficients read as 8-byte doubles. The coefficients are stored in column major order.
            /// The bytes in the stream should use little-endian ordering. IEEE floating point is assumed for the coefficients.
            /// </remarks>
            /// <returns>
            /// <i>C</i>, the vector of cosine coefficients and <i>S</i>, the vector of sine coefficients.
            /// </returns>
            public static Coeff FromStream(Stream stream,
                ref int N, ref int M, bool truncate = false)
            {
                if (truncate)
                {
                    if (!((N >= M && M >= 0) || (N == -1 && M == -1)))
                        // The last condition is that M = -1 implies N = -1.
                        throw new GeographicException($"Bad requested degree and order {N} {M}");
                }

                Span<int> nm = stackalloc int[2];
                Utility.ReadArray(stream, nm);

                int N0 = nm[0], M0 = nm[1];
                if (!((N0 >= M0 && M0 >= 0) || (N0 == -1 && M0 == -1)))
                    // The last condition is that M0 = -1 implies N0 = -1.
                    throw new GeographicException($"Bad degree and order {N0} {M0}");
                N = truncate ? Min(N, N0) : N0;
                M = truncate ? Min(M, M0) : M0;

                Memory<double> C = new double[Csize(N, M)],
                               S = new double[Ssize(N, M)];

                int skip = (Csize(N0, M0) -
                            Csize(N0, M)) * sizeof(double);
                if (N == N0)
                {
                    Utility.ReadArray(stream, C.Span);
                    if (skip != 0) stream.Seek(skip, SeekOrigin.Current);
                    Utility.ReadArray(stream, S.Span);
                    if (skip != 0) stream.Seek(skip, SeekOrigin.Current);
                }
                else
                {
                    for (int m = 0, k = 0; m <= M; ++m)
                    {
                        Utility.ReadArray(stream, C.Span.Slice(k));
                        stream.Seek((N0 - N) * sizeof(double), SeekOrigin.Current);
                        k += N + 1 - m;
                    }
                    if (skip != 0) stream.Seek(skip, SeekOrigin.Current);
                    for (int m = 1, k = 0; m <= M; ++m)
                    {
                        Utility.ReadArray(stream, S.Span.Slice(k));
                        stream.Seek((N0 - N) * sizeof(double), SeekOrigin.Current);
                        k += N + 1 - m;
                    }
                    if (skip != 0) stream.Seek(skip, SeekOrigin.Current);
                }

                return new Coeff(C, S, N, N, M);
            }
        }
    }
}
