using System;
using System.Collections.Generic;
using System.Numerics;
using System.Text;

using static System.Math;
using static GeographicLib.Macros;
using static GeographicLib.MathEx;

namespace GeographicLib
{
    /// <summary>
    /// Exact geodesic calculations.
    /// </summary>
    /// <remarks>
    /// <para>
    /// The equations for geodesics on an ellipsoid can be expressed in terms of
    /// incomplete elliptic integrals.The Geodesic class expands these integrals
    /// in a series in the flattening <i>f</i> and this provides an accurate solution
    /// for <i>f</i> ∈ [-0.01, 0.01].  The GeodesicExact class computes the
    /// ellitpic integrals directly and so provides a solution which is valid for
    /// all <i>f</i>.However, in practice, its use should be limited to about
    /// <i>b</i>/<i>a</i> ∈ [0.01, 100] or <i>f</i> ∈ [-99, 0.99].
    /// </para>
    /// For the WGS84 ellipsoid, these classes are 2-3 times slower than the
    /// series solution and 2-3 times less accurate(because it's less easy
    /// to control round-off errors with the elliptic integral formulation); i.e.,
    /// the error is about 40 nm(40 nanometers) instead of 15 nm.However the
    /// error in the series solution scales as <i>f</i>^7 while the
    /// error in the elliptic integral solution depends weakly on <i>f</i>.If the
    /// quarter meridian distance is 10000 km and the ratio <i>b</i>/<i>a</i> = 1
    /// - <i>f</i> is varied then the approximate maximum error(expressed as a
    /// distance) is
    /// <code>
    ///       1 - f error(nm)
    ///       1/128     387
    ///       1/64      345
    ///       1/32      269
    ///       1/16      210
    ///       1/8       115
    ///       1/4        69
    ///       1/2        36
    ///         1        15
    ///         2        25
    ///         4        96
    ///         8       318
    ///        16       985
    ///        32      2352
    ///        64      6008
    ///       128     19024
    /// </code>
    /// The area in this classes is computing by finding an accurate approximation
    /// to the area integrand using a discrete sine transform fitting <i>N</i> equally
    /// spaced points in σ. <i>N</i> chosen to ensure full accuracy for
    /// <i>b</i>/<i>a</i> ∈ [0.01, 100] or <i>f</i> ∈ [-99, 0.99].
    /// <para>
    /// See <a href="https://geographiclib.sourceforge.io/C++/doc/geodesic.html#geodellip">
    /// Geodesics in terms of elliptic integrals</a>
    /// for the formulation.See the documentation on the
    /// <see cref="Geodesic"/> class for additional information on the geodesic problems.
    /// </para>
    /// </remarks>
    public partial class GeodesicExact : GeodesicBase
    {
        private const int ndiv = 100;
        private static readonly ReadOnlyMemory<byte> narr = new byte[]{
          22,21,19,18,17,17,16,15,15,15,14,14,14,13,13,13,13,13,13,12,12,12,12,12,
          12,11,11,11,11,11,11,11,11,11,11,10,10,10,10,10,10,10,10,9,9,9,9,9,9,9,9,
          9,9,9,9,9,9,8,8,8,8,8,8,8,8,8,8,7,7,7,7,7,7,7,7,7,7,7,7,7,7,6,6,6,6,6,6,
          6,6,5,5,5,5,5,5,5,5,4,4,3,2,3,4,4,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,7,7,
          7,7,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,
          9,9,9,9,9,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,12,12,
          12,12,12,13,13,13,13,13,14,14,15,15,16,17,18,19
        }; /* 2*ndiv+1 elements */

        private const uint
            maxit1_ = 20,
            maxit2_ = maxit1_ + DBL_MANT_DIG + 10;

        private const double tol0_ = DBL_EPSILON;
        internal static readonly double
            tiny_ = Sqrt(DBL_MIN),
            tol1_ = 200 * tol0_,
            tol2_ = Sqrt(tol0_),
            tolb_ = tol0_, // Check on bisection interval
            xthresh_ = 1000 * tol2_;

        internal readonly double _a, _f, _f1, _e2, _ep2, _n, _b, _c2, _etol2;
        internal readonly int _nC4;
        internal readonly DST _fft;

        /// <summary>
        /// Initialize a new <see cref="GeodesicExact"/> instance with specified ellipoid.
        /// </summary>
        /// <param name="ellipsoid">Source <see cref="IEllipsoid"/> object.</param>
        public GeodesicExact(IEllipsoid ellipsoid) : this(ellipsoid.EquatorialRadius, ellipsoid.Flattening) { }

        /// <summary>
        /// Initialize a new <see cref="GeodesicExact"/> instance with specified equatorial radius and flattening of the ellipsoid.
        /// </summary>
        /// <param name="a">equatorial radius (meters).</param>
        /// <param name="f">flattening of ellipsoid.  Setting <i>f</i> = 0 gives a sphere.</param>
        public GeodesicExact(double a, double f)
        {
            _a = a;
            _f = f;
            _f1 = 1 - f;
            _e2 = _f * (2 - _f);
            _ep2 = _e2 / Sq(_f1); // e2 / (1 - e2)
            _n = _f / (2 - _f);
            _b = _a * _f1;

            // The Geodesic class substitutes atanh(sqrt(e2)) for asinh(sqrt(ep2)) in
            // the definition of _c2.  The latter is more accurate for very oblate
            // ellipsoids (which the Geodesic class does not attempt to handle).
            _c2 = (Sq(_a) + Sq(_b) *
                  (_f == 0 ? 1 :
                   (_f > 0 ? Asinh(Sqrt(_ep2)) : Atan(Sqrt(-_e2))) /
                   Sqrt(Abs(_e2)))) / 2; // authalic radius squared

            // The sig12 threshold for "really short".  Using the auxiliary sphere
            // solution with dnm computed at (bet1 + bet2) / 2, the relative error in
            // the azimuth consistency check is sig12^2 * abs(f) * min(1, 1-f/2) / 2.
            // (Error measured for 1/100 < b/a < 100 and abs(f) >= 1/1000.  For a
            // given f and sig12, the max error occurs for lines near the pole.  If
            // the old rule for computing dnm = (dn1 + dn2)/2 is used, then the error
            // increases by a factor of 2.)  Setting this equal to epsilon gives
            // sig12 = etol2.  Here 0.1 is a safety factor (error decreased by 100)
            // and max(0.001, abs(f)) stops etol2 getting too large in the nearly
            // spherical case.
            _etol2 = 0.1 * tol2_ /
                    Sqrt(Max(0.001, Abs(_f)) * Min(1d, 1 - _f / 2) / 2);

            if (!(IsFinite(_a) && _a > 0))
                throw new GeographicException("Equatorial radius is not positive");

            if (!(IsFinite(_b) && _b > 0))
                throw new GeographicException("Polar semi-axis is not positive");

            // Required number of terms in DST for full accuracy for all precisions as
            // a function of n in [-0.99, 0.99].  Values determined by running
            // develop/AreaEst compiled with GEOGRAPHICLIB_PRECISION = 5.  For
            // precision 4 and 5, GEOGRAPHICLIB_DIGITS was set to, resp., 384 and 768.
            // The error criterion is relative error less than or equal to epsilon/2 =
            // 0.5^digits, with digits = 24, 53, 64, 113, 256.  The first 4 are the the
            // "standard" values for float, double, long double, and float128; the last
            // is the default for GeographicLib + mpfr.  Also listed is the value of
            // alp0 resulting in the most error for the given N.
            //
            //          float       double   long double    quad         mpfr
            //  n       N   alp0    N   alp0    N   alp0    N   alp0     N   alp0
            // -0.99  1024  0.09  3072  0.05  4096  0.04  8192 43.50  16384 22.01
            // -0.98   512  0.18  1536  0.10  2048  0.09  4096  0.06   8192  0.04
            // -0.97   384  0.25  1024  0.16  1536  0.13  3072  0.09   6144  0.06
            // -0.96   256  0.36   768  0.21  1024  0.18  2048  0.13   4096  0.09
            // -0.95   192  0.47   768  0.23   768  0.23  1536  0.17   4096  0.10
            // -0.94   192  0.51   512  0.31   768  0.26  1536  0.18   3072  0.13
            // -0.93   192  0.55   384  0.39   512  0.34  1024  0.24   3072  0.14
            // -0.92   128  0.73   384  0.42   512  0.37  1024  0.26   2048  0.18
            // -0.91   128  0.77   384  0.45   384  0.45   768  0.32   2048  0.19
            // -0.90    96  0.94   256  0.58   384  0.47   768  0.34   2048  0.21
            // -0.89    96  0.99   256  0.61   384  0.50   768  0.35   1536  0.25
            // -0.88    96  1.04   256  0.64   384  0.52   768  0.37   1536  0.26
            // -0.87    96  1.09   192  0.77   256  0.67   512  0.47   1536  0.27
            // -0.86    64  1.38   192  0.80   256  0.69   512  0.49   1536  0.28
            // -0.85    64  1.43   192  0.83   256  0.72   512  0.51   1024  0.36
            // -0.84    64  1.49   192  0.86   256  0.75   384  0.61   1024  0.37
            // -0.83    64  1.54   192  0.89   192  0.89   384  0.63   1024  0.39
            // -0.82    48  1.82   192  0.92   192  0.92   384  0.65   1024  0.40
            // -0.81    48  1.88   128  1.16   192  0.95   384  0.67   1024  0.41
            // -0.80    48  1.94   128  1.19   192  0.97   384  0.69    768  0.49
            // -0.79    48  1.99   128  1.23   192  1.00   384  0.71    768  0.50
            // -0.78    48  2.04   128  1.26   192  1.03   384  0.73    768  0.51
            // -0.77    48  2.10   128  1.29   192  1.05   256  0.91    768  0.53
            // -0.76    48  2.15   128  1.32   128  1.32   256  0.93    768  0.54
            // -0.75    48  2.20    96  1.56   128  1.35   256  0.96    768  0.55
            // -0.74    32  2.74    96  1.60   128  1.38   256  0.98    768  0.57
            // -0.73    32  2.81    96  1.63   128  1.41   256  1.00    768  0.58
            // -0.72    32  2.87    96  1.67   128  1.44   256  1.02    512  0.72
            // -0.71    32  2.93    96  1.70   128  1.47   192  1.20    512  0.74
            // -0.70    32  2.99    96  1.73    96  1.73   192  1.23    512  0.75
            // -0.69    32  3.05    96  1.77    96  1.77   192  1.25    512  0.77
            // -0.68    32  3.11    96  1.80    96  1.80   192  1.28    512  0.78
            // -0.67    24  3.64    96  1.84    96  1.84   192  1.30    512  0.80
            // -0.66    24  3.71    96  1.87    96  1.87   192  1.32    512  0.81
            // -0.65    24  3.77    64  2.33    96  1.90   192  1.35    384  0.95
            // -0.64    24  3.84    64  2.37    96  1.93   192  1.37    384  0.97
            // -0.63    24  3.90    64  2.41    96  1.97   192  1.39    384  0.98
            // -0.62    24  3.97    64  2.45    96  2.00   192  1.42    384  1.00
            // -0.61    24  4.04    64  2.49    96  2.03   192  1.44    384  1.02
            // -0.60    24  4.10    64  2.53    96  2.06   192  1.46    384  1.03
            // -0.59    24  4.16    64  2.57    64  2.57   128  1.82    384  1.05
            // -0.58    24  4.23    64  2.60    64  2.60   128  1.84    384  1.07
            // -0.57    24  4.29    48  3.05    64  2.64   128  1.87    384  1.08
            // -0.56    24  4.36    48  3.10    64  2.68   128  1.90    384  1.10
            // -0.55    16  5.38    48  3.14    64  2.72   128  1.93    384  1.11
            // -0.54    16  5.46    48  3.19    64  2.76   128  1.96    384  1.13
            // -0.53    16  5.54    48  3.23    64  2.80   128  1.98    256  1.40
            // -0.52    16  5.61    48  3.27    64  2.84   128  2.01    256  1.42
            // -0.51    16  5.69    48  3.32    64  2.88   128  2.04    256  1.44
            // -0.50    16  5.77    48  3.36    64  2.92    96  2.38    256  1.46
            // -0.49    16  5.85    48  3.41    48  3.41    96  2.42    256  1.48
            // -0.48    16  5.92    48  3.45    48  3.45    96  2.45    256  1.50
            // -0.47    16  6.00    48  3.50    48  3.50    96  2.48    256  1.52
            // -0.46    16  6.08    48  3.54    48  3.54    96  2.51    256  1.54
            // -0.45    12  7.06    48  3.59    48  3.59    96  2.54    256  1.56
            // -0.44    12  7.15    48  3.63    48  3.63    96  2.57    256  1.58
            // -0.43    12  7.24    32  4.49    48  3.68    96  2.61    256  1.60
            // -0.42    12  7.33    32  4.55    48  3.72    96  2.64    192  1.87
            // -0.41    12  7.42    32  4.60    48  3.77    96  2.67    192  1.89
            // -0.40    12  7.51    32  4.66    48  3.81    96  2.70    192  1.91
            // -0.39    12  7.60    32  4.71    48  3.86    96  2.73    192  1.94
            // -0.38    12  7.69    32  4.77    48  3.90    96  2.77    192  1.96
            // -0.37    12  7.77    32  4.82    48  3.95    96  2.80    192  1.98
            // -0.36    12  7.86    32  4.88    48  3.99    96  2.83    192  2.00
            // -0.35    12  7.95    32  4.94    32  4.94    64  3.50    192  2.03
            // -0.34    12  8.04    32  4.99    32  4.99    64  3.54    192  2.05
            // -0.33    12  8.13    24  5.81    32  5.05    64  3.58    192  2.07
            // -0.32    12  8.22    24  5.88    32  5.10    64  3.62    192  2.10
            // -0.31    12  8.31    24  5.94    32  5.16    64  3.66    192  2.12
            // -0.30     8 10.16    24  6.01    32  5.22    64  3.70    192  2.14
            // -0.29     8 10.27    24  6.07    32  5.27    64  3.74    192  2.17
            // -0.28     8 10.38    24  6.14    32  5.33    64  3.78    128  2.68
            // -0.27     8 10.49    24  6.20    32  5.39    64  3.82    128  2.71
            // -0.26     8 10.60    24  6.27    32  5.45    64  3.87    128  2.74
            // -0.25     8 10.72    24  6.34    32  5.50    48  4.51    128  2.77
            // -0.24     8 10.83    24  6.40    24  6.40    48  4.55    128  2.80
            // -0.23     8 10.94    24  6.47    24  6.47    48  4.60    128  2.83
            // -0.22     8 11.05    24  6.54    24  6.54    48  4.65    128  2.86
            // -0.21     6 12.72    24  6.60    24  6.60    48  4.70    128  2.89
            // -0.20     6 12.85    24  6.67    24  6.67    48  4.75    128  2.92
            // -0.19     6 12.97    16  8.21    24  6.74    48  4.80    128  2.95
            // -0.18     6 13.10    16  8.29    24  6.81    48  4.85     96  3.44
            // -0.17     6 13.23    16  8.37    24  6.88    48  4.89     96  3.47
            // -0.16     6 13.36    16  8.46    24  6.95    48  4.94     96  3.51
            // -0.15     6 13.49    16  8.54    24  7.02    48  5.00     96  3.54
            // -0.14     6 13.62    16  8.62    24  7.09    48  5.05     96  3.58
            // -0.13     6 13.75    16  8.71    24  7.16    48  5.10     96  3.62
            // -0.12     6 13.88    16  8.80    16  8.80    32  6.28     96  3.65
            // -0.11     6 14.01    12 10.19    16  8.88    32  6.35     96  3.69
            // -0.10     4 16.85    12 10.28    16  8.97    32  6.41     96  3.73
            // -0.09     4 17.01    12 10.38    16  9.06    32  6.47     96  3.76
            // -0.08     4 17.17    12 10.48    16  9.14    32  6.54     96  3.80
            // -0.07     4 17.32    12 10.58    16  9.23    32  6.60     64  4.69
            // -0.06     4 17.48    12 10.69    12 10.69    24  7.67     64  4.74
            // -0.05     4 17.64    12 10.79    12 10.79    24  7.75     64  4.79
            // -0.04     4 17.80    12 10.89    12 10.89    24  7.82     64  4.84
            // -0.03     4 17.96     8 13.26    12 10.99    24  7.90     48  5.63
            // -0.02     4 18.11     8 13.38    12 11.10    24  7.97     48  5.68
            // -0.01     4 18.27     6 15.36     8 13.51    16  9.78     48  5.74
            //  0.00     4  1.00     4  1.00     4  1.00     4  1.00      4  1.00
            //  0.01     4 18.57     6 15.62     8 13.75    16  9.96     48  5.85
            //  0.02     4 18.70     8 13.86    12 11.51    24  8.28     48  5.91
            //  0.03     4 18.83     8 13.97    12 11.61    24  8.36     48  5.97
            //  0.04     4 18.96    12 11.71    12 11.71    24  8.44     64  5.23
            //  0.05     4 19.09    12 11.81    12 11.81    24  8.52     64  5.28
            //  0.06     4 19.22    12 11.92    12 11.92    24  8.60     64  5.33
            //  0.07     4 19.36    12 12.02    16 10.52    32  7.55     64  5.39
            //  0.08     4 19.49    12 12.13    16 10.61    32  7.63     64  5.44
            //  0.09     4 19.62    12 12.23    16 10.71    32  7.70     96  4.50
            //  0.10     4 19.76    12 12.34    16 10.80    32  7.77     96  4.54
            //  0.11     4 19.89    12 12.45    16 10.90    32  7.85     96  4.59
            //  0.12     6 17.01    16 11.00    16 11.00    32  7.92     96  4.63
            //  0.13     6 17.14    16 11.10    16 11.10    32  8.00     96  4.68
            //  0.14     6 17.27    16 11.20    24  9.26    48  6.64     96  4.73
            //  0.15     6 17.40    16 11.30    24  9.35    48  6.70     96  4.77
            //  0.16     6 17.53    16 11.40    24  9.44    48  6.77     96  4.82
            //  0.17     6 17.67    16 11.51    24  9.53    48  6.84     96  4.87
            //  0.18     6 17.80    16 11.61    24  9.62    48  6.90     96  4.92
            //  0.19     6 17.94    16 11.72    24  9.71    48  6.97    128  4.31
            //  0.20     6 18.07    16 11.83    24  9.80    48  7.04    128  4.36
            //  0.21     6 18.21    24  9.90    24  9.90    48  7.11    128  4.40
            //  0.22     6 18.35    24  9.99    24  9.99    48  7.18    128  4.45
            //  0.23     6 18.49    24 10.09    24 10.09    48  7.26    128  4.49
            //  0.24     6 18.63    24 10.19    24 10.19    48  7.33    128  4.54
            //  0.25     6 18.77    24 10.28    24 10.28    48  7.41    128  4.59
            //  0.26     6 18.92    24 10.39    24 10.39    48  7.48    128  4.64
            //  0.27     8 17.00    24 10.49    32  9.17    64  6.58    128  4.69
            //  0.28     8 17.14    24 10.59    32  9.26    64  6.65    128  4.74
            //  0.29     8 17.28    24 10.69    32  9.35    64  6.72    192  3.92
            //  0.30     8 17.43    24 10.80    32  9.45    64  6.79    192  3.96
            //  0.31     8 17.57    24 10.91    32  9.54    64  6.86    192  4.01
            //  0.32     8 17.72    24 11.02    32  9.64    64  6.93    192  4.05
            //  0.33     8 17.87    24 11.13    32  9.74    64  7.01    192  4.09
            //  0.34     8 18.02    24 11.24    32  9.84    64  7.08    192  4.14
            //  0.35     8 18.17    24 11.36    32  9.95    64  7.16    192  4.19
            //  0.36     8 18.32    24 11.47    32 10.05    64  7.24    192  4.23
            //  0.37     8 18.48    32 10.16    32 10.16    64  7.32    192  4.28
            //  0.38     8 18.63    32 10.27    32 10.27    96  6.08    192  4.33
            //  0.39     8 18.79    32 10.38    48  8.58    96  6.15    192  4.38
            //  0.40     8 18.95    32 10.49    48  8.68    96  6.22    192  4.43
            //  0.41     8 19.11    32 10.60    48  8.78    96  6.30    192  4.49
            //  0.42    12 16.45    32 10.72    48  8.88    96  6.37    192  4.54
            //  0.43    12 16.61    32 10.84    48  8.98    96  6.45    192  4.59
            //  0.44    12 16.77    32 10.96    48  9.08    96  6.52    256  4.04
            //  0.45    12 16.93    32 11.09    48  9.19    96  6.60    256  4.09
            //  0.46    12 17.10    32 11.21    48  9.30    96  6.68    256  4.14
            //  0.47    12 17.26    32 11.34    48  9.41    96  6.77    256  4.19
            //  0.48    12 17.44    32 11.47    48  9.52    96  6.85    256  4.24
            //  0.49    12 17.61    48  9.64    48  9.64    96  6.94    256  4.30
            //  0.50    12 17.79    48  9.76    48  9.76    96  7.03    256  4.35
            //  0.51    12 17.97    48  9.88    48  9.88    96  7.12    256  4.41
            //  0.52    12 18.15    48 10.00    48 10.00    96  7.21    256  4.47
            //  0.53    12 18.34    48 10.13    48 10.13   128  6.36    256  4.53
            //  0.54    12 18.53    48 10.26    48 10.26   128  6.45    256  4.59
            //  0.55    12 18.72    48 10.40    64  9.10   128  6.53    384  3.82
            //  0.56    12 18.92    48 10.53    64  9.22   128  6.63    384  3.87
            //  0.57    12 19.12    48 10.68    64  9.35   128  6.72    384  3.93
            //  0.58    12 19.33    48 10.82    64  9.48   128  6.82    384  3.98
            //  0.59    12 19.54    48 10.97    64  9.61   128  6.92    384  4.05
            //  0.60    12 19.75    48 11.13    64  9.75   128  7.02    384  4.11
            //  0.61    12 19.97    48 11.28    64  9.89   128  7.13    384  4.17
            //  0.62    12 20.20    48 11.45    64 10.04   128  7.24    384  4.24
            //  0.63    16 18.39    48 11.62    64 10.19   192  6.05    384  4.31
            //  0.64    16 18.62    64 10.35    64 10.35   192  6.15    384  4.38
            //  0.65    16 18.86    64 10.51    96  8.71   192  6.25    384  4.45
            //  0.66    16 19.10    64 10.68    96  8.85   192  6.36    384  4.53
            //  0.67    16 19.35    64 10.86    96  9.00   192  6.47    512  4.00
            //  0.68    16 19.60    64 11.04    96  9.16   192  6.58    512  4.08
            //  0.69    16 19.87    64 11.23    96  9.32   192  6.70    512  4.15
            //  0.70    16 20.14    64 11.43    96  9.49   192  6.83    512  4.23
            //  0.71    16 20.42    64 11.63    96  9.67   192  6.96    512  4.32
            //  0.72    16 20.71    64 11.85    96  9.85   192  7.10    512  4.40
            //  0.73    16 21.01    96 10.04    96 10.04   192  7.25    512  4.50
            //  0.74    16 21.32    96 10.25    96 10.25   256  6.44    768  3.76
            //  0.75    16 21.65    96 10.46    96 10.46   256  6.58    768  3.84
            //  0.76    16 21.99    96 10.68   128  9.36   256  6.73    768  3.93
            //  0.77    24 19.41    96 10.92   128  9.57   256  6.89    768  4.03
            //  0.78    24 19.76    96 11.17   128  9.79   256  7.06    768  4.13
            //  0.79    24 20.13    96 11.44   128 10.03   256  7.24    768  4.24
            //  0.80    24 20.51    96 11.72   128 10.29   384  6.11    768  4.35
            //  0.81    24 20.92    96 12.02   128 10.56   384  6.28    768  4.48
            //  0.82    24 21.35    96 12.34   192  8.99   384  6.46   1024  4.00
            //  0.83    24 21.81   128 11.16   192  9.26   384  6.66   1024  4.13
            //  0.84    24 22.29   128 11.50   192  9.55   384  6.88   1024  4.26
            //  0.85    24 22.82   128 11.86   192  9.87   384  7.11   1024  4.41
            //  0.86    24 23.38   128 12.26   192 10.21   384  7.37   1024  4.58
            //  0.87    24 24.00   128 12.70   192 10.59   512  6.67   1536  3.90
            //  0.88    24 24.67   192 11.01   192 11.01   512  6.95   1536  4.06
            //  0.89    24 25.41   192 11.48   256 10.07   512  7.26   1536  4.25
            //  0.90    24 26.24   192 12.00   256 10.54   768  6.27   1536  4.47
            //  0.91    24 27.17   192 12.61   256 11.09   768  6.62   2048  4.10
            //  0.92    24 28.23   192 13.30   384  9.74   768  7.02   2048  4.35
            //  0.93    24 29.45   256 12.46   384 10.38   768  7.50   3072  3.82
            //  0.94    24 30.86   256 13.36   384 11.16  1024  7.05   3072  4.13
            //  0.95    24 32.53   384 12.14   512 10.67  1024  7.72   3072  4.53
            //  0.96    24 34.51   384 13.42   512 11.83  1536  7.09   4096  4.40
            //  0.97    24 36.88   512 13.45   768 11.24  2048  7.11   6144  4.16
            //  0.98    16 41.78   768 13.48  1024 11.88  3072  7.12   8192  4.42
            //  0.99     8 44.82  1024 16.00  1536 13.51  6144  7.14  16384  4.43
            // Encode N as small integer: 2,3,4,6,8,12... -> 0,1,2,3,4,5...
            // using this awk script
            //
            // {
            //   n = $1;
            //   if (n % 3 == 0) {
            //     s = 1;
            //     n = n/3;
            //   } else {
            //     s = 0;
            //     n = n/2;
            //   }
            //   p = int( log(n)/log(2)+0.5 );
            //   printf "%d\n", 2*p+s;
            // }
            //
            // A couple of changes have been made: (1) the decrease in N for float and
            // n > 0.97 has been removed; (2) entrys of n=+/-1 have been included
            // (incrementing the previous code value by 1).


            var n = ndiv * _n;         // n in [-ndiv, ndiv]
            int j = ndiv + (int)(n < 0 ? Floor(n) : Ceiling(n)); // j in [0, 2*ndiv]
            int N = narr.Span[j];
            // Decode 0,1,2,3,4,5... -> 2,3,4,6,8,12...
            N = (N % 2 == 0 ? 2 : 3) * (1 << (N / 2));

            _fft = new DST(N);
            _nC4 = N;
        }

        /// <summary>
        /// A global instantiation of <see cref="GeodesicExact"/> with the parameters for the WGS84 ellipsoid.
        /// </summary>
        public static GeodesicExact WGS84 { get; } = new GeodesicExact(Ellipsoid.WGS84);

        /// <inheritdoc/>
        public override double EquatorialRadius => _a;

        /// <inheritdoc/>
        public override double Flattening => _f;

        /// <inheritdoc/>
        public override double EllipsoidArea => 4 * PI * _c2;

        #region Private methods

        private static double Astroid(double x, double y)
        {
            // Solve k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0 for positive root k.
            // This solution is adapted from Geocentric::Reverse.
            double k;
            double
              p = Sq(x),
              q = Sq(y),
              r = (p + q - 1) / 6;
            if (!(q == 0 && r <= 0))
            {
                double
                  // Avoid possible division by zero when r = 0 by multiplying equations
                  // for s and t by r^3 and r, resp.
                  S = p * q / 4,            // S = r^3 * s
                  r2 = Sq(r),
                  r3 = r * r2,
                  // The discriminant of the quadratic equation for T3.  This is zero on
                  // the evolute curve p^(1/3)+q^(1/3) = 1
                  disc = S * (S + 2 * r3);
                var u = r;
                if (disc >= 0)
                {
                    var T3 = S + r3;
                    // Pick the sign on the sqrt to maximize abs(T3).  This minimizes loss
                    // of precision due to cancellation.  The result is unchanged because
                    // of the way the T is used in definition of u.
                    T3 += T3 < 0 ? -Sqrt(disc) : Sqrt(disc); // T3 = (r * t)^3
                                                             // N.B. cbrt always returns the real root.  cbrt(-8) = -2.
                    var T = Cbrt(T3); // T = r * t
                                      // T can be zero; but then r2 / T -> 0.
                    u += T + (T != 0 ? r2 / T : 0);
                }
                else
                {
                    // T is complex, but the way u is defined the result is real.
                    var ang = Atan2(Sqrt(-disc), -(S + r3));
                    // There are three possible cube roots.  We choose the root which
                    // avoids cancellation.  Note that disc < 0 implies that r < 0.
                    u += 2 * r * Cos(ang / 3);
                }
                double
                  v = Sqrt(Sq(u) + q),    // guaranteed positive
                                          // Avoid loss of accuracy when u < 0.
                  uv = u < 0 ? q / (v - u) : u + v, // u+v, guaranteed positive
                  w = (uv - q) / (2 * v);           // positive?
                                                    // Rearrange expression for k to avoid loss of accuracy due to
                                                    // subtraction.  Division by 0 not possible because uv > 0, w >= 0.
                k = uv / (Sqrt(uv + Sq(w)) + w);   // guaranteed positive
            }
            else
            {               // q == 0 && r <= 0
                            // y = 0 with |x| <= 1.  Handle this case directly.
                            // for y small, positive root is k = abs(y)/sqrt(1-x^2)
                k = 0;
            }
            return k;
        }

        private void Lengths(EllipticFunction E,
                 double sig12,
                 double ssig1, double csig1, double dn1,
                 double ssig2, double csig2, double dn2,
                 double cbet1, double cbet2, GeodesicFlags outmask,
                 out double s12b, out double m12b, out double m0,
                 out double M12, out double M21)
        {
            // Return m12b = (reduced length)/_b; also calculate s12b = distance/_b,
            // and m0 = coefficient of secular term in expression for reduced length.
            s12b = m12b = m0 = M12 = M21 = double.NaN;

            outmask = outmask.Flags();
            // outmask & DISTANCE: set s12b
            // outmask & REDUCEDLENGTH: set m12b & m0
            // outmask & GEODESICSCALE: set M12 & M21

            // It's OK to have repeated dummy arguments,
            // e.g., s12b = m0 = M12 = M21 = dummy

            if (outmask.HasAny(GeodesicFlags.Distance))
                // Missing a factor of _b
                s12b = E.E() / (PI / 2) *
                  (sig12 + (E.DeltaE(ssig2, csig2, dn2) - E.DeltaE(ssig1, csig1, dn1)));
            if (outmask.HasAny(GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale))
            {
                double
                  m0x = -E.K2 * E.D() / (PI / 2),
                  J12 = m0x *
                  (sig12 + (E.DeltaD(ssig2, csig2, dn2) - E.DeltaD(ssig1, csig1, dn1)));
                if (outmask.HasAny(GeodesicFlags.ReducedLength))
                {
                    m0 = m0x;
                    // Missing a factor of _b.  Add parens around (csig1 * ssig2) and
                    // (ssig1 * csig2) to ensure accurate cancellation in the case of
                    // coincident points.
                    m12b = dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) -
                      csig1 * csig2 * J12;
                }
                if (outmask.HasAny(GeodesicFlags.GeodesicScale))
                {
                    var csig12 = csig1 * csig2 + ssig1 * ssig2;
                    var t = _ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2);
                    M12 = csig12 + (t * ssig2 - csig2 * J12) * ssig1 / dn1;
                    M21 = csig12 - (t * ssig1 - csig1 * J12) * ssig2 / dn2;
                }
            }
        }

        private double InverseStart(EllipticFunction E,
                  double sbet1, double cbet1, double dn1,
                  double sbet2, double cbet2, double dn2,
                  double lam12, double slam12, double clam12,
                  out double salp1, out double calp1,
                  out double salp2, out double calp2, out double dnm)
        {
            salp2 = calp2 = dnm = double.NaN;
            // Return a starting point for Newton's method in salp1 and calp1 (function
            // value is -1).  If Newton's method doesn't need to be used, return also
            // salp2 and calp2 and function value is sig12.
            double
              sig12 = -1,               // Return value
                                        // bet12 = bet2 - bet1 in [0, pi); bet12a = bet2 + bet1 in (-pi, 0]
              sbet12 = sbet2 * cbet1 - cbet2 * sbet1,
              cbet12 = cbet2 * cbet1 + sbet2 * sbet1;
            var sbet12a = sbet2 * cbet1 + cbet2 * sbet1;
            bool shortline = cbet12 >= 0 && sbet12 < 0.5 &&
              cbet2 * lam12 < 0.5;
            double somg12, comg12;
            if (shortline)
            {
                var sbetm2 = Sq(sbet1 + sbet2);
                // sin((bet1+bet2)/2)^2
                // =  (sbet1 + sbet2)^2 / ((sbet1 + sbet2)^2 + (cbet1 + cbet2)^2)
                sbetm2 /= sbetm2 + Sq(cbet1 + cbet2);
                dnm = Sqrt(1 + _ep2 * sbetm2);
                var omg12 = lam12 / (_f1 * dnm);
                somg12 = Sin(omg12); comg12 = Cos(omg12);
            }
            else
            {
                somg12 = slam12; comg12 = clam12;
            }

            salp1 = cbet2 * somg12;
            calp1 = comg12 >= 0 ?
              sbet12 + cbet2 * sbet1 * Sq(somg12) / (1 + comg12) :
              sbet12a - cbet2 * sbet1 * Sq(somg12) / (1 - comg12);

            double
              ssig12 = Hypot(salp1, calp1),
              csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12;

            if (shortline && ssig12 < _etol2)
            {
                // really short lines
                salp2 = cbet1 * somg12;
                calp2 = sbet12 - cbet1 * sbet2 *
                  (comg12 >= 0 ? Sq(somg12) / (1 + comg12) : 1 - comg12);
                Norm(ref salp2, ref calp2);
                // Set return value
                sig12 = Atan2(ssig12, csig12);
            }
            else if (Abs(_n) > 0.1 || // Skip astroid calc if too eccentric
                     csig12 >= 0 ||
                     ssig12 >= 6 * Abs(_n) * PI * Sq(cbet1))
            {
                // Nothing to do, zeroth order spherical approximation is OK
            }
            else
            {
                // Scale lam12 and bet2 to x, y coordinate system where antipodal point
                // is at origin and singular point is at y = 0, x = -1.
                double x, y, lamscale, betscale;

                var lam12x = Atan2(-slam12, -clam12); // lam12 - pi
                if (_f >= 0)
                {            // In fact f == 0 does not get here
                             // x = dlong, y = dlat
                    {
                        var k2 = Sq(sbet1) * _ep2;
                        E.Reset(-k2, -_ep2, 1 + k2, 1 + _ep2);
                        lamscale = _e2 / _f1 * cbet1 * 2 * E.H();
                    }
                    betscale = lamscale * cbet1;

                    x = lam12x / lamscale;
                    y = sbet12a / betscale;
                }
                else
                {                  // _f < 0
                                   // x = dlat, y = dlong
                    double
                      cbet12a = cbet2 * cbet1 - sbet2 * sbet1,
                      bet12a = Atan2(sbet12a, cbet12a);
                    // In the case of lon12 = 180, this repeats a calculation made in
                    // Inverse.
                    Lengths(E, PI + bet12a,
                            sbet1, -cbet1, dn1, sbet2, cbet2, dn2,
                            cbet1, cbet2, GeodesicFlags.ReducedLength, out _, out var m12b, out var m0, out _, out _);
                    x = -1 + m12b / (cbet1 * cbet2 * m0 * PI);
                    betscale = x < -0.01 ? sbet12a / x :
                      -_f * Sq(cbet1) * PI;
                    lamscale = betscale / cbet1;
                    y = lam12x / lamscale;
                }

                if (y > -tol1_ && x > -1 - xthresh_)
                {
                    // strip near cut
                    // Need real(x) here to cast away the volatility of x for min/max
                    if (_f >= 0)
                    {
                        salp1 = Min(1, -x); calp1 = -Sqrt(1 - Sq(salp1));
                    }
                    else
                    {
                        calp1 = Max(x > -tol1_ ? 0d : -1d, x);
                        salp1 = Sqrt(1 - Sq(calp1));
                    }
                }
                else
                {
                    // Estimate alp1, by solving the astroid problem.
                    //
                    // Could estimate alpha1 = theta + pi/2, directly, i.e.,
                    //   calp1 = y/k; salp1 = -x/(1+k);  for _f >= 0
                    //   calp1 = x/(1+k); salp1 = -y/k;  for _f < 0 (need to check)
                    //
                    // However, it's better to estimate omg12 from astroid and use
                    // spherical formula to compute alp1.  This reduces the mean number of
                    // Newton iterations for astroid cases from 2.24 (min 0, max 6) to 2.12
                    // (min 0 max 5).  The changes in the number of iterations are as
                    // follows:
                    //
                    // change percent
                    //    1       5
                    //    0      78
                    //   -1      16
                    //   -2       0.6
                    //   -3       0.04
                    //   -4       0.002
                    //
                    // The histogram of iterations is (m = number of iterations estimating
                    // alp1 directly, n = number of iterations estimating via omg12, total
                    // number of trials = 148605):
                    //
                    //  iter    m      n
                    //    0   148    186
                    //    1 13046  13845
                    //    2 93315 102225
                    //    3 36189  32341
                    //    4  5396      7
                    //    5   455      1
                    //    6    56      0
                    //
                    // Because omg12 is near pi, estimate work with omg12a = pi - omg12
                    var k = Astroid(x, y);
                    var
                      omg12a = lamscale * (_f >= 0 ? -x * k / (1 + k) : -y * (1 + k) / k);
                    somg12 = Sin(omg12a); comg12 = -Cos(omg12a);
                    // Update spherical estimate of alp1 using omg12 instead of lam12
                    salp1 = cbet2 * somg12;
                    calp1 = sbet12a - cbet2 * sbet1 * Sq(somg12) / (1 - comg12);
                }
            }
            // Sanity check on starting guess.  Backwards check allows NaN through.
            if (!(salp1 <= 0))
                Norm(ref salp1, ref calp1);
            else
            {
                salp1 = 1; calp1 = 0;
            }
            return sig12;
        }

        private double Lambda12(double sbet1, double cbet1, double dn1,
                      double sbet2, double cbet2, double dn2,
                      double salp1, double calp1, double slam120, double clam120,
                      out double salp2, out double calp2, out double sig12,
                      out double ssig1, out double csig1, out double ssig2, out double csig2,
                      EllipticFunction E,
                      out double domg12, bool diffp, out double dlam12)
        {
            dlam12 = double.NaN;

            if (sbet1 == 0 && calp1 == 0)
                // Break degeneracy of equatorial line.  This case has already been
                // handled.
                calp1 = -tiny_;

            double
              // sin(alp1) * cos(bet1) = sin(alp0)
              salp0 = salp1 * cbet1,
              calp0 = Hypot(calp1, salp1 * sbet1); // calp0 > 0

            double somg1, comg1, somg2, comg2, somg12, comg12, cchi1, cchi2, lam12;
            // tan(bet1) = tan(sig1) * cos(alp1)
            // tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1)
            ssig1 = sbet1; somg1 = salp0 * sbet1;
            csig1 = comg1 = calp1 * cbet1;
            // Without normalization we have schi1 = somg1.
            cchi1 = _f1 * dn1 * comg1;
            Norm(ref ssig1, ref csig1);
            // Math::norm(somg1, comg1); -- don't need to normalize!
            // Math::norm(schi1, cchi1); -- don't need to normalize!

            // Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
            // about this case, since this can yield singularities in the Newton
            // iteration.
            // sin(alp2) * cos(bet2) = sin(alp0)
            salp2 = cbet2 != cbet1 ? salp0 / cbet2 : salp1;
            // calp2 = sqrt(1 - sq(salp2))
            //       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
            // and subst for calp0 and rearrange to give (choose positive sqrt
            // to give alp2 in [0, pi/2]).
            calp2 = cbet2 != cbet1 || Abs(sbet2) != -sbet1 ?
              Sqrt(Sq(calp1 * cbet1) +
                   (cbet1 < -sbet1 ?
                    (cbet2 - cbet1) * (cbet1 + cbet2) :
                    (sbet1 - sbet2) * (sbet1 + sbet2))) / cbet2 :
              Abs(calp1);
            // tan(bet2) = tan(sig2) * cos(alp2)
            // tan(omg2) = sin(alp0) * tan(sig2).
            ssig2 = sbet2; somg2 = salp0 * sbet2;
            csig2 = comg2 = calp2 * cbet2;
            // Without normalization we have schi2 = somg2.
            cchi2 = _f1 * dn2 * comg2;
            Norm(ref ssig2, ref csig2);
            // Math::norm(somg2, comg2); -- don't need to normalize!
            // Math::norm(schi2, cchi2); -- don't need to normalize!

            // sig12 = sig2 - sig1, limit to [0, pi]
            sig12 = Atan2(Max(0, csig1 * ssig2 - ssig1 * csig2),
                                       csig1 * csig2 + ssig1 * ssig2);

            // omg12 = omg2 - omg1, limit to [0, pi]
            somg12 = Max(0, comg1 * somg2 - somg1 * comg2);
            comg12 = comg1 * comg2 + somg1 * somg2;
            var k2 = Sq(calp0) * _ep2;
            E.Reset(-k2, -_ep2, 1 + k2, 1 + _ep2);
            // chi12 = chi2 - chi1, limit to [0, pi]
            double
              schi12 = Max(0, cchi1 * somg2 - somg1 * cchi2),
              cchi12 = cchi1 * cchi2 + somg1 * somg2;
            // eta = chi12 - lam120
            var eta = Atan2(schi12 * clam120 - cchi12 * slam120,
                             cchi12 * clam120 + schi12 * slam120);
            var deta12 = -_e2 / _f1 * salp0 * E.H() / (PI / 2) *
              (sig12 + (E.DeltaH(ssig2, csig2, dn2) - E.DeltaH(ssig1, csig1, dn1)));
            lam12 = eta + deta12;
            // domg12 = deta12 + chi12 - omg12
            domg12 = deta12 + Atan2(schi12 * comg12 - cchi12 * somg12,
                                    cchi12 * comg12 + schi12 * somg12);
            if (diffp)
            {
                if (calp2 == 0)
                    dlam12 = -2 * _f1 * dn1 / sbet1;
                else
                {
                    Lengths(E, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                            cbet1, cbet2, GeodesicFlags.ReducedLength,
                            out _, out dlam12, out _, out _, out _);
                    dlam12 *= _f1 / (calp2 * cbet2);
                }
            }

            return lam12;
        }

        private double GenInverse(double lat1, double lon1, double lat2, double lon2,
                        GeodesicFlags outmask, out double s12,
                        out double salp1, out double calp1, out double salp2, out double calp2,
                        out double m12, out double M12, out double M21, out double S12)
        {
            m12 = M12 = M21 = s12 = S12 = double.NaN;
            salp1 = calp1 = salp2 = calp2 = double.NaN;

            // Compute longitude difference (AngDiff does this carefully).  Result is
            // in [-180, 180] but -180 is only for west-going geodesics.  180 is for
            // east-going and meridional geodesics.
            var lon12 = AngDiff(lon1, lon2, out var lon12s);
            // Make longitude difference positive.
            int lonsign = SignBit(lon12) ? -1 : 1;
            lon12 *= lonsign; lon12s *= lonsign;
            double
              lam12 = lon12 * Degree,
              slam12, clam12;

            // Calculate sincos of lon12 + error (this applies AngRound internally).
            SinCosde(lon12, lon12s, out slam12, out clam12);
            // the supplementary longitude difference
            lon12s = (HD - lon12) - lon12s;

            // If really close to the equator, treat as on equator.
            lat1 = AngRound(LatFix(lat1));
            lat2 = AngRound(LatFix(lat2));
            // Swap points so that point with higher (abs) latitude is point 1
            // If one latitude is a nan, then it becomes lat1.
            int swapp = Abs(lat1) < Abs(lat2) || double.IsNaN(lat2) ? -1 : 1;
            if (swapp < 0)
            {
                lonsign *= -1;
                Swap(ref lat1, ref lat2);
            }
            // Make lat1 <= -0
            int latsign = SignBit(lat1) ? 1 : -1;
            lat1 *= latsign;
            lat2 *= latsign;
            // Now we have
            //
            //     0 <= lon12 <= 180
            //     -90 <= lat1 <= -0
            //     lat1 <= lat2 <= -lat1
            //
            // longsign, swapp, latsign register the transformation to bring the
            // coordinates to this canonical form.  In all cases, 1 means no change was
            // made.  We make these transformations so that there are few cases to
            // check, e.g., on verifying quadrants in atan2.  In addition, this
            // enforces some symmetries in the results returned.

            double s12x = 0, m12x = 0;
            // Initialize for the meridian.  No longitude calculation is done in this
            // case to let the parameter default to 0.
            var E = new EllipticFunction(-_ep2);

            SinCosd(lat1, out var sbet1, out var cbet1); sbet1 *= _f1;
            // Ensure cbet1 = +epsilon at poles; doing the fix on beta means that sig12
            // will be <= 2*tiny for two points at the same pole.
            Norm(ref sbet1, ref cbet1); cbet1 = Max(tiny_, cbet1);

            SinCosd(lat2, out var sbet2, out var cbet2); sbet2 *= _f1;
            // Ensure cbet2 = +epsilon at poles
            Norm(ref sbet2, ref cbet2); cbet2 = Max(tiny_, cbet2);

            // If cbet1 < -sbet1, then cbet2 - cbet1 is a sensitive measure of the
            // |bet1| - |bet2|.  Alternatively (cbet1 >= -sbet1), abs(sbet2) + sbet1 is
            // a better measure.  This logic is used in assigning calp2 in Lambda12.
            // Sometimes these quantities vanish and in that case we force bet2 = +/-
            // bet1 exactly.  An example where is is necessary is the inverse problem
            // 48.522876735459 0 -48.52287673545898293 179.599720456223079643
            // which failed with Visual Studio 10 (Release and Debug)

            if (cbet1 < -sbet1)
            {
                if (cbet2 == cbet1)
                    sbet2 = CopySign(sbet1, sbet2);
            }
            else
            {
                if (Abs(sbet2) == -sbet1)
                    cbet2 = cbet1;
            }

            double
              dn1 = (_f >= 0 ? Sqrt(1 + _ep2 * Sq(sbet1)) :
                     Sqrt(1 - _e2 * Sq(cbet1)) / _f1),
              dn2 = (_f >= 0 ? Sqrt(1 + _ep2 * Sq(sbet2)) :
                     Sqrt(1 - _e2 * Sq(cbet2)) / _f1);

            double a12=double.NaN, sig12;

            bool meridian = lat1 == -QD || slam12 == 0;

            if (meridian)
            {

                // Endpoints are on a single full meridian, so the geodesic might lie on
                // a meridian.

                calp1 = clam12; salp1 = slam12; // Head to the target longitude
                calp2 = 1; salp2 = 0;           // At the target we're heading north

                double
                  // tan(bet) = tan(sig) * cos(alp)
                  ssig1 = sbet1, csig1 = calp1 * cbet1,
                  ssig2 = sbet2, csig2 = calp2 * cbet2;

                // sig12 = sig2 - sig1
                sig12 = Atan2(Max(0, csig1 * ssig2 - ssig1 * csig2),
                                           csig1 * csig2 + ssig1 * ssig2);
                {
                    Lengths(E, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                            cbet1, cbet2, outmask | GeodesicFlags.ReducedLength,
                            out s12x, out m12x, out _, out M12, out M21);
                }
                // Add the check for sig12 since zero length geodesics might yield m12 <
                // 0.  Test case was
                //
                //    echo 20.001 0 20.001 0 | GeodSolve -i
                //
                // In fact, we will have sig12 > pi/2 for meridional geodesic which is
                // not a shortest path.
                if (sig12 < 1 || m12x >= 0)
                {
                    // Need at least 2, to handle 90 0 90 180
                    if (sig12 < 3 * tiny_ ||
                            // Prevent negative s12 or m12 for short lines
                            (sig12 < tol0_ && (s12x < 0 || m12x < 0)))
                        sig12 = m12x = s12x = 0;
                    m12x *= _b;
                    s12x *= _b;
                    a12 = sig12 / Degree;
                }
                else
                    // m12 < 0, i.e., prolate and too close to anti-podal
                    meridian = false;
            }

            // somg12 == 2 marks that it needs to be calculated
            double omg12 = 0, somg12 = 2, comg12 = 0;
            if (!meridian &&
                sbet1 == 0 &&   // and sbet2 == 0
                (_f <= 0 || lon12s >= _f * HD))
            {

                // Geodesic runs along equator
                calp1 = calp2 = 0; salp1 = salp2 = 1;
                s12x = _a * lam12;
                sig12 = omg12 = lam12 / _f1;
                m12x = _b * Sin(sig12);
                if (outmask.HasAny(GeodesicFlags.GeodesicScale))
                    M12 = M21 = Cos(sig12);
                a12 = lon12 / _f1;

            }
            else if (!meridian)
            {

                // Now point1 and point2 belong within a hemisphere bounded by a
                // meridian and geodesic is neither meridional or equatorial.

                // Figure a starting point for Newton's method
                sig12 = InverseStart(E, sbet1, cbet1, dn1, sbet2, cbet2, dn2,
                                     lam12, slam12, clam12,
                                     out salp1, out calp1, out salp2, out calp2, out var dnm);

                if (sig12 >= 0)
                {
                    // Short lines (InverseStart sets salp2, calp2, dnm)
                    s12x = sig12 * _b * dnm;
                    m12x = Sq(dnm) * _b * Sin(sig12 / dnm);
                    if (outmask.HasAny(GeodesicFlags.GeodesicScale))
                        M12 = M21 = Cos(sig12 / dnm);
                    a12 = sig12 / Degree;
                    omg12 = lam12 / (_f1 * dnm);
                }
                else
                {

                    // Newton's method.  This is a straightforward solution of f(alp1) =
                    // lambda12(alp1) - lam12 = 0 with one wrinkle.  f(alp) has exactly one
                    // root in the interval (0, pi) and its derivative is positive at the
                    // root.  Thus f(alp) is positive for alp > alp1 and negative for alp <
                    // alp1.  During the course of the iteration, a range (alp1a, alp1b) is
                    // maintained which brackets the root and with each evaluation of
                    // f(alp) the range is shrunk, if possible.  Newton's method is
                    // restarted whenever the derivative of f is negative (because the new
                    // value of alp1 is then further from the solution) or if the new
                    // estimate of alp1 lies outside (0,pi); in this case, the new starting
                    // guess is taken to be (alp1a + alp1b) / 2.
                    //
                    // initial values to suppress warnings (if loop is executed 0 times)
                    double ssig1 = 0, csig1 = 0, ssig2 = 0, csig2 = 0, domg12 = 0;
                    uint numit = 0;
                    // Bracketing range
                    double salp1a = tiny_, calp1a = 1, salp1b = tiny_, calp1b = -1;
                    for (bool tripn = false, tripb = false; ; ++numit)
                    {
                        // 1/4 meridian = 10e6 m and random input.  max err is estimated max
                        // error in nm (checking solution of inverse problem by direct
                        // solution).  iter is mean and sd of number of iterations
                        //
                        //           max   iter
                        // log2(b/a) err mean  sd
                        //    -7     387 5.33 3.68
                        //    -6     345 5.19 3.43
                        //    -5     269 5.00 3.05
                        //    -4     210 4.76 2.44
                        //    -3     115 4.55 1.87
                        //    -2      69 4.35 1.38
                        //    -1      36 4.05 1.03
                        //     0      15 0.01 0.13
                        //     1      25 5.10 1.53
                        //     2      96 5.61 2.09
                        //     3     318 6.02 2.74
                        //     4     985 6.24 3.22
                        //     5    2352 6.32 3.44
                        //     6    6008 6.30 3.45
                        //     7   19024 6.19 3.30
                        double dv;
                        var v = Lambda12(sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1,
                                          slam12, clam12,
                                          out salp2, out calp2, out sig12, out ssig1, out csig1, out ssig2, out csig2,
                                          E, out domg12, numit < maxit1_, out dv);
                        if (tripb ||
                            // Reversed test to allow escape with NaNs
                            !(Abs(v) >= (tripn ? 8 : 1) * tol0_) ||
                            // Enough bisections to get accurate result
                            numit == maxit2_)
                            break;
                        // Update bracketing values
                        if (v > 0 && (numit > maxit1_ || calp1 / salp1 > calp1b / salp1b))
                        { salp1b = salp1; calp1b = calp1; }
                        else if (v < 0 && (numit > maxit1_ || calp1 / salp1 < calp1a / salp1a))
                        { salp1a = salp1; calp1a = calp1; }
                        if (numit < maxit1_ && dv > 0)
                        {
                            var
                              dalp1 = -v / dv;
                            // |dalp1| < pi test moved earlier because GEOGRAPHICLIB_PRECISION
                            // = 5 can result in dalp1 = 10^(10^8).  Then sin(dalp1) takes ages
                            // (because of the need to do accurate range reduction).
                            if (Abs(dalp1) < PI)
                            {
                                double
                                  sdalp1 = Sin(dalp1), cdalp1 = Cos(dalp1),
                                  nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
                                if (nsalp1 > 0)
                                {
                                    calp1 = calp1 * cdalp1 - salp1 * sdalp1;
                                    salp1 = nsalp1;
                                    Norm(ref salp1, ref calp1);
                                    // In some regimes we don't get quadratic convergence because
                                    // slope -> 0.  So use convergence conditions based on epsilon
                                    // instead of sqrt(epsilon).
                                    tripn = Abs(v) <= 16 * tol0_;
                                    continue;
                                }
                            }
                        }
                        // Either dv was not positive or updated value was outside legal
                        // range.  Use the midpoint of the bracket as the next estimate.
                        // This mechanism is not needed for the WGS84 ellipsoid, but it does
                        // catch problems with more eccentric ellipsoids.  Its efficacy is
                        // such for the WGS84 test set with the starting guess set to alp1 =
                        // 90deg:
                        // the WGS84 test set: mean = 5.21, sd = 3.93, max = 24
                        // WGS84 and random input: mean = 4.74, sd = 0.99
                        salp1 = (salp1a + salp1b) / 2;
                        calp1 = (calp1a + calp1b) / 2;
                        Norm(ref salp1, ref calp1);
                        tripn = false;
                        tripb = (Abs(salp1a - salp1) + (calp1a - calp1) < tolb_ ||
                                 Abs(salp1 - salp1b) + (calp1 - calp1b) < tolb_);
                    }
                    {
                        Lengths(E, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                                cbet1, cbet2, outmask, out s12x, out m12x, out _, out M12, out M21);
                    }
                    m12x *= _b;
                    s12x *= _b;
                    a12 = sig12 / Degree;
                    if (outmask.HasAny(GeodesicFlags.GeodesicScale))
                    {
                        // omg12 = lam12 - domg12
                        double sdomg12 = Sin(domg12), cdomg12 = Cos(domg12);
                        somg12 = slam12 * cdomg12 - clam12 * sdomg12;
                        comg12 = clam12 * cdomg12 + slam12 * sdomg12;
                    }
                }
            }

            if (outmask.HasAny(GeodesicFlags.Distance))
                s12 = 0d + s12x;           // Convert -0 to 0

            if (outmask.HasAny(GeodesicFlags.ReducedLength))
                m12 = 0d + m12x;           // Convert -0 to 0

            if (outmask.HasAny(GeodesicFlags.Area))
            {
                double
                  // From Lambda12: sin(alp1) * cos(bet1) = sin(alp0)
                  salp0 = salp1 * cbet1,
                  calp0 = Hypot(calp1, salp1 * sbet1); // calp0 > 0
                double alp12,
                  // Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0).
                  A4 = Sq(_a) * calp0 * salp0 * _e2;
                if (A4 != 0)
                {
                    double
                      k2 = Sq(calp0) * _ep2,
                      // From Lambda12: tan(bet) = tan(sig) * cos(alp)
                      ssig1 = sbet1, csig1 = calp1 * cbet1,
                      ssig2 = sbet2, csig2 = calp2 * cbet2;
                    Norm(ref ssig1, ref csig1);
                    Norm(ref ssig2, ref csig2);

                    var i4 = new I4Integrand(_ep2, k2);
                    Span<double> C4a = stackalloc double[_nC4];
                    _fft.Transform(ref i4, C4a);

                    double
                      B41 = DST.Integral(ssig1, csig1, C4a),
                      B42 = DST.Integral(ssig2, csig2, C4a);
                    S12 = A4 * (B42 - B41);
                }
                else
                    // Avoid problems with indeterminate sig1, sig2 on equator
                    S12 = 0;

                if (!meridian && somg12 == 2)
                {
                    somg12 = Sin(omg12); comg12 = Cos(omg12);
                }

                if (!meridian &&
                    // omg12 < 3/4 * pi
                    comg12 > -0.7071 &&     // Long difference not too big
                    sbet2 - sbet1 < 1.75)
                { // Lat difference not too big
                  // Use tan(Gamma/2) = tan(omg12/2)
                  // * (tan(bet1/2)+tan(bet2/2))/(1+tan(bet1/2)*tan(bet2/2))
                  // with tan(x/2) = sin(x)/(1+cos(x))
                    double domg12 = 1 + comg12, dbet1 = 1 + cbet1, dbet2 = 1 + cbet2;
                    alp12 = 2 * Atan2(somg12 * (sbet1 * dbet2 + sbet2 * dbet1),
                                       domg12 * (sbet1 * sbet2 + dbet1 * dbet2));
                }
                else
                {
                    // alp12 = alp2 - alp1, used in atan2 so no need to normalize
                    double
                      salp12 = salp2 * calp1 - calp2 * salp1,
                      calp12 = calp2 * calp1 + salp2 * salp1;
                    // The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
                    // salp12 = -0 and alp12 = -180.  However this depends on the sign
                    // being attached to 0 correctly.  The following ensures the correct
                    // behavior.
                    if (salp12 == 0 && calp12 < 0)
                    {
                        salp12 = tiny_ * calp1;
                        calp12 = -1;
                    }
                    alp12 = Atan2(salp12, calp12);
                }
                S12 += _c2 * alp12;
                S12 *= swapp * lonsign * latsign;
                // Convert -0 to 0
                S12 += 0;
            }

            // Convert calp, salp to azimuth accounting for lonsign, swapp, latsign.
            if (swapp < 0)
            {
                Swap(ref salp1, ref salp2);
                Swap(ref calp1, ref calp2);
                if (outmask.HasAny(GeodesicFlags.GeodesicScale))
                    Swap(ref M12, ref M21);
            }

            salp1 *= swapp * lonsign; calp1 *= swapp * latsign;
            salp2 *= swapp * lonsign; calp2 *= swapp * latsign;

            // Returned value in [0, 180]
            return a12;
        }

        #endregion

        /// <inheritdoc/>
        public override double GenInverse(double lat1, double lon1, double lat2, double lon2,
                       GeodesicFlags outmask,
                       out double s12, out double azi1, out double azi2,
                       out double m12, out double M12, out double M21, out double S12)
        {
            outmask = outmask.Flags();
            var a12 = GenInverse(lat1, lon1, lat2, lon2,
                               outmask, out s12, out var salp1, out var calp1, out var salp2, out var calp2,
                               out m12, out M12, out M21, out S12);
            if (outmask.HasAny(GeodesicFlags.Azimuth))
            {
                azi1 = Atan2d(salp1, calp1);
                azi2 = Atan2d(salp2, calp2);
            }
            else
            {
                azi1 = azi2 = double.NaN;
            }

            return a12;
        }

        /// <inheritdoc/>
        public override double GenDirect(double lat1, double lon1, double azi1,
                         bool arcmode, double s12_a12, GeodesicFlags outmask,
                         out double lat2, out double lon2, out double azi2,
                         out double s12, out double m12, out double M12, out double M21,
                         out double S12)
        {
            // Automatically supply DISTANCE_IN if necessary
            if (!arcmode) outmask |= GeodesicFlags.DistanceIn;

            return new GeodesicLineExact(this, lat1, lon1, azi1, outmask)
              .                         // Note the dot!
              GenPosition(arcmode, s12_a12, outmask,
                          out lat2, out lon2, out azi2, out s12, out m12, out M12, out M21, out S12);
        }

        /// <inheritdoc/>
        public override IGeodesicLine GenDirectLine(double lat1, double lon1, double azi1, bool arcmode, double s12_a12, GeodesicFlags caps = GeodesicFlags.All)
        {
            azi1 = AngNormalize(azi1);

            // Guard against underflow in salp0.  Also -0 is converted to +0.
            SinCosd(AngRound(azi1), out var salp1, out var calp1);

            // Automatically supply DISTANCE_IN if necessary
            if (!arcmode) caps |= GeodesicFlags.DistanceIn;

            return new GeodesicLineExact(this, lat1, lon1, azi1, salp1, calp1,
                                     caps, arcmode, s12_a12);
        }

        /// <inheritdoc/>
        public override IGeodesicLine Line(double lat1, double lon1, double azi1, GeodesicFlags caps = GeodesicFlags.All)
            => new GeodesicLineExact(this, lat1, lon1, azi1, caps);

        /// <inheritdoc/>
        public override IGeodesicLine InverseLine(double lat1, double lon1, double lat2, double lon2, GeodesicFlags caps = GeodesicFlags.All)
        {
            double  
                a12 = GenInverse(lat1, lon1, lat2, lon2,
                   // No need to specify AZIMUTH here
                   0u,out _, out var salp1, out var calp1, out var salp2, out var calp2,
                   out _, out _, out _, out _),
                azi1 = Atan2d(salp1, calp1);
            // Ensure that a12 can be converted to a distance
            if (caps.Flags().HasAny(GeodesicFlags.DistanceIn)) caps |= GeodesicFlags.Distance;
            return new GeodesicLineExact(this, lat1, lon1, azi1, salp1, calp1, caps, true, a12);
        }

        internal readonly struct I4Integrand
        {
            private readonly double X, tX, tdX, sX, sX1, sXX1, asinhsX, _k2;

            public I4Integrand(double ep2, double k2)
            {
                X = ep2;
                tX = t(X);
                tdX = td(X);
                _k2 = k2;

                sX = Sqrt(Abs(X));     // ep
                sX1 = Sqrt(1 + X);     // 1/(1-f)
                sXX1 = sX * sX1;
                asinhsX = X > 0 ? Asinh(sX) : Asin(sX);
            }

            /// <summary>
            /// return asinh(sqrt(x))/sqrt(x)
            /// </summary>
            /// <param name="x"></param>
            /// <returns></returns>
            private static double asinhsqrt(double x)
                => x == 0 ? 1 :
                      (x > 0 ? Asinh(Sqrt(x)) / Sqrt(x) :
                       Asin(Sqrt(-x)) / Sqrt(-x)); // NaNs end up here

            /// <summary>
            /// This differs by from t as defined following Eq 61 in Karney (2013) by
            /// the final subtraction of 1.  This changes nothing since Eq 61 uses the
            /// difference of two evaluations of t and improves the accuracy(?).
            /// </summary>
            /// <param name="x"></param>
            /// <returns></returns>
            private static double t(double x)
                // Group terms to minimize roundoff
                // with x = ep2, this is the same as
                // e2/(1-e2) + (atanh(e)/e - 1)
                => x + (Sqrt(1 + x) * asinhsqrt(x) - 1);

            /// <summary>
            /// d t(x) / dx
            /// </summary>
            /// <param name="x"></param>
            /// <returns></returns>
            private static double td(double x)
                => x == 0 ? 4 / 3d :
                  // Group terms to minimize roundoff
                  1 + (1 - asinhsqrt(x) / Sqrt(1 + x)) / (2 * x);

            private double DtX(double y)
            {
                if (X == y) return tdX;
                if (X * y <= 0) return (tX - t(y)) / (X - y);
                double
                  sy = Sqrt(Abs(y)), sy1 = Sqrt(1 + y),
                  z = (X - y) / (sX * sy1 + sy * sX1),
                  d1 = 2 * sX * sy,
                  d2 = 2 * (X * sy * sy1 + y * sXX1);
                return X > 0 ?
                  (1 + (Asinh(z) / z) / d1 - (asinhsX + Asinh(sy)) / d2) :
                  // NaNs fall through to here
                  (1 - (Asin(z) / z) / d1 - (asinhsX + Asin(sy)) / d2);
            }

            public double Evaluate(double sig)
            {
                var ssig = Sin(sig);
                return -DtX(_k2 * Sq(ssig)) * ssig / 2;
            }
        }
    }
}
