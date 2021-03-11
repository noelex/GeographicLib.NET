using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

using static System.Math;
using static GeographicLib.MathEx;

namespace GeographicLib
{
    /// <summary>
    /// An accumulator for sums.
    /// </summary>
    /// <remarks>
    /// This allows many numbers of floating point type T to be added together with twice the normal precision.
    /// Thus if T is double, the effective precision of the sum is 106 bits or about 32 decimal places.
    /// <para>
    /// The implementation follows J. R. Shewchuk,
    /// <a href="https://doi.org/10.1007/PL00009321">Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates</a>,
    /// Discrete &amp; Computational Geometry 18(3) 305–363 (1997).
    /// </para>
    /// Approximate timings (summing a vector&lt;double>)
    /// <list type="bullet">
    /// <item>double: 2ns</item>
    /// <item>Accumulator&lt;double>: 23ns</item>
    /// </list>
    /// In the documentation of the member functions, <i>sum</i> stands for the value currently held in the accumulator.
    /// </remarks>
    public class Accumulator
    {
        private double _s, _t;

        /// <summary>
        /// Construct from a T. This is not declared explicit, so that you can write Accumulator a = 5;
        /// </summary>
        /// <param name="y">set <i>sum</i> = <i>y</i>.</param>
        public Accumulator(double y = 0)
        {
            (_s, _t) = (y, 0);
        }

        /// <summary>
        /// Add a number to the accumulator.
        /// </summary>
        /// <param name="acc"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static Accumulator operator +(Accumulator acc, double y)
        {
            acc.Add(y);
            return acc;
        }

        /// <summary>
        /// Subtract a number from the accumulator.
        /// </summary>
        /// <param name="acc"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static Accumulator operator -(Accumulator acc, double y)
        {
            acc.Add(-y);
            return acc;
        }

        /// <summary>
        /// Multiply accumulator by an integer.
        /// To avoid loss of accuracy, use only integers such that <i>n</i> × <i>T</i> is exactly representable as a <i>T</i> (i.e., ± powers of two).
        /// Use <i>n</i> = −1 to negate sum.
        /// </summary>
        /// <param name="acc"></param>
        /// <param name="n">set <i>sum</i> *= <i>n</i>.</param>
        /// <returns></returns>
        public static Accumulator operator *(Accumulator acc, int n)
        {
            acc._s *= n;
            acc._t *= n;
            return acc;
        }

        /// <summary>
        /// Multiply accumulator by a number.
        /// The fma (fused multiply and add) instruction is used (if available) in order to maintain accuracy.
        /// </summary>
        /// <param name="acc"></param>
        /// <param name="y">set <i>sum</i> *= <i>y</i>.</param>
        /// <returns></returns>
        public static Accumulator operator *(Accumulator acc, double y)
        {
            var d = acc._s; acc._s *= y;

            d = FusedMultiplyAdd(y, d, -acc._s); // the error in the first multiplication
            acc._t = FusedMultiplyAdd(y, acc._t, d); // add error to the second term

            return acc;
        }

        /// <summary>
        /// Return the value held in the accumulator.
        /// </summary>
        /// <param name="acc"></param>
        public static implicit operator double(Accumulator acc) => acc._s;

        /// <summary>
        /// Set the accumulator to a number.
        /// </summary>
        /// <param name="y">set <i>sum</i> = <i>y</i>.</param>
        public static implicit operator Accumulator(double y) => new Accumulator(y);

        /// <summary>
        /// Return the result of adding a number to sum (but don't change sum).
        /// </summary>
        /// <param name="y">the number to be added to the sum.</param>
        /// <returns><i>sum</i> + <i>y</i>.</returns>
        public double Sum(double y)
        {
            Accumulator a = this;
            a.Add(y);
            return a._s;
        }

        /// <summary>
        /// Reduce accumulator to the range [-y/2, y/2].
        /// </summary>
        /// <param name="y">the modulus</param>
        /// <returns></returns>
        public Accumulator Remainder(double y)
        {
            _s = IEEERemainder(_s, y);
            Add(0);
            return this;
        }

        private void Add(double y)
        {
            // Here's Shewchuk's solution...
            // hold exact sum as [s, t, u]
            // Accumulate starting at least significant end

            y = MathEx.Sum(y, _t, out var u);
            _s = MathEx.Sum(y, _s, out _t);
            // Start is _s, _t decreasing and non-adjacent.  Sum is now (s + t + u)
            // exactly with s, t, u non-adjacent and in decreasing order (except for
            // possible zeros).  The following code tries to normalize the result.
            // Ideally, we want _s = round(s+t+u) and _u = round(s+t+u - _s).  The
            // following does an approximate job (and maintains the decreasing
            // non-adjacent property).  Here are two "failures" using 3-bit floats:
            //
            // Case 1: _s is not equal to round(s+t+u) -- off by 1 ulp
            // [12, -1] - 8 -> [4, 0, -1] -> [4, -1] = 3 should be [3, 0] = 3
            //
            // Case 2: _s+_t is not as close to s+t+u as it shold be
            // [64, 5] + 4 -> [64, 8, 1] -> [64,  8] = 72 (off by 1)
            //                    should be [80, -7] = 73 (exact)
            //
            // "Fixing" these problems is probably not worth the expense.  The
            // representation inevitably leads to small errors in the accumulated
            // values.  The additional errors illustrated here amount to 1 ulp of the
            // less significant word during each addition to the Accumulator and an
            // additional possible error of 1 ulp in the reported sum.
            //
            // Incidentally, the "ideal" representation described above is not
            // canonical, because _s = round(_s + _t) may not be true.  For example,
            // with 3-bit floats:
            //
            // [128, 16] + 1 -> [160, -16] -- 160 = round(145).
            // But [160, 0] - 16 -> [128, 16] -- 128 = round(144).
            //
            if (_s == 0)              // This implies t == 0,
                _s = u;                 // so result is u
            else
                _t += u;                // otherwise just accumulate u to t.
        }


    }
}
