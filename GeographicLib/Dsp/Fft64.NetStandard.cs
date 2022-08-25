#if !NET5_0_OR_GREATER
using System;
using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace GeographicLib.Dsp
{
    internal partial class Fft64
    {
        private void FFT(ReadOnlyComplexSpan intput, ComplexSpan output, ComplexSpan dftCoefficient)
        {
            Debug.Assert(intput.Real.Length == windowSize);
            Debug.Assert(intput.Imaginary.Length == windowSize);
            Debug.Assert(output.Real.Length == windowSize);
            Debug.Assert(output.Imaginary.Length == windowSize);

            // Cache variables onto the stack
            int depth = this.depth;

            // We need O(n) extra memory. Make sure the last layer is assigned to the output array.
            // Lower layer will be the final output.
            Span<double> tempReal = stackalloc double[windowSize];
            Span<double> tempImaginary = stackalloc double[windowSize];
            ComplexSpan temp = new ComplexSpan(tempReal, tempImaginary);
            ComplexSpan passInput = depth % 2 == 0 ? output : temp;
            ComplexSpan passOutput = depth % 2 == 0 ? temp : output;

            // Assign the initial pass.
            intput.Real.CopyTo(passInput.Real);
            intput.Imaginary.CopyTo(passInput.Imaginary);

            // Calculate the remaining passes.
            // i.e. calculating the ith pass based on the output of the (i - 1)th path.
            int windowSizeMask = windowSize - 1;
            int vectorCount = Vector<double>.Count;
            Span<double> dftCoefficientReal = dftCoefficient.Real;
            Span<double> dftCoefficientImaginary = dftCoefficient.Imaginary;
            for (int pass = 1; pass <= depth; pass++)
            {
                int halfCapacityOfTopBits = 1 << pass - 1; // 2 ^ (pass - 1)
                int capacityOfBottomBits = 1 << (depth - pass); // 2 ^ (depth - pass)
                int bottomBitsMask = capacityOfBottomBits - 1;
                int topBitsMask = windowSizeMask & ~bottomBitsMask;

                Span<double> passInputReal = passInput.Real;
                Span<double> passInputImaginary = passInput.Imaginary;
                Span<double> passOutputReal = passOutput.Real;
                Span<double> passOutputImaginary = passOutput.Imaginary;
                for (int i = 0; i <= windowSizeMask;) // i will increment in nested loops
                {
                    int topIndex = i >> (depth - pass); // i % 2 ^ (depth - #pass), aka top #pass bits
                    if (capacityOfBottomBits < vectorCount)
                    {
                        if (topIndex < halfCapacityOfTopBits)
                        {
                            int coefficientIndex = topIndex << depth - pass;
                            double oddCoefficientReal = dftCoefficientReal[coefficientIndex];
                            double oddCoefficientImaginary = dftCoefficientImaginary[coefficientIndex];
                            int butterflyIndex = coefficientIndex << 1;
                            for (int bottomIndex = 0; bottomIndex < capacityOfBottomBits; bottomIndex++, i++)
                            {
                                int evenIndex = bottomIndex + butterflyIndex;
                                double evenReal = passInputReal[evenIndex];
                                double evenImaginary = passInputImaginary[evenIndex];
                                int oddIndex = bottomIndex + butterflyIndex + capacityOfBottomBits;
                                double oddReal = passInputReal[oddIndex];
                                double oddImaginary = passInputImaginary[oddIndex];
                                passOutputReal[i] = evenReal + oddCoefficientReal * oddReal - oddCoefficientImaginary * oddImaginary;
                                passOutputImaginary[i] = evenImaginary + oddCoefficientReal * oddImaginary + oddCoefficientImaginary * oddReal;
                            }
                        }
                        else
                        {
                            int coefficientIndex = topIndex - halfCapacityOfTopBits << depth - pass;
                            double oddCoefficientReal = dftCoefficientReal[coefficientIndex];
                            double oddCoefficientImaginary = dftCoefficientImaginary[coefficientIndex];
                            int butterflyIndex = coefficientIndex << 1;
                            for (int bottomIndex = 0; bottomIndex < capacityOfBottomBits; bottomIndex++, i++)
                            {
                                int evenIndex = bottomIndex + butterflyIndex;
                                double evenReal = passInputReal[evenIndex];
                                double evenImaginary = passInputImaginary[evenIndex];
                                int oddIndex = bottomIndex + butterflyIndex + capacityOfBottomBits;
                                double oddReal = passInputReal[oddIndex];
                                double oddImaginary = passInputImaginary[oddIndex];
                                passOutputReal[i] = evenReal - oddCoefficientReal * oddReal + oddCoefficientImaginary * oddImaginary;
                                passOutputImaginary[i] = evenImaginary - oddCoefficientReal * oddImaginary - oddCoefficientImaginary * oddReal;
                            }
                        }
                    }
                    else
                    {
                        if (topIndex < halfCapacityOfTopBits)
                        {
                            int coefficientIndex = topIndex << depth - pass;
                            Vector<double> oddCoefficientReal = new Vector<double>(dftCoefficientReal[coefficientIndex]);
                            Vector<double> oddCoefficientImaginary = new Vector<double>(dftCoefficientImaginary[coefficientIndex]);
                            int butterflyIndex = coefficientIndex << 1;
                            for (int bottomIndex = 0; bottomIndex < capacityOfBottomBits; bottomIndex += vectorCount, i += vectorCount)
                            {
                                int evenIndex = bottomIndex + butterflyIndex;
                                Copy(passInputReal.Slice(evenIndex), out Vector<double> evenReal);
                                Copy(passInputImaginary.Slice(evenIndex), out Vector<double> evenImaginary);

                                int oddIndex = bottomIndex + butterflyIndex + capacityOfBottomBits;
                                Copy(passInputReal.Slice(oddIndex), out Vector<double> oddReal);
                                Copy(passInputImaginary.Slice(oddIndex), out Vector<double> oddImaginary);

                                Copy(evenReal + oddCoefficientReal * oddReal - oddCoefficientImaginary * oddImaginary, passOutputReal.Slice(i));
                                Copy(evenImaginary + oddCoefficientReal * oddImaginary + oddCoefficientImaginary * oddReal, passOutputImaginary.Slice(i));
                            }
                        }
                        else
                        {
                            int coefficientIndex = topIndex - halfCapacityOfTopBits << depth - pass;
                            Vector<double> oddCoefficientReal = new Vector<double>(dftCoefficientReal[coefficientIndex]);
                            Vector<double> oddCoefficientImaginary = new Vector<double>(dftCoefficientImaginary[coefficientIndex]);
                            int butterflyIndex = coefficientIndex << 1;
                            for (int bottomIndex = 0; bottomIndex < capacityOfBottomBits; bottomIndex += vectorCount, i += vectorCount)
                            {
                                int evenIndex = bottomIndex + butterflyIndex;
                                Copy(passInputReal.Slice(evenIndex), out Vector<double> evenReal);
                                Copy(passInputImaginary.Slice(evenIndex), out Vector<double> evenImaginary);

                                int oddIndex = bottomIndex + butterflyIndex + capacityOfBottomBits;
                                Copy(passInputReal.Slice(oddIndex), out Vector<double> oddReal);
                                Copy(passInputImaginary.Slice(oddIndex), out Vector<double> oddImaginary);

                                Copy(evenReal - oddCoefficientReal * oddReal + oddCoefficientImaginary * oddImaginary, passOutputReal.Slice(i));
                                Copy(evenImaginary - oddCoefficientReal * oddImaginary - oddCoefficientImaginary * oddReal, passOutputImaginary.Slice(i));
                            }
                        }
                    }
                }
                ComplexSpan exchange = passOutput;
                passOutput = passInput;
                passInput = exchange;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void Copy<T>(Span<T> src, out Vector<T> dest) where T : struct
#if NETSTANDARD2_0
            => dest = Unsafe.ReadUnaligned<Vector<T>>(ref Unsafe.As<T, byte>(ref MemoryMarshal.GetReference(src)));
#else
            => dest = new Vector<T>(src);
#endif

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void Copy<T>(Vector<T> src, Span<T> dest) where T : struct
            => Unsafe.WriteUnaligned(ref Unsafe.As<T, byte>(ref MemoryMarshal.GetReference(dest)), src);
    }
}
#endif
