#if NET5_0_OR_GREATER
using System;
using System.Diagnostics;
using System.Numerics;
using System.Runtime.Intrinsics.X86;
using System.Runtime.Intrinsics;

namespace GeographicLib.Dsp
{
    internal partial class Fft64
    {
        private void FFT(ReadOnlyComplexSpan intput, ComplexSpan output, ComplexSpan coefficient)
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
            int vectorCount = Vector<double>.Count;
            for (int pass = 1; pass <= depth; pass++)
            {
                if (depth - pass >= 4 && Fma.IsSupported)
                {
                    FftPassFma(pass, passInput, passOutput, coefficient);
                }
                else
                {
                    FftPassScalar(pass, passInput, passOutput, coefficient);
                }
                ComplexSpan exchange = passOutput;
                passOutput = passInput;
                passInput = exchange;
            }
        }

        private unsafe void FftPassFma(int pass, ReadOnlyComplexSpan passInput, ComplexSpan passOutput, ReadOnlyComplexSpan coefficient)
        {
            int windowSizeMask = windowSize - 1;
            int halfCapacityOfTopBits = 1 << pass - 1; // 2 ^ (pass - 1)
            int capacityOfBottomBits = 1 << (depth - pass); // 2 ^ (depth - pass)
            int bottomBitsMask = capacityOfBottomBits - 1;
            int topBitsMask = windowSizeMask & ~bottomBitsMask;
            int vectorCount = Vector256<double>.Count;

            fixed (double*
                passInputReal = passInput.Real,
                passInputImaginary = passInput.Imaginary,
                passOutputReal = passOutput.Real,
                passOutputImaginary = passOutput.Imaginary,
                coefficientReal = coefficient.Real,
                coefficientImaginary = coefficient.Imaginary)
            {
                for (int i = 0; i <= windowSizeMask;) // i will increment in nested loops
                {
                    int topIndex = i >> (depth - pass); // i % 2 ^ (depth - #pass), aka top #pass bits
                    if (topIndex < halfCapacityOfTopBits)
                    {
                        int coefficientIndex = topIndex << depth - pass;
                        double oddCoefficientRealValue = *(coefficientReal + coefficientIndex);
                        double* oddCoefficientRealData = stackalloc double[]
                        {
                            oddCoefficientRealValue,
                            oddCoefficientRealValue,
                            oddCoefficientRealValue,
                            oddCoefficientRealValue
                        };
                        Vector256<double> oddCoefficientReal = Avx2.LoadVector256(oddCoefficientRealData);
                        double oddCoefficientImaginaryValue = *(coefficientImaginary + coefficientIndex);
                        double* oddCoefficientImaginaryData = stackalloc double[]
                        {
                            oddCoefficientImaginaryValue,
                            oddCoefficientImaginaryValue,
                            oddCoefficientImaginaryValue,
                            oddCoefficientImaginaryValue
                        };
                        Vector256<double> oddCoefficientImaginary = Avx2.LoadVector256(oddCoefficientImaginaryData);

                        int butterflyIndex = coefficientIndex << 1;
                        for (int bottomIndex = 0; bottomIndex < capacityOfBottomBits; bottomIndex += vectorCount, i += vectorCount)
                        {
                            int evenIndex = bottomIndex + butterflyIndex;
                            Vector256<double> evenReal = Avx2.LoadVector256(passInputReal + evenIndex);
                            Vector256<double> evenImaginary = Avx2.LoadVector256(passInputImaginary + evenIndex);
                            int oddIndex = bottomIndex + butterflyIndex + capacityOfBottomBits;
                            Vector256<double> oddReal = Avx2.LoadVector256(passInputReal + oddIndex);
                            Vector256<double> oddImaginary = Avx2.LoadVector256(passInputImaginary + oddIndex);
                            Avx2.Store(passOutputReal + i, Fma.MultiplyAdd(oddCoefficientReal, oddReal, Fma.MultiplyAddNegated(oddCoefficientImaginary, oddImaginary, evenReal)));
                            Avx2.Store(passOutputImaginary + i, Fma.MultiplyAdd(oddCoefficientReal, oddImaginary, Fma.MultiplyAdd(oddCoefficientImaginary, oddReal, evenImaginary)));
                        }
                    }
                    else
                    {
                        int coefficientIndex = topIndex - halfCapacityOfTopBits << depth - pass;
                        double oddCoefficientRealValue = *(coefficientReal + coefficientIndex);
                        double* oddCoefficientRealData = stackalloc double[]
                        {
                            oddCoefficientRealValue,
                            oddCoefficientRealValue,
                            oddCoefficientRealValue,
                            oddCoefficientRealValue
                        };
                        Vector256<double> oddCoefficientReal = Avx2.LoadVector256(oddCoefficientRealData);
                        double oddCoefficientImaginaryValue = *(coefficientImaginary + coefficientIndex);
                        double* oddCoefficientImaginaryData = stackalloc double[]
                        {
                            oddCoefficientImaginaryValue,
                            oddCoefficientImaginaryValue,
                            oddCoefficientImaginaryValue,
                            oddCoefficientImaginaryValue
                        };
                        Vector256<double> oddCoefficientImaginary = Avx2.LoadVector256(oddCoefficientImaginaryData);

                        int butterflyIndex = coefficientIndex << 1;
                        for (int bottomIndex = 0; bottomIndex < capacityOfBottomBits; bottomIndex += vectorCount, i += vectorCount)
                        {
                            int evenIndex = bottomIndex + butterflyIndex;
                            Vector256<double> evenReal = Avx2.LoadVector256(passInputReal + evenIndex);
                            Vector256<double> evenImaginary = Avx2.LoadVector256(passInputImaginary + evenIndex);
                            int oddIndex = bottomIndex + butterflyIndex + capacityOfBottomBits;
                            Vector256<double> oddReal = Avx2.LoadVector256(passInputReal + oddIndex);
                            Vector256<double> oddImaginary = Avx2.LoadVector256(passInputImaginary + oddIndex);
                            Avx2.Store(passOutputReal + i, Fma.MultiplyAddNegated(oddCoefficientReal, oddReal, Fma.MultiplyAdd(oddCoefficientImaginary, oddImaginary, evenReal)));
                            Avx2.Store(passOutputImaginary + i, Fma.MultiplyAddNegated(oddCoefficientReal, oddImaginary, Fma.MultiplyAddNegated(oddCoefficientImaginary, oddReal, evenImaginary)));
                        }
                    }
                }
            }
        }

        private void FftPassScalar(int pass, ReadOnlyComplexSpan passInput, ComplexSpan passOutput, ReadOnlyComplexSpan coefficient)
        {
            int windowSizeMask = windowSize - 1;
            int halfCapacityOfTopBits = 1 << pass - 1; // 2 ^ (pass - 1)
            int capacityOfBottomBits = 1 << (depth - pass); // 2 ^ (depth - pass)
            int bottomBitsMask = capacityOfBottomBits - 1;
            int topBitsMask = windowSizeMask & ~bottomBitsMask;

            ReadOnlySpan<double> passInputReal = passInput.Real;
            ReadOnlySpan<double> passInputImaginary = passInput.Imaginary;
            Span<double> passOutputReal = passOutput.Real;
            Span<double> passOutputImaginary = passOutput.Imaginary;
            ReadOnlySpan<double> coefficientReal = coefficient.Real;
            ReadOnlySpan<double> coefficientImaginary = coefficient.Imaginary;
            for (int i = 0; i <= windowSizeMask;) // i will increment in nested loops
            {
                int topIndex = i >> (depth - pass); // i % 2 ^ (depth - #pass), aka top #pass bits
                if (topIndex < halfCapacityOfTopBits)
                {
                    int coefficientIndex = topIndex << depth - pass;
                    double oddCoefficientReal = coefficientReal[coefficientIndex];
                    double oddCoefficientImaginary = coefficientImaginary[coefficientIndex];
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
                    double oddCoefficientReal = coefficientReal[coefficientIndex];
                    double oddCoefficientImaginary = coefficientImaginary[coefficientIndex];
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
        }
    }
}
#endif