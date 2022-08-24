﻿// SPDX-FileCopyrightText: 2018 - 2020 Martin Chloride
// SPDX-License-Identifier: MIT

using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib.Dsp
{
    ref struct ComplexSpan
    {
        public readonly Span<double> Real;
        public readonly Span<double> Imaginary;

        public ComplexSpan(Span<double> real, Span<double> imaginary)
        {
            Real = real;
            Imaginary = imaginary;
        }

        public Complex this[int index]
        {
            get => new Complex(Real[index], Imaginary[index]);
            set
            {
                Real[index] = value.Real;
                Imaginary[index] = value.Imaginary;
            }
        }
    }

    ref struct ReadOnlyComplexSpan
    {
        public readonly ReadOnlySpan<double> Real;
        public readonly ReadOnlySpan<double> Imaginary;

        public ReadOnlyComplexSpan(ReadOnlySpan<double> real, ReadOnlySpan<double> imaginary)
        {
            Real = real;
            Imaginary = imaginary;
        }

        public Complex this[int index]
        {
            get => new Complex(Real[index], Imaginary[index]);
        }

        public static implicit operator ReadOnlyComplexSpan(ComplexSpan rhs)
        {
            return new ReadOnlyComplexSpan(rhs.Real, rhs.Imaginary);
        }
    }
}