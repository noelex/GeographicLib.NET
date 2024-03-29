﻿using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Text;

namespace GeographicLib
{
    /// <summary>
    /// Some utility routines for GeographicLib.NET.
    /// </summary>
    public static class Utility
    {
#if NETSTANDARD2_0
        internal static int IndexOf(this string str, char c, StringComparison stringComparison)
            => str.IndexOf(c.ToStringInvariant(), stringComparison);
#endif
        /// <summary>
        /// Swap byte order of elements in the specified array.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="array"></param>
        internal static void Swab<T>(Span<T> array) where T : struct
        {
            for (int i = 0; i < array.Length; i++)
            {
                MemoryMarshal.Cast<T, byte>(array.Slice(i, 1)).Reverse();
            }
        }

        internal static string ToFixedString(this double x, int p = -1)
        {
            if (!MathEx.IsFinite(x))
                return x < 0 ? "-inf" :
                    (x > 0 ? "inf" : "nan");
#if NETSTANDARD2_0
            // In .NET Framework and .NET Core up to .NET Core 2.0, MidpointRounding.AwayFromZero is the default.
            // See https://docs.microsoft.com/en-us/dotnet/standard/base-types/standard-numeric-format-strings
            if (p >= 0)
            {
                x = Math.Round(x, p, MidpointRounding.ToEven);
            }

            // double.ToString does not output sign for zero in early versions of .NET.
            var sign = x == 0 && MathEx.SignBit(x) ? "-" : "";
            return string.Format("{0}{1}", sign, p >= 0 ? x.ToStringInvariant($"F{p}") : x.ToStringInvariant());
#else
            return p >= 0 ? x.ToStringInvariant($"F{p}") : x.ToStringInvariant();
#endif
        }

        internal static string ToStringInvariant<T>(this T value, string format)
            where T : IFormattable
            => value.ToString(format, CultureInfo.InvariantCulture);

        internal static string ToStringInvariant(this double value)
            => value.ToString(CultureInfo.InvariantCulture);

        internal static string ToStringInvariant(this long value)
            => value.ToString(CultureInfo.InvariantCulture);

        internal static string ToStringInvariant(this int value)
            => value.ToString(CultureInfo.InvariantCulture);

        internal static string ToStringInvariant(this char value)
            => value.ToString(CultureInfo.InvariantCulture);

        internal static double ParseFract(this string s) => ParseFract(s.AsSpan());

        internal static double ParseFract(this ReadOnlySpan<char> s)
        {
            var delim = s.IndexOf('/');
            return
              !(delim != -1 && delim >= 1 && delim + 2 <= s.Length) ?
                s.ParseDoubleInfNan() :
                // delim in [1, size() - 2]
                s.Slice(0, delim).ParseDoubleInfNan() / s.Slice(delim + 1).ParseDoubleInfNan();
        }

        internal static double ParseDoubleInfNan(this string s) => ParseDoubleInfNan(s.AsSpan());

        internal static double ParseDoubleInfNan(this ReadOnlySpan<char> s)
        {
            double x;

            var t = s.Trim();
            var sep = t.IndexOfAny(" \t\n\v\f\r".AsSpan());
            if (sep > 0)
            {
                throw new GeographicException($"Extra text {t.Slice(sep + 1).ToString()} at end of {t.ToString()}");
            }

            if (!t.TryParseDouble(out x) && (x = t.NumMatch()) == 0)
            {
                throw new GeographicException($"Cannot decode {t.ToString()}");
            }

#if NETSTANDARD2_0
            // double.TryParse ignores sign of zero in early versions of .NET.
            if (t[0] == '-')
            {
                x = MathEx.CopySign(x, -1);
            }
#endif
            return x;
        }

        internal static double NumMatch(this string s) => s.AsSpan().NumMatch();

        internal static double NumMatch(this ReadOnlySpan<char> s)
        {
            if (s.Length < 3) return 0;

            Span<char> t = stackalloc char[s.Length];
            s.ToUpper(t, CultureInfo.CurrentCulture);

            var sign = t[0] == '-' ? -1 : 1;
            var p0 = t[0] == '-' || t[0] == '+' ? 1 : 0;
            var p1 = t.FindLastNotOf("0");

            if (p1 == -1 || p1 + 1 < p0 + 3)
                return 0;

            if (t.SequenceEqual("NAN".AsSpan()) || t.SequenceEqual("1.#QNAN".AsSpan()) ||
               t.SequenceEqual("1.#SNAN".AsSpan()) || t.SequenceEqual("1.#IND".AsSpan()) || t.SequenceEqual("1.#R".AsSpan()))
            {
                return double.NaN;
            }

            t = t.Slice(p0);
            if (t.SequenceEqual("INF".AsSpan()) || t.SequenceEqual("1.#INF".AsSpan()) ||
                t.SequenceEqual("∞".AsSpan()) || t.SequenceEqual("INFINITY".AsSpan()))
            {
                return sign * double.PositiveInfinity;
            }

            return 0;
        }

        /// <summary>
        /// Read data of type <see cref="byte"/> from a binary stream to an array of type <typeparamref name="IntT"/>.
        /// <para>
        /// The data in the file is in (<paramref name="bigendp"/> ? big : little)-endian format.
        /// </para>
        /// </summary>
        /// <typeparam name="IntT">the type of the objects in the array (internal).</typeparam>
        /// <param name="stream">the input stream containing the data of type <see cref="byte"/> (external).</param>
        /// <param name="array">the output array of type <typeparamref name="IntT"/> (internal).</param>
        /// <param name="bigendp"><see langword="true"/> if the external storage format is big-endian.</param>
        internal static void ReadArray<IntT>(Stream stream, Span<IntT> array, bool bigendp = false)
            where IntT : struct
        {
            var buffer = MemoryMarshal.Cast<IntT, byte>(array);

            var c = stream.Read(buffer);
            if (c != buffer.Length)
                throw new GeographicException("Failure reading data");

            if (bigendp == BitConverter.IsLittleEndian)
            {
                Swab(array);
            }
        }

        internal static bool ParseLine(
            ReadOnlySpan<char> line, out string key, out string value, 
            char equals = '\0', char comment = '#')
        {
            key = value = null;

            var n = line.IndexOf(comment);
            var linea = n==-1?line: line.Slice(0, n).Trim();
            if (linea.IsEmpty) return false;

            n = equals != 0 ? linea.IndexOf(equals)
                : linea.IndexOfAny(" \t\n\v\f\r".AsSpan());
            key = linea.Slice(0, n).Trim().ToString();
            if (string.IsNullOrEmpty(key)) return false;

            if (n != -1) value = linea.Slice(n + 1).Trim().ToString();

            return true;
        }

        internal static int FindFirstNotOf(this string source, string chars, int offset = 0) => source.AsSpan().FindFirstNotOf(chars, offset);

        internal static int FindFirstNotOf(this ReadOnlySpan<char> source, string chars, int offset = 0)
        {
            if (source == null) throw new ArgumentNullException("source");
            if (chars == null) throw new ArgumentNullException("chars");

            if (source.Length == 0) return -1;
            if (chars.Length == 0) return 0;

            for (int i = offset; i < source.Length; i++)
            {
                if (chars.IndexOf(source[i]) == -1) return i;
            }

            return -1;
        }

        internal static int FindFirstOf(this ReadOnlySpan<char> source, string chars, int offset = 0)
        {
            if (source == null) throw new ArgumentNullException("source");
            if (chars == null) throw new ArgumentNullException("chars");

            if (source.Length == 0) return -1;
            if (chars.Length == 0) return 0;

            for (int i = offset; i < source.Length; i++)
            {
                if (chars.Contains(source[i])) return i;
            }

            return -1;
        }

        internal static int FindLastNotOf(this ReadOnlySpan<char> source, string chars, int offset = 0)
        {
            if (source == null) throw new ArgumentNullException("source");
            if (chars == null) throw new ArgumentNullException("chars");

            if (source.Length == 0) return -1;
            if (chars.Length == 0) return 0;

            for (int i = source.Length - 1; i >= offset; i--)
            {
                if (!chars.Contains(source[i])) return i;
            }

            return -1;
        }

        internal static int FindLastNotOf(this Span<char> source, string chars, int offset = 0) => FindLastNotOf((ReadOnlySpan<char>)source, chars, offset);

        /// <summary>
        /// Convert a string representing a date to a fractional year.
        /// </summary>
        /// <param name="s">the string to be converted.</param>
        /// <returns>the fractional year.</returns>
        /// <remarks>
        /// The string is first read as an ordinary number (e.g., 2010 or 2012.5);
        /// if this is successful, the value is returned.  Otherwise the string
        /// should be of the form yyyy-mm or yyyy-mm-dd and this is converted to a
        /// number with 2010-01-01 giving 2010.0 and 2012-07-03 giving 2012.5.  The
        /// string "now" is interpreted as the present date.
        /// </remarks>
        public static double FractionalYear(string s)
        {
            if (s.TryParseDouble(out var result))
            {
                return result;
            }

            var ymd = s == "now" ? DateTime.Now : s.ParseDateTime();
            return ymd.Year + Math.Round(ymd.DayOfYear / (double)(new DateTime(ymd.Year + 1, 1, 1) - new DateTime(ymd.Year, 1, 1)).Days, 1);
        }

        internal static DateTime ParseDateTime(this string input)
            => DateTime.Parse(input, CultureInfo.InvariantCulture);

        internal static double ParseDouble(this string input)
            => double.Parse(input, NumberStyles.Float | NumberStyles.AllowThousands, CultureInfo.InvariantCulture);

        internal static bool TryParseDouble(this string input, out double result)
            => double.TryParse(input, NumberStyles.Float | NumberStyles.AllowThousands, CultureInfo.InvariantCulture, out result);

        internal static int ParseInt32(this string input)
            => int.Parse(input, NumberStyles.Integer, CultureInfo.InvariantCulture);

        internal static bool TryParseInt32(this string input, out int result)
            => int.TryParse(input, NumberStyles.Integer, CultureInfo.InvariantCulture, out result);

        internal static uint ParseUInt32(this string input)
            => uint.Parse(input, NumberStyles.Integer, CultureInfo.InvariantCulture);

        internal static bool TryParseUInt32(this string input, out uint result)
            => uint.TryParse(input, NumberStyles.Integer, CultureInfo.InvariantCulture, out result);

        internal static long ParseInt64(this string input)
            => long.Parse(input, NumberStyles.Integer, CultureInfo.InvariantCulture);

        internal static bool TryParseInt64(this string input, out long result)
            => long.TryParse(input, NumberStyles.Integer, CultureInfo.InvariantCulture, out result);

#if NETSTANDARD2_0
        internal static DateTime ParseDateTime(this ReadOnlySpan<char> input)
            => ParseDateTime(input.ToString());

        internal static double ParseDouble(this ReadOnlySpan<char> input)
            => ParseDouble(input.ToString());

        internal static bool TryParseDouble(this ReadOnlySpan<char> input, out double result)
            => TryParseDouble(input.ToString(), out result);

        internal static int ParseInt32(this ReadOnlySpan<char> input)
            => int.Parse(input.ToString(), NumberStyles.Integer, CultureInfo.InvariantCulture);

        internal static bool TryParseInt32(this ReadOnlySpan<char> input, out int result)
            => int.TryParse(input.ToString(), NumberStyles.Integer, CultureInfo.InvariantCulture, out result);

        internal static uint ParseUInt32(this ReadOnlySpan<char> input)
            => uint.Parse(input.ToString(), NumberStyles.Integer, CultureInfo.InvariantCulture);

        internal static bool TryParseUInt32(this ReadOnlySpan<char> input, out uint result)
            => uint.TryParse(input.ToString(), NumberStyles.Integer, CultureInfo.InvariantCulture, out result);

        internal static long ParseInt64(this ReadOnlySpan<char> input)
            => ParseInt64(input.ToString());

        internal static bool TryParseInt64(this ReadOnlySpan<char> input, out long result)
            => TryParseInt64(input.ToString(), out result);
#else
        internal static DateTime ParseDateTime(this ReadOnlySpan<char> input)
            => DateTime.Parse(input, CultureInfo.InvariantCulture);

        internal static double ParseDouble(this ReadOnlySpan<char> input)
            => double.Parse(input, NumberStyles.Float | NumberStyles.AllowThousands, CultureInfo.InvariantCulture);

        internal static bool TryParseDouble(this ReadOnlySpan<char> input, out double result)
            => double.TryParse(input, NumberStyles.Float | NumberStyles.AllowThousands, CultureInfo.InvariantCulture, out result);

        internal static int ParseInt32(this ReadOnlySpan<char> input)
            => int.Parse(input, NumberStyles.Integer, CultureInfo.InvariantCulture);

        internal static bool TryParseInt32(this ReadOnlySpan<char> input, out int result)
            => int.TryParse(input, NumberStyles.Integer, CultureInfo.InvariantCulture, out result);

        internal static uint ParseUInt32(this ReadOnlySpan<char> input)
            => uint.Parse(input, NumberStyles.Integer, CultureInfo.InvariantCulture);

        internal static bool TryParseUInt32(this ReadOnlySpan<char> input, out uint result)
            => uint.TryParse(input, NumberStyles.Integer, CultureInfo.InvariantCulture, out result);

        internal static long ParseInt64(this ReadOnlySpan<char> input)
            => long.Parse(input, NumberStyles.Integer, CultureInfo.InvariantCulture);

        internal static bool TryParseInt64(this ReadOnlySpan<char> input, out long result)
            => long.TryParse(input, NumberStyles.Integer, CultureInfo.InvariantCulture, out result);
#endif
    }
}