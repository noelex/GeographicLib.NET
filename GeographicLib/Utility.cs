using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Text;

namespace GeographicLib
{
    internal static class Utility
    {
        private readonly static FieldInfo charPosField = typeof(StreamReader).GetField("_charPos", BindingFlags.NonPublic | BindingFlags.Instance | BindingFlags.DeclaredOnly);
        private readonly static FieldInfo charLenField = typeof(StreamReader).GetField("_charLen", BindingFlags.NonPublic | BindingFlags.Instance | BindingFlags.DeclaredOnly);
        private readonly static FieldInfo charBufferField = typeof(StreamReader).GetField("_byteBuffer", BindingFlags.NonPublic | BindingFlags.Instance | BindingFlags.DeclaredOnly);

        public static long Position(this StreamReader reader)
        {
            var byteBuffer = (byte[])charBufferField.GetValue(reader);
            var charLen = (int)charLenField.GetValue(reader);
            var charPos = (int)charPosField.GetValue(reader);

            return reader.BaseStream.Position - byteBuffer.Length - charPos;
            // reader.CurrentEncoding.GetByteCount(charBuffer, charPos, charLen - charPos);
        }

        public static bool IsInteger<T>()
        {
            return typeof(T) == typeof(sbyte) ||
                typeof(T) == typeof(byte) ||
                typeof(T) == typeof(short) ||
                typeof(T) == typeof(ushort) ||
                typeof(T) == typeof(int) ||
                typeof(T) == typeof(uint) ||
                typeof(T) == typeof(long) ||
                typeof(T) == typeof(ulong);
        }

#if NETSTANDARD2_0
        public static int IndexOf(this string str, char c, StringComparison stringComparison)
            => str.IndexOf(c.ToString(), stringComparison);
#endif
        /// <summary>
        /// Swap byte order of elements in the specified array.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="array"></param>
        public static void Swab<T>(Span<T> array) where T : struct
        {
            for (int i = 0; i < array.Length; i++)
            {
                MemoryMarshal.Cast<T, byte>(array.Slice(i, 1)).Reverse();
            }
        }

        public static string ToFixedString(this double x, int p = -1)
        {
            if (!MathEx.IsFinite(x))
                return x < 0 ? "-inf" :
                    (x > 0 ? "inf" : "nan");
            return p >= 0 ? x.ToString($"F{p}") : x.ToString();
        }

        public static double ParseDouble(this string s) => ParseDouble(s.AsSpan());

        public static double ParseDouble(this ReadOnlySpan<char> s)
        {
            double x;

            var t = s.Trim();
            var sep = t.IndexOfAny(" \t\n\v\f\r".AsSpan());
            if (sep > 0)
            {
                throw new GeographicException($"Extra text {t.Slice(sep + 1).ToString()} at end of {t.ToString()}");
            }


            if (!double.TryParse(t.ToString(), out x) && (x = t.NumMatch()) != 0)
            {
                throw new GeographicException($"Cannot decode {t.ToString()}");
            }

            return x;
        }

        public static double NumMatch(this string s) => s.AsSpan().NumMatch();

        public static double NumMatch(this ReadOnlySpan<char> s)
        {
            if (s.Length < 3) return 0;

            Span<char> t = stackalloc char[s.Length];
            s.ToUpper(t, System.Globalization.CultureInfo.CurrentCulture);

            var sign = t[0] == '-' ? -1 : 1;
            var p0 = t[0] == '-' || t[0] == '+' ? 1 : 0;
            var p1 = t.FindLastNotOf("0");

            if (p1 == -1 || p1 + 1 < p0 + 3)
                return 0;

            if (s.SequenceEqual("NAN".AsSpan()) || s.SequenceEqual("1.#QNAN".AsSpan()) ||
               s.SequenceEqual("1.#SNAN".AsSpan()) || s.SequenceEqual("1.#IND".AsSpan()) || s.SequenceEqual("1.#R".AsSpan()))
            {
                return double.NaN;
            }

            if (t.SequenceEqual("INF".AsSpan()) || t.SequenceEqual("1.#INF".AsSpan()) || t.SequenceEqual("∞".AsSpan()))
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
        public static void ReadArray<IntT>(Stream stream, Span<IntT> array, bool bigendp = false)
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

        public static bool ParseLine(ReadOnlySpan<char> line, out string key, out string value, char delim = '\0')
        {
            key = value = null;

            var n = line.IndexOf('#');
            var linea = line.Slice(0, n).Trim();
            if (linea.IsEmpty) return false;

            n = delim != 0 ? linea.IndexOf(delim) : linea.IndexOfAny(" \t\n\v\f\r".AsSpan());
            key = linea.Slice(0, n).Trim().ToString();
            if (string.IsNullOrEmpty(key)) return false;

            if (n != -1) value = linea.Slice(n + 1).ToString();

            return true;
        }

        public static int FindFirstNotOf(this string source, string chars, int offset = 0) => source.AsSpan().FindFirstNotOf(chars, offset);

        public static int FindFirstNotOf(this ReadOnlySpan<char> source, string chars, int offset = 0)
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

        public static int FindLastNotOf(this ReadOnlySpan<char> source, string chars, int offset = 0)
        {
            if (source == null) throw new ArgumentNullException("source");
            if (chars == null) throw new ArgumentNullException("chars");

            if (source.Length == 0) return -1;
            if (chars.Length == 0) return 0;

            for (int i = source.Length - 1; i >= offset; i--)
            {
                if (chars.IndexOf(source[i]) == -1) return i;
            }

            return -1;
        }

        public static int FindLastNotOf(this Span<char> source, string chars, int offset = 0) => FindLastNotOf((ReadOnlySpan<char>)source, chars, offset);
    }
}