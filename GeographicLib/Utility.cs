using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.InteropServices;
using System.Text;

namespace GeographicLib
{
    internal static class Utility
    {
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

        public static bool ParseLine(ReadOnlySpan<char> line, out string key, out string value, char delim= '\0')
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

        public static int FindFirstNotOf(this ReadOnlySpan<char> source, string chars)
        {
            if (source == null) throw new ArgumentNullException("source");
            if (chars == null) throw new ArgumentNullException("chars");

            if (source.Length == 0) return -1;
            if (chars.Length == 0) return 0;

            for (int i = 0; i < source.Length; i++)
            {
                if (chars.IndexOf(source[i]) == -1) return i;
            }

            return -1;
        }
    }
}