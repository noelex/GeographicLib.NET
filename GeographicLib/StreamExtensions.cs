#if NETSTANDARD2_0
using System;
using System.Buffers;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib
{
    static class StreamExtensions
    {
        public static int Read(this Stream source, Span<byte> buffer)
        {
            var arr = ArrayPool<byte>.Shared.Rent(buffer.Length);
            try
            {
                var count=source.Read(arr, 0, buffer.Length);
                arr.AsSpan().Slice(0, count).CopyTo(buffer);
                return count;
            }
            finally
            {
                ArrayPool<byte>.Shared.Return(arr);
            }
        }
    }
}
#endif