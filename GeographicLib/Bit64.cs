using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace GeographicLib
{
    /// <summary>
    /// Provides a uniform way to reinterpret 64-bit binary data as <see cref="System.Double"/>, <see cref="System.Int64"/> or <see cref="System.UInt64"/>.
    /// </summary>
    [StructLayout(LayoutKind.Explicit)]
    internal struct Bit64
    {
        [FieldOffset(0)]
        public double Double;

        [FieldOffset(0)]
        public long Int64;

        [FieldOffset(0)]
        public ulong UInt64;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static implicit operator Bit64(double v) => new Bit64 { Double = v };

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static implicit operator Bit64(long v) => new Bit64 { Int64 = v };

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static implicit operator Bit64(ulong v) => new Bit64 { UInt64 = v };

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static explicit operator double(Bit64 v) => v.Double;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static explicit operator long(Bit64 v) => v.Int64;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static explicit operator ulong(Bit64 v) => v.UInt64;
    }
}

