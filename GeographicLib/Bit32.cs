using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace GeographicLib
{
    /// <summary>
    /// Provides a uniform way to reinterpret 32-bit binary data as <see cref="System.Single"/>, <see cref="System.Int32"/> or <see cref="System.UInt32"/>.
    /// </summary>
    [StructLayout(LayoutKind.Explicit)]
    internal struct Bit32
    {
        [FieldOffset(0)]
        public float Single;

        [FieldOffset(0)]
        public int Int32;

        [FieldOffset(0)]
        public uint UInt32;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static implicit operator Bit32(float v) => new Bit32 { Single = v };

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static implicit operator Bit32(int v) => new Bit32 { Int32 = v };

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static implicit operator Bit32(uint v) => new Bit32 { UInt32 = v };

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static explicit operator float(Bit32 v) => v.Single;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static explicit operator int(Bit32 v) => v.Int32;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static explicit operator uint(Bit32 v) => v.UInt32;
    }
}