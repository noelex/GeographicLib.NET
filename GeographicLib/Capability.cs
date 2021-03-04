using System;
using System.Collections.Generic;
using System.Text;

namespace GeographicLib
{
    [Flags]
    internal enum Capability : uint
    {
        None = 0,
        C1 = 1 << 0,
        C1p = 1 << 1,
        C2 = 1U << 2,
        C3 = 1U << 3,
        C4 = 1U << 4,
        All = 0x1FU,
        Mask = All,
        OutAll = 0x7F80U,
        OutMask = 0xFF80U,       // Includes LONG_UNROLL,

        E = C1,
        D = C2,
        H = C3
    }
}
