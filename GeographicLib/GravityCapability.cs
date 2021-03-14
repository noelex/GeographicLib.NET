using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib
{
    [Flags]
    enum GravityCapability:uint
    {
        None = 0U,
        G = 1U << 0,       // implies potentials W and V
        T = 1U << 1,
        Delta = 1U << 2 | T, // delta implies T?
        C = 1U << 3,
        Gamma0 = 1U << 4,
        Gamma = 1U << 5,
        All = 0x3FU,
    }
}
