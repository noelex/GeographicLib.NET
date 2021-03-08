using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib
{
    /// <summary>
    /// Indicator for trailing units on an angle.
    /// </summary>
    public enum TrailingUnit
    {
        /// <summary>
        /// Trailing unit is degrees.
        /// </summary>
        Degree = 0,

        /// <summary>
        /// Trailing unit is arc minutes.
        /// </summary>
        Minute = 1,

        /// <summary>
        /// Trailing unit is arc seconds.
        /// </summary>
        Second = 2
    }
}
