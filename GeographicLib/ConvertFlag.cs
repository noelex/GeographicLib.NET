using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib
{
    /// <summary>
    /// Flags indicating conversions between heights above the geoid and heights above the ellipsoid.
    /// </summary>
    public enum ConvertFlag
    {
        /// <summary>
        /// The multiplier for converting from heights above the ellipsoid to heights above the geoid.
        /// </summary>
        EllipsoidToGeoid = -1,

        /// <summary>
        /// No conversion.
        /// </summary>
        None = 0,

        /// <summary>
        /// The multiplier for converting from heights above the geoid to heights above the ellipsoid.
        /// </summary>
        GeoidToEllipsoid = 1
    }
}
