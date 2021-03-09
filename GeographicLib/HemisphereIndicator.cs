using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib
{
    /// <summary>
    /// Indicator for presence of hemisphere indicator (N/S/E/W) on latitudes and longitudes.
    /// </summary>
    public enum HemisphereIndicator
    {
        /// <summary>
        /// No indicator present.
        /// </summary>
        None = 0,

        /// <summary>
        /// Latitude indicator (N/S) present.
        /// </summary>
        Latitude = 1,

        /// <summary>
        /// Longitude indicator (E/W) present.
        /// </summary>
        Longitude = 2,

        /// <summary>
        /// Used in <see cref="DMS.Encode(double, int, HemisphereIndicator, char)"/> to indicate output of an azimuth in [000, 360) with no letter indicator.
        /// </summary>
        Azimuth = 3,

        /// <summary>
        /// Used in <see cref="DMS.Encode(double, int, HemisphereIndicator, char)"/> to indicate output of a plain number.
        /// </summary>
        Number = 4
    }
}
