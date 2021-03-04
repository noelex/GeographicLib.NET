using System;
using System.Collections.Generic;
using System.Text;

namespace GeographicLib
{
    /// <summary>
    /// In this class we bring together the UTM and UPS coordinates systems.
    /// The UTM divides the earth between latitudes −80° and 84° into 60 zones numbered 1 thru 60.
    /// Zone assign zone number 0 to the UPS regions, covering the two poles.
    /// Within <see cref="UTMUPS"/>, non-negative zone numbers refer to one of the "physical" zones,
    /// 0 for UPS and [1, 60] for UTM. Negative "pseudo-zone" numbers are used to select one of the physical zones.
    /// </summary>
    public enum ZoneSpec
    {
        /// <summary>
        /// The smallest pseudo-zone number.
        /// </summary>
        MinPseudoZone = -4,

        /// <summary>
        /// A marker for an undefined or invalid zone. Equivalent to NaN.
        /// </summary>
        Invalid = -4,

        /// <summary>
        /// If a coordinate already include zone information (e.g., it is an <see cref="Geocodes.MGRS"/> coordinate),
        /// use that, otherwise apply the <see cref="Standard"/> rules.
        /// </summary>
        Match = -3,

        /// <summary>
        /// Apply the standard rules for UTM zone assigment extending the UTM zone to each pole to give a zone number in [1, 60].
        /// For example, use UTM zone 38 for longitude in [42°, 48°). The rules include the Norway and Svalbard exceptions.
        /// </summary>
        UTM = -2,

        /// <summary>
        /// Apply the standard rules for zone assignment to give a zone number in [0, 60].
        /// If the latitude is not in [−80°, 84°), then use <see cref="UPS"/> = 0,
        /// otherwise apply the rules for <see cref="UTM"/>.
        /// The tests on latitudes and longitudes are all closed on the lower end open on the upper.
        /// Thus for UTM zone 38, latitude is in [−80°, 84°) and longitude is in [42°, 48°).
        /// </summary>
        Standard = -1,

        /// <summary>
        /// The largest pseudo-zone number.
        /// </summary>
        MaxPseudoZone = -1,

        /// <summary>
        /// The smallest physical zone number.
        /// </summary>
        MinZone = 0,

        /// <summary>
        /// The zone number used for UPS.
        /// </summary>
        UPS = 0,

        /// <summary>
        /// The smallest UTM zone number.
        /// </summary>
        MinUTMZone = 1,

        /// <summary>
        /// The largest UTM zone number.
        /// </summary>
        MaxUTMZone = 60,

        /// <summary>
        /// The largest physical zone number.
        /// </summary>
        MaxZone = 60
    }
}
