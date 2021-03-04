using GeographicLib.Geocodes;
using System;
using System.Collections.Generic;
using System.Text;

namespace GeographicLib
{
    /// <summary>
    /// Provide conversion between geographic coordinates.
    /// </summary>
    /// <remarks>
    /// This class stores a geographic position which may be set via the constructors or <see cref="Reset"/> via
    /// <list type="bullet">
    /// <item>latitude and longitude</item>
    /// <item>UTM or UPS coordinates</item>
    /// <item>a string representation of these or an <see cref="MGRS"/> coordinate string</item>
    /// </list>
    /// The state consists of the latitude and longitude and the supplied UTM or UPS coordinates (possibly derived from the <see cref="MGRS"/> coordinates).
    /// If latitude and longitude were given then the UTM/UPS coordinates follows the standard conventions.
    /// <para>
    /// The mutable state consists of the UTM or UPS coordinates for a alternate zone.
    /// A property <see cref="AltZone"/> is provided to set the alternate UPS/UTM zone.
    /// </para>
    /// <para>
    /// Methods are provided to return the geographic coordinates, the input UTM or UPS coordinates (and associated meridian convergence and scale),
    /// or alternate UTM or UPS coordinates (and their associated meridian convergence and scale).
    /// </para>
    /// <para>
    /// Once the input string has been parsed, you can print the result out in any of the formats, decimal degrees, degrees minutes seconds,
    /// <see cref="MGRS"/>, UTM/UPS.
    /// </para>
    /// </remarks>
    public class GeoCoords
    {
    }
}
