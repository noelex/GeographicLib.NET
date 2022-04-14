using System;
using System.Collections.Generic;
using System.Text;

namespace GeographicLib
{
    /// <summary>
    /// Bit masks for what calculations to do.  These masks do double duty.
    /// They signify to the <see cref="GeodesicLine"/> constructor and to
    /// <see cref="Geodesic.Line"/> what capabilities should be included in the <see cref="GeodesicLine"/>
    /// object.  They also specify which results to return in the general
    /// routines <see cref="Geodesic.GenDirect"/> and <see cref="Geodesic.GenInverse(double, double, double, double, GeodesicFlags, out double, out double, out double, out double, out double, out double, out double)"/> routines.
    /// </summary>
    [Flags]
    public enum GeodesicFlags : uint
    {
        /// <summary>
        /// No capabilities, no output.
        /// </summary>
        None = 0,

        /// <summary>
        /// Calculate latitude <i>lat2</i>.  (It's not necessary to include this as a
        /// capability to <see cref="GeodesicLine"/> or <see cref="GeodesicLineExact"/> because this is included by default.)
        /// </summary>
        Latitude = 1 << 7 | GeodesicCapability.None,

        /// <summary>
        /// Calculate latitude <i>lon2</i>.
        /// </summary>
        Longitude = 1 << 8 | GeodesicCapability.C3,

        /// <summary>
        /// Calculate azimuths <i>azi1</i> and <i>azi2</i>.
        /// (It's not necessary to include this as a capability to <see cref="GeodesicLine"/> or <see cref="GeodesicLineExact"/>
        /// because this is included by default.)
        /// </summary>
        Azimuth = 1 << 9 | GeodesicCapability.None,

        /// <summary>
        /// Calculate distance <i>s12</i>.
        /// </summary>
        Distance = 1 << 10 | GeodesicCapability.C1,

        /// <summary>
        /// A combination of the common capabilities:
        /// <see cref="Latitude"/>, <see cref="Longitude"/>,
        /// <see cref="Azimuth"/> and <see cref="Distance"/>.
        /// </summary>
        Standard = Latitude | Longitude | Azimuth | Distance,

        /// <summary>
        /// Allow distance <i>s12</i> to be used as input in the direct geodesic problem.
        /// </summary>
        DistanceIn = 1 << 11 | GeodesicCapability.C1 | GeodesicCapability.C1p,

        /// <summary>
        /// Calculate reduced length <i>m12</i>.
        /// </summary>
        ReducedLength = 1 << 12 | GeodesicCapability.C1 | GeodesicCapability.C2,

        /// <summary>
        /// Calculate geodesic scales <i>M12</i> and <i>M21</i>.
        /// </summary>
        GeodesicScale = 1 << 13 | GeodesicCapability.C1 | GeodesicCapability.C2,

        /// <summary>
        /// Calculate area <i>S12</i>.
        /// </summary>
        Area = 1 << 14 | GeodesicCapability.C4,

        /// <summary>
        /// Unroll <i>lon2</i> in the direct calculation.
        /// </summary>
        LongUnroll = 1 << 15,

        /// <summary>
        /// All capabilities, calculate everything.  (<see cref="LongUnroll"/> is not included in this mask.)
        /// </summary>
        All = GeodesicCapability.OutAll | GeodesicCapability.All
    }

    static class GeodesicFlagsExtesions
    {
        public static GeodesicFlags Flags(this GeodesicFlags src) => src & (GeodesicFlags)GeodesicCapability.OutMask;

        public static GeodesicCapability Capabilities(this GeodesicFlags src) => (GeodesicCapability)src & GeodesicCapability.All;

        /// <summary>
        /// Verifies whether <paramref name="src"/> and <paramref name="flags"/> has any common flag. 
        /// </summary>
        /// <param name="src"></param>
        /// <param name="flags"></param>
        /// <returns></returns>
        public static bool HasAny(this GeodesicFlags src, GeodesicFlags flags) =>  (uint)(src & flags) != 0;

        /// <summary>
        /// Verifies whether <paramref name="src"/> and <paramref name="caps"/> has any common flag. 
        /// </summary>
        /// <param name="src"></param>
        /// <param name="caps"></param>
        /// <returns></returns>
        public static bool HasAny(this GeodesicFlags src, GeodesicCapability caps) => (uint)(src & (GeodesicFlags)caps) != 0;

        // public static bool HasCapabilities(this GeodesicFlags src, Capability flags) => (uint)(src & (GeodesicFlags)flags) != 0;
    }
}
