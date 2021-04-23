using GeographicLib.Geocodes;
using System;
using System.Collections.Generic;
using System.Text;

using static System.Math;
using static GeographicLib.MathEx;
using GeographicLib.Projections;

namespace GeographicLib
{
    /// <summary>
    /// Convert between geographic coordinates and UTM/UPS.
    /// </summary>
    /// <remarks>
    /// UTM and UPS are defined
    /// <list type="bullet">
    /// <item>
    /// J. W. Hager, J. F. Behensky, and B. W. Drew,
    /// <a href="https://web.archive.org/web/20161214054445/http://earth-info.nga.mil/GandG/publications/tm8358.2/TM8358_2.pdf">The Universal Grids: Universal Transverse Mercator (UTM) and Universal Polar Stereographic (UPS)</a>,
    /// Defense Mapping Agency, Technical Manual TM8358.2 (1989).
    /// </item>
    /// </list>
    /// Section 2-3 defines UTM and section 3-2.4 defines UPS.
    /// This document also includes approximate algorithms for the computation of the underlying transverse Mercator and polar stereographic projections.
    /// Here we substitute much more accurate algorithms given by <see cref="Projections.TransverseMercator"/> and <see cref="Projections.PolarStereographic"/>.
    /// These are the algorithms recommended by the NGA document
    /// <list type="bullet">
    /// <item>
    /// <a href="https://earth-info.nga.mil/coordsys/coord-download.php?file=website/NGA.SIG.0012_2.0.0_UTMUPS.pdf">The Universal Grids and the Transverse Mercator and Polar Stereographic Map Projections</a>,
    /// NGA.SIG.0012_2.0.0_UTMUPS (2014).
    /// </item>
    /// </list>
    /// In this implementation, the conversions are closed, i.e., output from Forward is legal input for Reverse and vice versa.
    /// The error is about 5nm in each direction.
    /// However, the conversion from legal UTM/UPS coordinates to geographic coordinates and back might throw an error if the initial point
    /// is within 5nm of the edge of the allowed range for the UTM/UPS coordinates.
    /// <para>
    /// The simplest way to guarantee the closed property is to define allowed ranges for the eastings and northings for UTM and UPS coordinates.
    /// The UTM boundaries are the same for all zones.
    /// (The only place the exceptional nature of the zone boundaries is evident is when converting to UTM/UPS coordinates requesting the standard zone.)
    /// The MGRS lettering scheme imposes natural limits on UTM/UPS coordinates which may be converted into MGRS coordinates.
    /// For the conversion to/from geographic coordinates these ranges have been extended by 100km in order to provide a generous overlap
    /// between UTM and UPS and between UTM zones.
    /// </para>
    /// <para>
    /// The <a href="http://www.nga.mil/">NGA</a> software package <a href="https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84#tab_geotrans">geotrans</a> also provides conversions to and from UTM and UPS.
    /// Version 2.4.2 (and earlier) suffers from some drawbacks:
    /// </para>
    /// <list type="bullet">
    /// <item>Inconsistent rules are used to determine the whether a particular UTM or UPS coordinate is legal. A more systematic approach is taken here.</item>
    /// <item>The underlying projections are not very accurately implemented.</item>
    /// </list>
    /// The <see cref="EncodeZone"/> encodes the UTM zone and hemisphere to allow UTM/UPS coordinated to be displayed as,
    /// for example, "<c>38N 444500 3688500</c>". According to NGA.SIG.0012_2.0.0_UTMUPS the use of "<c>N</c>" to denote "north" in the context is not 
    /// allowed (since a upper case letter in this context denotes the MGRS latitude band).
    /// Consequently, as of version 1.36, <see cref="EncodeZone"/> uses the lower case letters "<c>n</c>" and "<c>s</c>" to denote the hemisphere.
    /// In addition <see cref="EncodeZone"/> accepts an optional final argument <i>abbrev</i>, which, if <see langword="false"/>,
    /// results in the hemisphere being spelled out as in "<c>38north</c>".
    /// </remarks>
    public static class UTMUPS
    {
        private static readonly ReadOnlyMemory<int>
            falseeasting_ = new[] { MGRS.upseasting_ * MGRS.tile_, MGRS.upseasting_ * MGRS.tile_,
                                    MGRS.utmeasting_ * MGRS.tile_, MGRS.utmeasting_ * MGRS.tile_ },
            falsenorthing_ = new[] { MGRS.upseasting_ * MGRS.tile_, MGRS.upseasting_ * MGRS.tile_,
                                     MGRS.maxutmSrow_ * MGRS.tile_, MGRS.minutmNrow_ * MGRS.tile_ },
            mineasting_ = new[] { MGRS.minupsSind_ * MGRS.tile_, MGRS.minupsNind_ * MGRS.tile_,
                                  MGRS.minutmcol_ * MGRS.tile_, MGRS.minutmcol_ * MGRS.tile_ },
            maxeasting_ = new[] { MGRS.maxupsSind_ * MGRS.tile_, MGRS.maxupsNind_ * MGRS.tile_,
                                  MGRS.maxutmcol_ * MGRS.tile_, MGRS.maxutmcol_ * MGRS.tile_ },
            minnorthing_ = new[] { MGRS.minupsSind_ * MGRS.tile_, MGRS.minupsNind_ * MGRS.tile_,
                                   MGRS.minutmSrow_ * MGRS.tile_, (MGRS.minutmNrow_ + MGRS.minutmSrow_ - MGRS.maxutmSrow_) * MGRS.tile_ },
            maxnorthing_ = new[] { MGRS.maxupsSind_ * MGRS.tile_, MGRS.maxupsNind_ * MGRS.tile_,
                                  (MGRS.maxutmSrow_ + MGRS.maxutmNrow_ - MGRS.minutmNrow_) * MGRS.tile_, MGRS.maxutmNrow_ * MGRS.tile_ };

        private const int epsg01N = 32601, // EPSG code for UTM 01N
                          epsg60N = 32660, // EPSG code for UTM 60N
                          epsgN = 32661, // EPSG code for UPS   N
                          epsg01S = 32701, // EPSG code for UTM 01S
                          epsg60S = 32760, // EPSG code for UTM 60S
                          epsgS = 32761; // EPSG code for UPS   S

        private static double CentralMeridian(int zone) => 6.0 * zone - 183;

        private static bool CheckCoords(bool utmp, bool northp, double x, double y,
                            bool mgrslimits = false, bool throwp = true)
        {
            // Limits are all multiples of 100km and are all closed on the both ends.
            // Failure tests are such that NaNs succeed.
            var slop = mgrslimits ? 0 : MGRS.tile_;
            int ind = (utmp ? 2 : 0) + (northp ? 1 : 0);
            if (x < mineasting_.Span[ind] - slop || x > maxeasting_.Span[ind] + slop)
            {
                if (!throwp) return false;
                throw new GeographicException($"Easting {x / 1000}km not in {(mgrslimits ? "MGRS/" : "")}{(utmp ? "UTM" : "UPS")} range for " +
                    $"{(northp ? "N" : "S")}  hemisphere [{(mineasting_.Span[ind] - slop) / 1000}km, {(maxeasting_.Span[ind] + slop) / 1000}km]");
            }
            if (y < minnorthing_.Span[ind] - slop || y > maxnorthing_.Span[ind] + slop)
            {
                if (!throwp) return false;
                throw new GeographicException($"Northing {y / 1000}km not in {(mgrslimits ? "MGRS/" : "")}{(utmp ? "UTM" : "UPS")} range for " +
                    $"{(northp ? "N" : "S")} hemisphere [{(minnorthing_.Span[ind] - slop) / 1000}km, {(maxnorthing_.Span[ind] + slop) / 1000}km]");
            }
            return true;
        }

        /// <summary>
        /// Gets a value representing the equatorial radius of the WGS84 ellipsoid (meters).
        /// </summary>
        /// <remarks>
        /// The WGS84 value is returned because the UTM and UPS projections are based on this ellipsoid.
        /// </remarks>
        public static double EquatorialRadius => Constants.WGS84_a;

        /// <summary>
        /// Gets a value representing the flattening of the WGS84 ellipsoid.
        /// </summary>
        /// <remarks>
        /// The WGS84 value is returned because the UTM and UPS projections are based on this ellipsoid.
        /// </remarks>
        public static double Flattening => Constants.WGS84_f;

        /// <summary>
        /// Gets a value representing the shift (meters) necessary to align north and south halves of a UTM zone (10^7).
        /// </summary>
        public static double UTMShift => MGRS.utmNshift_;

        /// <summary>
        /// Gets the standard zone.
        /// </summary>
        /// <param name="lat">latitude (degrees).</param>
        /// <param name="lon">longitude (degrees).</param>
        /// <param name="setzone">
        /// zone override (optional). If omitted, use the standard rules for picking the zone.
        /// If setzone is given then use that zone if it is non-negative, otherwise apply the rules given in <see cref="ZoneSpec"/>.
        /// </param>
        /// <returns></returns>
        /// <remarks>
        /// This is exact.
        /// </remarks>
        public static int StandardZone(double lat, double lon, int setzone = (int)ZoneSpec.Standard)
        {
            if (!(setzone >= (int)ZoneSpec.MinPseudoZone && setzone <= (int)ZoneSpec.MaxZone))
                throw new GeographicException($"Illegal zone requested {setzone}");
            if (setzone >= (int)ZoneSpec.MinZone || setzone == (int)ZoneSpec.Invalid)
                return setzone;
            if (double.IsNaN(lat) || double.IsNaN(lon)) // Check if lat or lon is a NaN
                return (int)ZoneSpec.Invalid;
            if (setzone == (int)ZoneSpec.UTM || (lat >= -80 && lat < 84))
            {
                int ilon = (int)Floor(AngNormalize(lon));
                if (ilon == 180) ilon = -180; // ilon now in [-180,180)
                int zone = (ilon + 186) / 6;
                int band = MGRS.LatitudeBand(lat);
                if (band == 7 && zone == 31 && ilon >= 3) // The Norway exception
                    zone = 32;
                else if (band == 9 && ilon >= 0 && ilon < 42) // The Svalbard exception
                    zone = 2 * ((ilon + 183) / 12) + 1;
                return zone;
            }
            else
                return (int)ZoneSpec.UPS;
        }

        /// <summary>
        /// <see cref="Reverse(int, bool, double, double, out double, out double, bool)"/> without returning convergence and scale.
        /// </summary>
        /// <param name="zone">the UTM zone (zero means UPS).</param>
        /// <param name="northp">hemisphere (true means north, false means south).</param>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <param name="mgrslimits">if <see langword="true"/> enforce the stricter MGRS limits on the coordinates(default = <see langword="false"/>).</param>
        /// <returns>
        /// <i>lat</i>, latitude of point (degrees) and <i>lon</i>, longitude of point (degrees).
        /// </returns>
        public static (double lat, double lon) Reverse(int zone, bool northp, double x, double y, bool mgrslimits = false)
            => Reverse(zone, northp, x, y, out _, out _, mgrslimits);

        /// <summary>
        /// Reverse projection, from UTM/UPS to geographic.
        /// </summary>
        /// <param name="zone">the UTM zone (zero means UPS).</param>
        /// <param name="northp">hemisphere (true means north, false means south).</param>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <param name="gamma">meridian convergence at point (degrees).</param>
        /// <param name="k">scale of projection at point.</param>
        /// <param name="mgrslimits">if <see langword="true"/> enforce the stricter MGRS limits on the coordinates(default = <see langword="false"/>).</param>
        /// <returns>
        /// <i>lat</i>, latitude of point (degrees) and <i>lon</i>, longitude of point (degrees).
        /// </returns>
        /// <remarks>
        /// The accuracy of the conversion is about 5nm.
        /// <para>
        /// UTM eastings are allowed to be in the range [0km, 1000km], northings are allowed to be in in [0km, 9600km] for the northern hemisphere
        /// and in [900km, 10000km] for the southern hemisphere.
        /// However UTM northings can be continued across the equator.
        /// So the actual limits on the northings are [-9100km, 9600km] for the "northern" hemisphere and [900km, 19600km] for the "southern" hemisphere.
        /// </para>
        /// <para>
        /// UPS eastings and northings are allowed to be in the range [1200km, 2800km] in the northern hemisphere and in [700km, 3300km] in the southern hemisphere.
        /// </para>
        /// <para>
        /// These ranges are 100km larger than allowed for the conversions to <see cref="MGRS"/>.
        /// (100km is the maximum extra padding consistent with eastings remaining non-negative.)
        /// This allows generous overlaps between zones and UTM and UPS. If <paramref name="mgrslimits"/> = <see langword="true"/>,
        /// then all the ranges are shrunk by 100km so that they agree with the stricter <see cref="MGRS"/> ranges.
        /// No checks are performed besides these (e.g., to limit the distance outside the standard zone boundaries).
        /// </para>
        /// </remarks>
        public static (double lat, double lon) Reverse(
            int zone, bool northp, double x, double y, out double gamma, out double k, bool mgrslimits = false)
        {
            gamma = k = double.NaN;

            if (zone == (int)ZoneSpec.Invalid || double.IsNaN(x) || double.IsNaN(y))
            {
                return (double.NaN, double.NaN);
            }

            if (!(zone >= (int)ZoneSpec.MinZone && zone <= (int)ZoneSpec.MaxZone))
                throw new GeographicException($"Zone {zone} not in range [0, 60]");
            bool utmp = zone != (int)ZoneSpec.UPS;
            CheckCoords(utmp, northp, x, y, mgrslimits);
            int ind = (utmp ? 2 : 0) + (northp ? 1 : 0);
            x -= falseeasting_.Span[ind];
            y -= falsenorthing_.Span[ind];
            if (utmp)
                return TransverseMercator.UTM.Reverse(CentralMeridian(zone), x, y, out gamma, out k);
            else
                return PolarStereographic.UPS.Reverse(northp, x, y, out gamma, out k);
        }

        /// <summary>
        /// <see cref="Forward(double, double, out double, out double, int, bool)"/> without convergence and scale.
        /// </summary>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <param name="setzone">zone override (optional).</param>
        /// <param name="mgrslimits">if <see langword="true"/> enforce the stricter MGRS limits on the coordinates(default = <see langword="false"/>).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>zone</i>, the UTM zone (zero means UPS).</item>
        /// <item><i>northp</i>, hemisphere (<see langword="true"/> means north, <see langword="false"/> means south).</item>
        /// <item><i>x</i>, easting of point (meters).</item>
        /// <item><i>y</i>, northing of point (meters).</item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// If <paramref name="setzone"/> is omitted, use the standard rules for picking the zone.
        /// If <paramref name="setzone"/> is given then use that zone if it is non-negative, otherwise apply the rules given in <see cref="ZoneSpec"/>.
        /// The accuracy of the conversion is about 5nm.
        /// <para>
        /// The northing <i>y</i> jumps by <see cref="UTMShift"/> when crossing the equator in the southerly direction.
        /// Sometimes it is useful to remove this discontinuity in <i>y</i> by extending the "northern" hemisphere using <see cref="Transfer"/>.
        /// </para>
        /// </remarks>
        public static (int zone, bool northp, double x, double y) Forward(
            double lat, double lon, int setzone = (int)ZoneSpec.Standard, bool mgrslimits = false)
            => Forward(lat, lon, out _, out _, setzone, mgrslimits);

        /// <summary>
        /// Forward projection, from geographic to UTM/UPS.
        /// </summary>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <param name="setzone">zone override (optional).</param>
        /// <param name="gamma">meridian convergence at point (degrees).</param>
        /// <param name="k">scale of projection at point.</param>
        /// <param name="mgrslimits">if <see langword="true"/> enforce the stricter MGRS limits on the coordinates(default = <see langword="false"/>).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>zone</i>, the UTM zone (zero means UPS).</item>
        /// <item><i>northp</i>, hemisphere (<see langword="true"/> means north, <see langword="false"/> means south).</item>
        /// <item><i>x</i>, easting of point (meters).</item>
        /// <item><i>y</i>, northing of point (meters).</item>
        /// </list>
        /// </returns>
        public static (int zone, bool northp, double x, double y) Forward(
            double lat, double lon, out double gamma, out double k, int setzone = (int)ZoneSpec.Standard, bool mgrslimits = false)
        {
            if (Abs(lat) > 90)
                throw new GeographicException($"Latitude {lat}d not in [-90d, 90d]");
            bool northp1 = lat >= 0;
            int zone1 = StandardZone(lat, lon, setzone);
            if (zone1 == (int)ZoneSpec.Invalid)
            {
                gamma = k = double.NaN;
                return (zone1, northp1, double.NaN, double.NaN);
            }

            double x1, y1, gamma1, k1;
            bool utmp = zone1 != (int)ZoneSpec.UPS;

            if (utmp)
            {
                double
                  lon0 = CentralMeridian(zone1),
                  dlon = lon - lon0;
                dlon = Abs(dlon - 360 * Floor((dlon + 180) / 360));
                if (!(dlon <= 60))
                    // Check isn't really necessary because CheckCoords catches this case.
                    // But this allows a more meaningful error message to be given.
                    throw new GeographicException($"Longitude {lon}d more than 60d from center of UTM zone {zone1}");
                (x1, y1) = TransverseMercator.UTM.Forward(lon0, lat, lon, out gamma1, out k1);
            }
            else
            {
                if (Abs(lat) < 70)
                    // Check isn't really necessary ... (see above).
                    throw new GeographicException($"Latitude {lat}d more than 20d from {(northp1 ? "N" : "S")} pole");
                (x1, y1) = PolarStereographic.UPS.Forward(northp1, lat, lon, out gamma1, out k1);
            }

            int ind = (utmp ? 2 : 0) + (northp1 ? 1 : 0);
            x1 += falseeasting_.Span[ind];
            y1 += falsenorthing_.Span[ind];

            if (!CheckCoords(zone1 != (int)ZoneSpec.UPS, northp1, x1, y1, mgrslimits, false))
                throw new GeographicException($"Latitude {lat}, longitude {lon} out of legal range for {(utmp ? ("UTM zone " + zone1) : "UPS")}");

            gamma = gamma1;
            k = k1;

            return (zone1, northp1, x1, y1);
        }

        /// <summary>
        /// Transfer UTM/UPS coordinated from one zone to another.
        /// </summary>
        /// <param name="zonein">the UTM zone for <paramref name="xin"/> and <paramref name="yin"/> (or zero for UPS).</param>
        /// <param name="northpin">hemisphere for <paramref name="xin"/> and <paramref name="yin"/>
        /// (<see langword="true"/> means north, <see langword="false"/> means south).</param>
        /// <param name="xin">easting of point (meters) in <paramref name="zonein"/>.</param>
        /// <param name="yin">northing of point (meters) in <paramref name="zonein"/>.</param>
        /// <param name="zoneout">the requested UTM zone for <i>xout</i> and <i>yout</i> (or zero for UPS).</param>
        /// <param name="northpout">hemisphere for <i>xout</i> and <i>yout</i>.</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>zone</i>, the actual UTM zone for <i>xout</i> and <i>yout</i> (or zero for UPS);
        /// this equals <paramref name="zoneout"/> if <paramref name="zoneout"/> ≥ 0.</item>
        /// <item><i>xout</i>, easting of point (meters) in <paramref>zoneout</paramref>.</item>
        /// <item><i>yout</i>, northing of point (meters) in <paramref>zoneout</paramref>.</item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// <paramref name="zonein"/> must be in the range [0, 60] with <paramref name="zonein"/> = <see cref="ZoneSpec.UPS"/>, 0, indicating UPS.
        /// <paramref name="zonein"/> may also be <see cref="ZoneSpec.Invalid"/>.
        /// <para>
        /// <paramref name="zoneout"/> must be in the range [-4, 60]. If <paramref name="zoneout"/> &lt; <see cref="ZoneSpec.MinZone"/> then the rules
        /// give in the documentation of <see cref="ZoneSpec"/> are applied, and <i>zone</i> is set to the actual zone used for output.
        /// </para>
        /// <para>
        /// (<i>xout, yout</i>) can overlap with (<i>xin, yin</i>).
        /// </para>
        /// </remarks>
        public static (int zone, double xout, double yout) Transfer(int zonein, bool northpin, double xin, double yin, int zoneout, bool northpout)
        {
            int zone;
            double xout, yout;

            bool northp = northpin;
            if (zonein != zoneout)
            {
                // Determine lat, lon
                var (lat, lon) = Reverse(zonein, northpin, xin, yin);
                // Try converting to zoneout
                (zone, _, xout, yout) = Forward(lat, lon, (zoneout == (int)ZoneSpec.Match) ? zonein : zoneout);

                if (zone == 0 && northp != northpout)
                    throw new GeographicException("Attempt to transfer UPS coordinates between hemispheres");
            }
            else
            {
                if (zoneout == 0 && northp != northpout)
                    throw new GeographicException("Attempt to transfer UPS coordinates between hemispheres");

                (zone, xout, yout) = (zoneout, xin, yin);
            }

            if (northp != northpout)
                // Can't get here if UPS
                yout += (northpout ? -1 : 1) * MGRS.utmNshift_;

            return (zone, xout, yout);
        }

        /// <summary>
        /// Decode a UTM/UPS zone string.
        /// </summary>
        /// <param name="zonestr">string representation of zone and hemisphere.</param>
        /// <returns>
        /// <i>zone</i>, the UTM zone (zero means UPS) and <i>northp</i>, hemisphere (<see langword="true"/> means north, <see langword="false"/> means south).
        /// </returns>
        /// <remarks>
        /// For UTM, <paramref name="zonestr"/> has the form of a zone number in the range [1, 60] followed by a hemisphere letter,
        /// <c>n</c> or <c>s</c> (or "<c>north</c>" or "<c>south</c>" spelled out).
        /// For UPS, it consists just of the hemisphere letter (or the spelled out hemisphere).
        /// The returned value of <i>zone</i> is 0 for UPS.
        /// Note well that "<c>38s</c>" indicates the southern hemisphere of zone 38 and not latitude band <c>S</c>, 32° ≤ <i>lat</i> &lt; 40°.
        /// <c>n</c>, <c>01s</c>, <c>2n</c>, <c>38s</c>, <c>south</c>, <c>3north</c> are legal.
        /// <c>0n</c>, <c>001s</c>, <c>+3n</c>, <c>61n</c>, <c>38P</c> are illegal.
        /// <c>INV</c> is a special value for which the returned value of is <see cref="ZoneSpec.Invalid"/>.
        /// </remarks>
        public static (int zone, bool northp) DecodeZone(ReadOnlySpan<char> zonestr)
        {
            var zlen = zonestr.Length;
            if (zlen == 0)
                throw new GeographicException("Empty zone specification");

            // Longest zone spec is 32north, 42south, invalid = 7
            if (zlen > 7)
                throw new GeographicException("More than 7 characters in zone specification " + zonestr.ToString());

            var hemiStart = zonestr.FindFirstNotOf("1234567890");
            var startsWithNumber =
#if NETSTANDARD2_0
                int.TryParse(zonestr.Slice(0, hemiStart).ToString(), out var zone1);
#else
                int.TryParse(zonestr.Slice(0, hemiStart), out var zone1);
#endif
            // if (zone1 == 0) zone1 = UPS; (not necessary)

            if (zone1 == (int)ZoneSpec.UPS)
            {
                if (startsWithNumber)
                    // Don't allow 0n as an alternative to n for UPS coordinates
                    throw new GeographicException($"Illegal zone 0 in {zonestr.ToString()}, use just the hemisphere for UPS");
            }
            else if (!(zone1 >= (int)ZoneSpec.MinUTMZone && zone1 <= (int)ZoneSpec.MaxUTMZone))
                throw new GeographicException($"Zone {zone1} not in range [1, 60]");
            else if (!char.IsDigit(zonestr[0]))
                throw new GeographicException("Must use unsigned number for zone " + zone1);
            else if (hemiStart > 2)
                throw new GeographicException("More than 2 digits use to specify zone " + zone1);

            var hemi = zonestr.Slice(hemiStart);
            if (!startsWithNumber &&
                (hemi.CompareTo("inv".AsSpan(), StringComparison.OrdinalIgnoreCase) == 0 || hemi.CompareTo("invalid".AsSpan(), StringComparison.OrdinalIgnoreCase) == 0))
            {
                return ((int)ZoneSpec.Invalid, false);
            }

            bool northp1 = hemi.CompareTo("north".AsSpan(), StringComparison.OrdinalIgnoreCase) == 0 || hemi.CompareTo("n".AsSpan(), StringComparison.OrdinalIgnoreCase) == 0;
            if (!(northp1 || hemi.CompareTo("south".AsSpan(), StringComparison.OrdinalIgnoreCase) == 0 || hemi.CompareTo("s".AsSpan(), StringComparison.OrdinalIgnoreCase) == 0))
                throw new GeographicException($"Illegal hemisphere {hemi.ToString()} in {zonestr.ToString()}, specify north or south");

            return (zone1, northp1);
        }

        /// <summary>
        /// Encode a UTM/UPS zone string.
        /// </summary>
        /// <param name="zone">the UTM zone (zero means UPS).</param>
        /// <param name="northp">hemisphere (<see langword="true"/> means north, <see langword="false"/> means south).</param>
        /// <param name="abbrev">
        /// if <see langword="true"/> (the default) use abbreviated (<c>n</c>/<c>s</c>) notation for hemisphere;
        /// otherwise spell out the hemisphere (<c>north</c>/<c>south</c>).
        /// </param>
        /// <returns>string representation of zone and hemisphere.</returns>
        /// <remarks>
        /// <paramref name="zone"/> must be in the range [0, 60] with zone = 0, indicating UPS (but the resulting string does not contain "0").
        /// <paramref name="zone"/> may also be <see cref="ZoneSpec.Invalid"/>, in which case the returned string is "<c>inv</c>".
        /// This reverses <see cref="DecodeZone(ReadOnlySpan{char})"/>.
        /// </remarks>
        public static string EncodeZone(int zone, bool northp, bool abbrev = true)
        {
            if (zone == (int)ZoneSpec.Invalid)
                return abbrev ? "inv" : "invalid";
            if (!(zone >= (int)ZoneSpec.MinZone && zone <= (int)ZoneSpec.MaxZone))
                throw new GeographicException($"Zone {zone} not in range [0, 60]");

            return $"{(zone == (int)ZoneSpec.UPS ? "" : zone.ToString("D2"))}{(abbrev ? (northp ? "n" : "s") : (northp ? "north" : "south"))}";
        }

        /// <summary>
        /// Decode EPSG.
        /// </summary>
        /// <param name="epsg">the EPSG code.</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>zone</i>, the UTM zone (zero means UPS).</item>
        /// <item><i>northp</i>, hemisphere (<see langword="true"/> means north, <see langword="false"/> means south).</item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// EPSG (European Petroleum Survery Group) codes are a way to refer to many different projections.
        /// <see cref="DecodeEPSG"/> decodes those referring to UTM or UPS projections for the WGS84 ellipsoid.
        /// If the code does not refer to one of these projections, <i>zone</i> is set to <see cref="ZoneSpec.Invalid"/>.
        /// See <a href="https://www.spatialreference.org/ref/epsg/"></a>.
        /// </remarks>
        public static (int zone, bool northp) DecodeEPSG(int epsg)
        {
            var northp = false;
            var zone = (int)ZoneSpec.Invalid;

            if (epsg >= epsg01N && epsg <= epsg60N)
            {
                zone = (epsg - epsg01N) + (int)ZoneSpec.MinUTMZone;
                northp = true;
            }
            else if (epsg == epsgN)
            {
                zone = (int)ZoneSpec.UPS;
                northp = true;
            }
            else if (epsg >= epsg01S && epsg <= epsg60S)
            {
                zone = (epsg - epsg01S) + (int)ZoneSpec.MinUTMZone;
            }
            else if (epsg == epsgS)
            {
                zone = (int)ZoneSpec.UPS;
            }

            return (zone, northp);
        }

        /// <summary>
        /// Encode zone as EPSG.
        /// </summary>
        /// <param name="zone">the UTM zone (zero means UPS).</param>
        /// <param name="northp">hemisphere (true means north, false means south).</param>
        /// <returns>EPSG code (or -1 if zone is not in the range [0, 60])</returns>
        /// <remarks>
        /// Convert <paramref name="zone"/> and <paramref name="northp"/> to the corresponding EPSG (European Petroleum Survery Group) codes.
        /// </remarks>
        public static int EncodeEPSG(int zone, bool northp)
        {
            int epsg = -1;
            if (zone == (int)ZoneSpec.UPS)
                epsg = epsgS;
            else if (zone >= (int)ZoneSpec.MinUTMZone && zone <= (int)ZoneSpec.MaxUTMZone)
                epsg = (zone - (int)ZoneSpec.MinUTMZone) + epsg01S;
            if (epsg >= 0 && northp)
                epsg += epsgN - epsgS;
            return epsg;
        }
    }
}
