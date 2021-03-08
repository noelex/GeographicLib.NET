using GeographicLib.Geocodes;
using System;
using System.Collections.Generic;
using System.Text;

using static System.Math;
using static GeographicLib.MathEx;

namespace GeographicLib
{
    /// <summary>
    /// Provide conversion between geographic coordinates.
    /// </summary>
    /// <remarks>
    /// This class stores a geographic position which may be set via the constructors or <see cref="Reset(double, double, int)"/> via
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
    public class GeoCoords : IEllipsoid
    {
        private double _lat, _long, _easting, _northing, _gamma, _k;
        private bool _northp;
        private int _zone;                  // See UTMUPS::zonespec

        private double _alt_easting, _alt_northing, _alt_gamma, _alt_k;
        private int _alt_zone;

        /// <summary>
        /// Initialize a new <see cref="GeoCoords"/> instance from a string.
        /// </summary>
        /// <param name="s">1-element, 2-element, or 3-element string representation of the position.</param>
        /// <param name="centerp">governs the interpretation of <see cref="MGRS"/> coordinates (see below).</param>
        /// <param name="longfirst">governs the interpretation of geographic coordinates (see below).</param>
        /// <remarks>
        /// Parse as a string and interpret it as a geographic position.
        /// The input string is broken into space (or comma) separated pieces and Basic decision on which format is based on number of components
        /// <list type="number">
        /// <item>"Lat Long" or "Long Lat"</item>
        /// <item><see cref="MGRS"/></item>
        /// <item>"Zone Easting Northing" or "Easting Northing Zone"</item>
        /// </list>
        /// The following inputs are approximately the same (Ar Ramadi Bridge, Iraq)
        /// <list type="bullet">
        /// <item>
        /// Latitude and Longitude
        /// <code>
        /// 33.44 43.27, N33d26.4' E43d16.2', 43d16'12"E 33d26'24"N, 43:16:12E 33:26:24
        /// </code>
        /// </item>
        /// <item>
        /// <see cref="MGRS"/>
        /// <code>
        /// 38SLC30, 38SLC391014, 38SLC3918701405, 37SHT9708
        /// </code>
        /// </item>
        /// <item>
        /// UTM
        /// <code>
        /// 38n 339188 3701405, 897039 3708229 37n
        /// </code>
        /// </item>
        /// </list>
        /// <b>Latitude and Longitude parsing</b>:
        /// Latitude precedes longitude, unless a N, S, E, W hemisphere designator is used on one or both coordinates.
        /// If <paramref name="longfirst"/> = <see langword="true"/> (default is <see langword="false"/>),
        /// then longitude precedes latitude in the absence of a hemisphere designator.
        /// Thus (with <paramref name="longfirst"/> = <see langword="false"/>)
        /// <code>
        /// 40 -75, N40 W75, -75 N40, 75W 40N, E-75 -40S
        /// </code>
        /// are all the same position. The coordinates may be given in decimal degrees, degrees and decimal minutes,
        /// degrees, minutes, seconds, etc. Use d, ', and " to mark off the degrees, minutes and seconds.
        /// Various alternative symbols for degrees, minutes, and seconds are allowed.
        /// Alternatively, use : to separate these components. A single addition or subtraction is allowed.
        /// (See <see cref="DMS.Decode(string)"/> for details.) Thus
        /// <code>
        /// 40d30'30", 40d30'30, 40°30'30, 40d30.5', 40d30.5, 40:30:30, 40:30.5, 40.508333333, 40:30+0:0:30, 40:31-0:0.5
        /// </code>
        /// all specify the same angle.The leading sign applies to the following components so -1d30 is -(1+30/60) = −1.5.
        /// However, note that -1:30-0:0:15 is parsed as (-1:30) + (-0:0:15) = −(1+30/60) − (15/3600).
        /// Latitudes must be in the range[−90°, 90°]. Internally longitudes are reduced to the range[−180°, 180°].
        /// <para>
        /// <b>UTM/UPS parsing</b>: For UTM zones (−80° ≤ Lat &lt; 84°),
        /// the zone designator is made up of a zone number (for 1 to 60) and a hemisphere letter (n or s),
        /// e.g., 38n (38north can also be used).
        /// The latitude band designer ([C–M] in the southern hemisphere and [N–X] in the northern) should NOT be used.
        /// (This is part of the <see cref="MGRS"/> coordinate.)
        /// The zone designator for the poles (where UPS is employed) is a hemisphere letter by itself,
        /// i.e., n or s (north or south can also be used).
        /// </para>
        /// <para>
        /// <b>MGRS parsing</b> interprets the grid references as square area at the specified precision
        /// (1m, 10m, 100m, etc.). If <paramref name="centerp"/> = <see langword="true"/> (the default),
        /// the center of this square is then taken to be the precise position; thus:
        /// <list type="bullet">
        /// <item><c>38SMB = 38n 450000 3650000</c></item>
        /// <item><c>38SMB4484 = 38n 444500 3684500</c></item>
        /// <item><c>38SMB44148470 = 38n 444145 3684705</c></item>
        /// </list>
        /// Otherwise, the "south-west" corner of the square is used, i.e.,
        /// <list type="bullet">
        /// <item><c>38SMB = 38n 400000 3600000</c></item>
        /// <item><c>38SMB4484 = 38n 444000 3684000</c></item>
        /// <item><c>38SMB44148470 = 38n 444140 3684700</c></item>
        /// </list>
        /// </para>
        /// </remarks>
        public GeoCoords(string s, bool centerp = true, bool longfirst = false)
            => Reset(s, centerp, longfirst);

        /// <summary>
        /// Initialize a new <see cref="GeoCoords"/> instance from geographic coordinates.
        /// </summary>
        /// <param name="latitude">latitude in degrees.</param>
        /// <param name="longitude">longitude in degrees.</param>
        /// <param name="zone">
        /// if specified, force the UTM/UPS representation to use a specified zone using the rules given in <see cref="ZoneSpec"/>.
        /// </param>
        public GeoCoords(double latitude, double longitude, int zone = (int)ZoneSpec.Standard)
            => Reset(latitude, longitude, zone);

        /// <summary>
        /// Initialize a new <see cref="GeoCoords"/> instance from UTM/UPS coordinates.
        /// </summary>
        /// <param name="zone">UTM zone (zero means UPS).</param>
        /// <param name="northp">hemisphere (<see langword="true"/> means north, <see langword="false"/> means south).</param>
        /// <param name="easting">easting in meters.</param>
        /// <param name="northing">northing in meters.</param>
        public GeoCoords(int zone, bool northp, double easting, double northing)
            => Reset(zone, northp, easting, northing);

        /// <summary>
        /// Gets a value representing latitude of the coordinate in degrees.
        /// </summary>
        public double Latitude => _lat;

        /// <summary>
        /// Gets a value representing longitude of the coordinate in degrees.
        /// </summary>
        public double Longitude => _long;

        /// <summary>
        /// Gets a value representing easting of the coordinate in meters.
        /// </summary>
        public double Easting => _easting;

        /// <summary>
        /// Gets a value representing northing of the coordinate in meters.
        /// </summary>
        public double Northing => _northing;

        /// <summary>
        /// Gets a value representing meridian convergence (degrees) for the UTM/UPS projection.
        /// </summary>
        public double Convergence => _gamma;

        /// <summary>
        /// Gets a value representing scale for the UTM/UPS projection.
        /// </summary>
        public double Scale => _k;

        /// <summary>
        /// Gets a value representing whether current coordinate is in north hemiphere.
        /// </summary>
        public bool IsNorthHemisphere => _northp;

        /// <summary>
        /// Gets the hemisphere letter <c>n</c> or <c>s</c>.
        /// </summary>
        public char Hemisphere => _northp ? 'n' : 's';

        /// <summary>
        /// Gets a value representing the zone corresponding to the input (return 0 for UPS).
        /// </summary>
        public int Zone => _zone;

        /// <summary>
        /// Gets a value representing current alternate zone (return 0 for UPS).
        /// </summary>
        public int AltZone => _alt_zone;

        /// <summary>
        /// Gets a value representing easting (meters) for alternate zone.
        /// </summary>
        public double AltEasting => _alt_easting;

        /// <summary>
        /// Gets a value representing northing (meters) for alternate zone.
        /// </summary>
        public double AltNorthing => _alt_northing;

        /// <summary>
        /// Gets a value representing meridian convergence (degrees) for alternate zone.
        /// </summary>
        public double AltConvergence => _alt_gamma;

        /// <summary>
        /// Gets a value representing scale for alternate zone.
        /// </summary>
        public double AltScale => _alt_k;

        /// <inheritdoc/>
        public double EquatorialRadius => UTMUPS.EquatorialRadius;

        /// <inheritdoc/>
        public double Flattening => UTMUPS.Flattening;

        /// <summary>
        /// Reset the location in terms of geographic coordinates. See <see cref="GeoCoords.GeoCoords(double, double, int)"/>.
        /// </summary>
        /// <param name="latitude">latitude in degrees.</param>
        /// <param name="longitude">longitude in degrees.</param>
        /// <param name="zone">
        /// if specified, force the UTM/UPS representation to use a specified zone using the rules given in <see cref="ZoneSpec"/>.
        /// </param>
        public void Reset(double latitude, double longitude, int zone = (int)ZoneSpec.Standard)
        {
            (_zone, _northp, _easting, _northing) = UTMUPS.Forward(latitude, longitude, out _gamma, out _k, zone);
            _lat = latitude;
            _long = longitude;
            if (_long >= 180) _long -= 360;
            else if (_long < -180) _long += 360;
            CopyToAlt();
        }

        /// <summary>
        /// Reset the location in terms of UPS/UPS coordinates. See <see cref="GeoCoords.GeoCoords(int, bool, double, double)"/>.
        /// </summary>
        /// <param name="zone">UTM zone (zero means UPS).</param>
        /// <param name="northp">hemisphere (<see langword="true"/> means north, <see langword="false"/> means south).</param>
        /// <param name="easting">easting in meters.</param>
        /// <param name="northing">northing in meters.</param>
        public void Reset(int zone, bool northp, double easting, double northing)
        {
            _zone = zone;
            _northp = northp;
            _easting = easting;
            _northing = northing;
            FixHemisphere();
            CopyToAlt();
        }

        /// <summary>
        /// Reset the location from a string. See <see cref="GeoCoords(string, bool, bool)"/>.
        /// </summary>
        /// <param name="s">1-element, 2-element, or 3-element string representation of the position.</param>
        /// <param name="centerp">governs the interpretation of <see cref="MGRS"/> coordinates.</param>
        /// <param name="longfirst">governs the interpretation of geographic coordinates.</param>
        public void Reset(string s, bool centerp = true, bool longfirst = false)
        {
            var sa = s.Split(" \t\n\v\f\r,".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);

            if (sa.Length == 1)
            {
                int prec;
                (_zone, _northp, _easting, _northing, prec) = MGRS.Reverse(sa[0].AsSpan(), centerp);
                (_lat, _long) = UTMUPS.Reverse(_zone, _northp, _easting, _northing, out _gamma, out _k);
            }
            else if (sa.Length == 2)
            {
                (_lat, _long) = DMS.Decode(sa[0], sa[1], longfirst);
                _long = AngNormalize(_long);
                (_zone, _northp, _easting, _northing) = UTMUPS.Forward(_lat, _long , out _gamma, out _k);
            }
            else if (sa.Length == 3)
            {
                int zoneind, coordind;
                if (sa[0].Length > 0 && char.IsLetter(sa[0][sa[0].Length - 1]))
                {
                    zoneind = 0;
                    coordind = 1;
                }
                else if (sa[2].Length > 0 && char.IsLetter(sa[2][sa[2].Length - 1]))
                {
                    zoneind = 2;
                    coordind = 0;
                }
                else
                    throw new GeographicException("Neither " + sa[0] + " nor " + sa[2]
                                        + " of the form UTM/UPS Zone + Hemisphere"
                                        + " (ex: 38n, 09s, n)");

                (_zone, _northp) = UTMUPS.DecodeZone(sa[zoneind].AsSpan());
                for (int i = 0; i < 2; ++i)
                {
                    if (i != 0)
                    {
                        _northing = sa[coordind + i].ParseDouble();
                    }
                    else
                    {
                        _easting = sa[coordind + i].ParseDouble();
                    }
                }

                (_lat, _long) = UTMUPS.Reverse(_zone, _northp, _easting, _northing, out _gamma, out _k);
                FixHemisphere();
            }
            else
                throw new GeographicException("Coordinate requires 1, 2, or 3 elements");
            CopyToAlt();
        }

        /// <summary>
        /// Specify alternate zone number.
        /// </summary>
        /// <param name="zone">zone number for the alternate representation.</param>
        /// <remarks>
        /// See <see cref="ZoneSpec"/> for more information on the interpretation of zone.
        /// Note that <paramref name="zone"/> == <see cref="ZoneSpec.Standard"/> (the default) use the standard UPS or UTM zone,
        /// <see cref="ZoneSpec.Match"/> does nothing retaining the existing alternate representation.
        /// Before this is called the alternate zone is the input zone.
        /// </remarks>
        public void SetAltZone(int zone = (int)ZoneSpec.Standard)
        {
            if (zone == (int)ZoneSpec.Match)
                return;

            zone = UTMUPS.StandardZone(_lat, _long, zone);
            if (zone == _zone)
                CopyToAlt();
            else
            {
                (_alt_zone, _, _alt_easting, _alt_northing) = UTMUPS.Forward(_lat, _long, out _alt_gamma, out _alt_k, zone);
            }
        }

        /// <summary>
        /// String representation with latitude and longitude as signed decimal degrees.
        /// </summary>
        /// <param name="prec">precision (relative to about 1m).</param>
        /// <param name="longfirst">if <see langword="true"/> give longitude first (default = <see langword="false"/>)</param>
        /// <returns>decimal latitude/longitude string representation.</returns>
        /// <remarks>
        /// Precision specifies accuracy of representation as follows:
        /// <list type="bullet">
        /// <item><paramref name="prec"/> = −5 (min), 1°</item>
        /// <item><paramref name="prec"/> = 0, 10−5° (about 1m)</item>
        /// <item><paramref name="prec"/> = 3, 10−8°</item>
        /// <item><paramref name="prec"/> = 9 (max), 10−14°</item>
        /// </list>
        /// </remarks>
        public string ToGeoString(int prec = 0, bool longfirst = false)
        {
            prec = Max(0, Min(9 + 0 /*Math::extra_digits()*/, prec) + 5);
            var os = new StringBuilder();
            var format = $"F{prec}";//os << fixed << setprecision(prec);

            double a = longfirst ? _long : _lat;
            double b = longfirst ? _lat : _long;
            if (!double.IsNaN(a))
                os.AppendFormat(format, a);
            else
                os.Append("nan");
            os.Append(" ");
            if (!double.IsNaN(b))
                os.AppendFormat(format, b);
            else
                os.Append("nan");
            return os.ToString();
        }

        /// <summary>
        /// String representation with latitude and longitude as degrees, minutes, seconds, and hemisphere.
        /// </summary>
        /// <param name="prec">precision (relative to about 1m)</param>
        /// <param name="longfirst">if <see langword="true"/> give longitude first (default = <see langword="false"/>)</param>
        /// <param name="dmssep">if not <c>0</c>, use as the <see cref="DMS"/> separator character (instead of d, ', " delimiters).</param>
        /// <returns><see cref="DMS"/> latitude/longitude string representation.</returns>
        /// <remarks>
        /// Precision specifies accuracy of representation as follows:
        /// <list type="bullet">
        /// <item><paramref name="prec"/> = −5 (min), 1°</item>
        /// <item><paramref name="prec"/> = −4, 0.1°</item>
        /// <item><paramref name="prec"/> = −3, 1'</item>
        /// <item><paramref name="prec"/> = −2, 0.1'</item>
        /// <item><paramref name="prec"/> = −1, 1"</item>
        /// <item><paramref name="prec"/> = 0, 0.1" (about 3m)</item>
        /// <item><paramref name="prec"/> = 1, 0.01"</item>
        /// <item><paramref name="prec"/> = 10 (max), 10−11"</item>
        /// </list>
        /// </remarks>
        public string ToDMSString(int prec = 0, bool longfirst = false, char dmssep = '\0')
        {
            prec = Max(0, Min(10 + 0, prec) + 5);
            return DMS.Encode(longfirst ? _long : _lat, prec,
                               longfirst ? HemisphereIndicator.Longitude : HemisphereIndicator.Latitude, dmssep) +
              " " + DMS.Encode(longfirst ? _lat : _long, prec,
                               longfirst ? HemisphereIndicator.Latitude : HemisphereIndicator.Longitude, dmssep);
        }

        /// <summary>
        /// Gets <see cref="MGRS"/> string of the coordinate.
        /// </summary>
        /// <param name="prec">precision (relative to about 1m).</param>
        /// <returns>A <see cref="MGRS"/> string.</returns>
        /// <remarks>
        /// This gives the coordinates of the enclosing grid square with size given by the precision.
        /// Thus <c>38n 444180 3684790</c> converted to a <see cref="MGRS"/> coordinate at precision −2 (100m) is
        /// <c>38SMB441847</c> and not <c>38SMB442848</c>. <paramref name="prec"/> specifies the precision of
        /// the <see cref="MGRS"/> string as follows:
        /// <list type="bullet">
        /// <item><paramref name="prec"/> = −6 (min), only the grid zone is returned, e.g., <c>38S</c></item>
        /// <item><paramref name="prec"/> = −5, 100km, e.g., <c>38SMB</c></item>
        /// <item><paramref name="prec"/> = −4, 10km</item>
        /// <item><paramref name="prec"/> = −3, 1km</item>
        /// <item><paramref name="prec"/> = −2, 100m</item>
        /// <item><paramref name="prec"/> = −1, 10m</item>
        /// <item><paramref name="prec"/> = 0, 1m</item>
        /// <item><paramref name="prec"/> = 1, 0.1m</item>
        /// <item><paramref name="prec"/> = 6 (max), 1μm</item>
        /// </list>
        /// </remarks>
        public string ToMGRSString(int prec = 0)
            => MGRS.Forward(_zone, _northp, _easting, _northing, _lat, Max(-1, Min(6, prec) + 5));

        /// <summary>
        /// Gets UTM/UPS string of the coordinate.
        /// </summary>
        /// <param name="prec">precision (relative to about 1m).</param>
        /// <param name="abbrev">if <see langword="true"/> (the default) use abbreviated (<c>n</c>/<c>s</c>) notation for hemisphere;
        /// otherwise spell out the hemisphere (<c>north</c>/<c>south</c>).</param>
        /// <returns>UTM/UPS string representation: zone designator, easting, and northing.</returns>
        /// <remarks>
        /// Precision specifies accuracy of representation as follows:
        /// <list type="bullet">
        /// <item><paramref name="prec"/> = −5, 100km</item>
        /// <item><paramref name="prec"/> = −3, 1km</item>
        /// <item><paramref name="prec"/> = 0, 1m</item>
        /// <item><paramref name="prec"/> = 3, 1mm</item>
        /// <item><paramref name="prec"/> = 6, 1μm</item>
        /// <item><paramref name="prec"/> = 9 (max), 1nm</item>
        /// </list>
        /// </remarks>
        public string ToUTMUPSString(int prec = 0, bool abbrev = true)
            => UTMUPSString(_zone, _northp, _easting, _northing, prec, abbrev);

        /// <summary>
        /// Gets UTM/UPS string with hemisphere override.
        /// </summary>
        /// <param name="northp">hemisphere override</param>
        /// <param name="prec">precision (relative to about 1m)</param>
        /// <param name="abbrev">if <see langword="true"/> (the default) use abbreviated (<c>n</c>/<c>s</c>) notation for hemisphere;
        /// otherwise spell out the hemisphere (<c>north</c>/<c>south</c>).</param>
        /// <returns>UTM/UPS string representation: zone designator, easting, and northing.</returns>
        public string ToUTMUPSString(bool northp, int prec = 0, bool abbrev = true)
        {
            var (e, n, _) = UTMUPS.Transfer(_zone, _northp, _easting, _northing, _zone, northp);
            return UTMUPSString(_zone, northp, e, n, prec, abbrev);
        }

        /// <summary>
        /// Gets <see cref="MGRS"/> string for the alternate zone. See <see cref="ToMGRSString(int)"/>.
        /// </summary>
        /// <param name="prec">precision (relative to about 1m).</param>
        /// <returns>A <see cref="MGRS"/> string.</returns>
        public string ToAltMGRSString(int prec = 0)
            => MGRS.Forward(_alt_zone, _northp, _alt_easting, _alt_northing, _lat, Max(-1, Min(6, prec) + 5));

        /// <summary>
        /// Gets UTM/UPS string for the alternate zone. See <see cref="ToUTMUPSString(int, bool)"/>.
        /// </summary>
        /// <param name="prec">precision (relative to about 1m).</param>
        /// <param name="abbrev">if <see langword="true"/> (the default) use abbreviated (<c>n</c>/<c>s</c>) notation for hemisphere;
        /// otherwise spell out the hemisphere (<c>north</c>/<c>south</c>).</param>
        /// <returns>UTM/UPS string representation: zone designator, easting, and northing.</returns>
        public string ToAltUTMUPSString(int prec = 0, bool abbrev = true)
            => UTMUPSString(_alt_zone, _northp, _alt_easting, _alt_northing, prec, abbrev);

        /// <summary>
        /// UTM/UPS string for the alternate zone, with hemisphere override.
        /// </summary>
        /// <param name="northp">hemisphere override</param>
        /// <param name="prec">precision (relative to about 1m)</param>
        /// <param name="abbrev">if <see langword="true"/> (the default) use abbreviated (<c>n</c>/<c>s</c>) notation for hemisphere;
        /// otherwise spell out the hemisphere (<c>north</c>/<c>south</c>).</param>
        /// <returns>UTM/UPS string representation: zone designator, easting, and northing.</returns>
        public string ToAltUTMUPSString(bool northp, int prec = 0, bool abbrev = true)
        {
            var (e, n, _) = UTMUPS.Transfer(_alt_zone, _northp, _alt_easting, _alt_northing, _alt_zone, northp);
            return UTMUPSString(_alt_zone, northp, e, n, prec, abbrev);
        }

        static string UTMUPSString(int zone, bool northp,
                     double easting, double northing,
                     int prec, bool abbrev)
        {
            var os = new StringBuilder();
            prec = Max(-5, Min(9 + 0/*Math::extra_digits()*/, prec));
            // Need extra real because, since C++11, pow(float, int) returns double
            var scale = prec < 0 ? Pow(10, -prec) : 1;
            os.Append(UTMUPS.EncodeZone(zone, northp, abbrev));// << fixed << setfill('0');

            if (IsFinite(easting))
            {
                os.Append(" ");
                os.Append((easting / scale).ToFixedString(Max(0, prec)));
                if (prec < 0 && Abs(easting / scale) > 0.5)
                    os.Append("0".PadLeft(-prec, '0'));
            }
            else
                os.Append(" nan");

            if (IsFinite(northing))
            {
                os.Append((northing / scale).ToFixedString(Max(0, prec)));
                if (prec < 0 && Abs(northing / scale) > 0.5)
                    os.Append("0".PadLeft(-prec, '0'));
            }
            else
                os.Append(" nan");

            return os.ToString();
        }

        private void FixHemisphere()
        {
            if (_lat == 0 || (_northp && _lat >= 0) || (!_northp && _lat < 0) || double.IsNaN(_lat))
                // Allow either hemisphere for equator
                return;
            if (_zone != (int)ZoneSpec.UPS)
            {
                _northing += (_northp ? 1 : -1) * UTMUPS.UTMShift;
                _northp = !_northp;
            }
            else
                throw new GeographicException("Hemisphere mixup");
        }

        private void CopyToAlt() =>
            (_alt_easting, _alt_northing, _alt_gamma, _alt_k, _alt_zone) = (_easting, _northing, _gamma, _k, _zone);

    }
}
