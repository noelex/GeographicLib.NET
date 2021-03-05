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
        /// </remarks>
        public GeoCoords(string s, bool centerp = true, bool longfirst = false)
            => Reset(s, centerp, longfirst);

        public GeoCoords(double latitude, double longitude, int zone = (int)ZoneSpec.Standard)
            => Reset(latitude, longitude, zone);

        public GeoCoords(int zone, bool northp, double easting, double northing)
            => Reset(zone, northp, easting, northing);

        public double Latitude => _lat;

        public double Longitude => _long;

        public double Easting => _easting;

        public double Northing => _northing;

        public double Convergence => _gamma;

        public double Scale => _k;

        public bool IsNorthHemisphere => _northp;

        public char Hemisphere => _northp ? 'n' : 's';

        public int Zone => _zone;

        public int AltZone => _alt_zone;

        public double AltEasting => _alt_easting;

        public double AltNorthing => _alt_northing;

        public double AltConvergence => _alt_gamma;

        public double AltScale => _alt_k;

        /// <inheritdoc/>
        public double EquatorialRadius => UTMUPS.EquatorialRadius;

        /// <inheritdoc/>
        public double Flattening => UTMUPS.Flattening;

        public void Reset(double latitude, double longitude, int zone = (int)ZoneSpec.Standard)
        {
            _lat = latitude;
            _long = longitude;
            if (_long >= 180) _long -= 360;
            else if (_long < -180) _long += 360;
            CopyToAlt();
        }

        public void Reset(int zone, bool northp, double easting, double northing)
        {
            _zone = zone;
            _northp = northp;
            _easting = easting;
            _northing = northing;
            FixHemisphere();
            CopyToAlt();
        }

        public void Reset(string s, bool centerp = true, bool longfirst = false)
        {
            var sa = new List<string>();
            const string spaces = " \t\n\v\f\r,"; // Include comma as a space

            for (int pos0 = 0, pos1; pos0 != -1;)
            {
                pos1 = s.FindFirstNotOf(spaces, pos0);
                if (pos1 == -1)
                    break;
                pos0 = s.IndexOf(spaces, pos1);
                sa.Add(s.Substring(pos1, pos0 == -1 ? pos0 : pos0 - pos1));
            }
            if (sa.Count == 1)
            {
                int prec;
                (_zone, _northp, _easting, _northing, prec) = MGRS.Reverse(sa[0].AsSpan(), centerp);
                (_lat, _long) = UTMUPS.Reverse(_zone, _northp, _easting, _northing, out _gamma, out _k);
            }
            else if (sa.Count == 2)
            {
                throw new NotSupportedException("DMS string input is not supported yet.");
                //DMS::DecodeLatLon(sa[0], sa[1], _lat, _long, longfirst);
                //_long = Math::AngNormalize(_long);
                //UTMUPS::Forward(_lat, _long,
                //                 _zone, _northp, _easting, _northing, _gamma, _k);
            }
            else if (sa.Count == 3)
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

        public string ToDMSString(int prec = 0, bool longfirst = false, char dmssep = '\0')
            => throw new NotImplementedException();

        public string ToMGRSString(int prec = 0)
            => MGRS.Forward(_zone, _northp, _easting, _northing, _lat, Max(-1, Min(6, prec) + 5));

        public string ToUTMUPSString(int prec = 0, bool abbrev = true)
            => UTMUPSString(_zone, _northp, _easting, _northing, prec, abbrev);

        public string ToUTMUPSString(bool northp, int prec = 0, bool abbrev = true)
        {
            var (e, n, _) = UTMUPS.Transfer(_zone, _northp, _easting, _northing, _zone, northp);
            return UTMUPSString(_zone, northp, e, n, prec, abbrev);
        }

        public string ToAltMGRSString(int prec = 0)
            => MGRS.Forward(_alt_zone, _northp, _alt_easting, _alt_northing, _lat, Max(-1, Min(6, prec) + 5));

        public string ToAltUTMUPSString(int prec = 0, bool abbrev = true)
            => UTMUPSString(_alt_zone, _northp, _alt_easting, _alt_northing, prec, abbrev);

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
