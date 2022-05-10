using System;
using System.Collections.Generic;
using System.Text;

using static System.Math;
using static GeographicLib.MathEx;

namespace GeographicLib.Geocodes
{
    /// <summary>
    /// Conversions for geohashes.
    /// </summary>
    /// <remarks>
    /// Geohashes are described in
    /// <list type="bullet">
    /// <item><a href="https://en.wikipedia.org/wiki/Geohash"></a></item>
    /// <item><a href="http://geohash.org/"></a></item>
    /// </list>
    /// <para>
    /// They provide a compact string representation of a particular geographic location (expressed as latitude and longitude),
    /// with the property that if trailing characters are dropped from the string the geographic location remains nearby.
    /// The classes <see cref="Georef"/> and <see cref="GARS"/> implement similar compact representations.
    /// </para>
    /// </remarks>
    public static class Geohash
    {
        private const int maxlen_ = 18;
        private const ulong mask_ = 1UL << 45;
        private const string lcdigits_ = "0123456789bcdefghjkmnpqrstuvwxyz";
        private const string ucdigits_ = "0123456789BCDEFGHJKMNPQRSTUVWXYZ";

        private static readonly double shift = Ldexp(1, 45),
                                       loneps = HD / shift,
                                       lateps = QD / shift;

        /// <summary>
        /// Convert from geographic coordinates to a geohash.
        /// </summary>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <param name="len">the length of the resulting geohash.</param>
        /// <returns>the geohash.</returns>
        /// <remarks>
        /// Internally, <i>len</i> is first put in the range [0, 18]. (<i>len</i> = 18 provides approximately 1μm precision.)
        /// <para>
        /// If <i>lat</i> or <i>lon</i> is <see cref="double.NaN"/>, the returned geohash is "invalid".
        /// </para>
        /// </remarks>
        public static string Forward(double lat, double lon, int len)
        {
            if (Abs(lat) > QD)
                throw new GeographicException($"Latitude {lat}d not in [-{QD}d, {QD}d]");

            if (double.IsNaN(lat) || double.IsNaN(lon))
            {
                return "invalid";
            }

            if (lat == QD) lat -= lateps / 2;
            lon = AngNormalize(lon);
            if (lon == HD) lon = -HD; // lon now in [-180,180)
                                        // lon/loneps in [-2^45,2^45); lon/loneps + shift in [0,2^46)
                                        // similarly for lat
            len = Max(0, Min(maxlen_, len));

            ulong
              ulon = (ulong)(Floor(lon / loneps) + shift),
              ulat = (ulong)(Floor(lat / lateps) + shift);

            Span<char> geohash1 = stackalloc char[maxlen_];
            uint byte_ = 0;

            for (uint i = 0; i < 5 * len;)
            {
                if ((i & 1) == 0)
                {
                    byte_ = (byte_ << 1) + (uint)(((ulon & mask_) != 0) ? 1 : 0);
                    ulon <<= 1;
                }
                else
                {
                    byte_ = (byte_ << 1) + (uint)(((ulon & mask_) != 0) ? 1 : 0); ;
                    ulat <<= 1;
                }
                ++i;
                if (i % 5 == 0)
                {
                    geohash1[(int)((i / 5) - 1)] = lcdigits_[(int)byte_];
                    byte_ = 0;
                }
            }

            return geohash1.Slice(0, len).ToString();
        }

        /// <summary>
        /// Convert from a geohash to geographic coordinates.
        /// </summary>
        /// <param name="geohash">the geohash.</param>
        /// <param name="centerp">if <see langword="true"/> (the default) return the center of the geohash location, otherwise return the south-west corner.</param>
        /// <returns>
        /// <i>lat</i>, latitude of point (degrees), <i>lon</i>, longitude of point (degrees) and <i>len</i>, the length of the geohash.
        /// </returns>
        /// <remarks>
        /// Only the first 18 characters for <i>geohash</i> are considered. (18 characters provides approximately 1μm precision.) 
        /// The case of the letters in <i>geohash</i> is ignored.
        /// <para>
        /// If the first 3 characters of <i>geohash</i> are "inv", then <i>lat</i> and <i>lon</i> are set to <see cref="double.NaN"/> and len is unchanged. ("nan" is treated similarly.)
        /// </para>
        /// </remarks>
        public static (double lat, double lon, int len) Reverse(ReadOnlySpan<char> geohash, bool centerp = true)
        {
            int len1 = Min(maxlen_, geohash.Length);
            if (len1 >= 3 &&
                (geohash.StartsWith("INV".AsSpan(), StringComparison.OrdinalIgnoreCase) || geohash.StartsWith("NAN".AsSpan(), StringComparison.OrdinalIgnoreCase))
               )
            {
                return (double.NaN, double.NaN, 0);
            }

            ulong ulon = 0, ulat = 0;
            for (uint k = 0, j = 0; k < len1; ++k)
            {
                int byte_ = ucdigits_.IndexOf(geohash[(int)k]);
                if (byte_ < 0)
                    throw new GeographicException("Illegal character in geohash: " + geohash.ToString());
                for (uint m = 16; m != 0; m >>= 1)
                {
                    if (j == 0)
                        ulon = (ulon << 1) + (uint)(((byte_ & m) != 0) ? 1 : 0);
                    else
                        ulat = (ulat << 1) + (uint)(((byte_ & m) != 0) ? 1 : 0);
                    j ^= 1;
                }
            }
            ulon <<= 1; ulat <<= 1;
            if (centerp)
            {
                ulon += 1;
                ulat += 1;
            }
            int s = 5 * (maxlen_ - len1);
            ulon <<= (s / 2);
            ulat <<= s - (s / 2);

            return (ulat * lateps - QD, ulon * loneps - HD, len1);
        }

        /// <summary>
        /// Gets the latitude resolution of a geohash.
        /// </summary>
        /// <param name="len">the length of the geohash.</param>
        /// <returns>the latitude resolution (degrees).</returns>
        /// <remarks>
        /// Internally, <paramref name="len"/> is first put in the range [0, 18].
        /// </remarks>
        public static double LatitudeResolution(int len)
        {
            len = Max(0, Min(maxlen_, len));
            return Ldexp(HD, -(5 * len / 2));
        }

        /// <summary>
        /// Gets the longitude resolution of a geohash.
        /// </summary>
        /// <param name="len">the length of the geohash.</param>
        /// <returns>the longitude resolution (degrees).</returns>
        /// <remarks>
        /// Internally, <paramref name="len"/> is first put in the range [0, 18].
        /// </remarks>
        public static double LongitudeResolution(int len)
        {
            len = Max(0, Min(maxlen_, len));
            return Ldexp(TD, -(5 * len - 5 * len / 2));
        }

        /// <summary>
        /// Gets the geohash length required to meet a given geographic resolution.
        /// </summary>
        /// <param name="res">the minimum of resolution in latitude and longitude (degrees).</param>
        /// <returns>geohash length in the range [0, 18].</returns>
        public static double GeohashLength(int res)
        {
            res = Abs(res);
            for (int len = 0; len < maxlen_; ++len)
                if (LongitudeResolution(len) <= res)
                    return len;
            return maxlen_;
        }

        /// <summary>
        /// Gets the geohash length required to meet a given geographic resolution.
        /// </summary>
        /// <param name="latres">the resolution in latitude (degrees).</param>
        /// <param name="lonres">the resolution in longitude (degrees).</param>
        /// <returns>geohash length in the range [0, 18].</returns>
        public static double GeohashLength(int latres, int lonres)
        {
            latres = Abs(latres);
            lonres = Abs(lonres);
            for (int len = 0; len < maxlen_; ++len)
                if (LatitudeResolution(len) <= latres &&
                    LongitudeResolution(len) <= lonres)
                    return len;
            return maxlen_;
        }

        /// <summary>
        /// Gets the decimal geographic precision required to match a given geohash length.
        /// This is the number of digits needed after decimal point in a decimal degrees representation.
        /// </summary>
        /// <param name="len">the length of the geohash.</param>
        /// <returns>the decimal precision (may be negative).</returns>
        /// <remarks>
        /// Internally, <paramref name="len"/> is first put in the range [0, 18]. The returned decimal precision is in the range [−2, 12].
        /// </remarks>
        public static int DecimalPrecision(int len)
            => -(int)Floor(Log(LatitudeResolution(len)) / Log(10d));
    }
}
