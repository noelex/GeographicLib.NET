using System;
using System.Collections.Generic;
using System.Text;

using static System.Math;
using static GeographicLib.MathEx;

namespace GeographicLib.Geocodes
{
    /// <summary>
    /// Conversions for the World Geographic Reference System (georef).
    /// </summary>
    /// <remarks>
    /// The World Geographic Reference System is described in
    /// <list type="bullet">
    /// <item><a href="https://en.wikipedia.org/wiki/Georef"></a></item>
    /// <item><a href="http://earth-info.nga.mil/GandG/coordsys/grids/georef.pdf"></a></item>
    /// </list>
    /// <para>
    /// It provides a compact string representation of a geographic area (expressed as latitude and longitude).
    /// The classes <see cref="GARS"/> and <see cref="Geohash"/> implement similar compact representations.
    /// </para>
    /// </remarks>
    public static class Georef
    {
        private const string digits_ = "0123456789";
        private const string lontile_ = "ABCDEFGHJKLMNPQRSTUVWXYZ";
        private const string lattile_ = "ABCDEFGHJKLM";
        private const string degrees_ = "ABCDEFGHJKLMNPQ";

        private const int tile_ = 15,               // The size of tile in degrees
                          lonorig_ = -180,          // Origin for longitude
                          latorig_ = -90,           // Origin for latitude
                          base_ = 10,               // Base for minutes
                          baselen_ = 4,
                          maxprec_ = 11,            // approximately equivalent to MGRS class
                          maxlen_ = baselen_ + 2 * maxprec_;
        /// <summary>
        /// Convert from geographic coordinates to <see cref="Georef"/>.
        /// </summary>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <param name="prec">the precision of the resulting georef.</param>
        /// <returns>the georef string.</returns>
        /// <remarks>
        /// <i>prec</i> specifies the precision of <i>georef</i> as follows:
        /// <list type="bullet">
        /// <item><i>prec</i> = −1 (min), 15°</item>
        /// <item><i>prec</i> = 0, 1°</item>
        /// <item><i>prec</i> = 1, converted to <i>prec</i> = 2</item>
        /// <item><i>prec</i> = 2, 1'</item>
        /// <item><i>prec</i> = 3, 0.1'</item>
        /// <item><i>prec</i> = 4, 0.01'</item>
        /// <item><i>prec</i> = 5, 0.001'</item>
        /// <item>…</item>
        /// <item><i>prec</i> = 11 (max), 10^−9'</item>
        /// </list>
        /// <para>
        /// If <paramref name="lat"/> or <paramref name="lon"/> is <see cref="double.NaN"/>, then georef is set to "INVALID".
        /// </para>
        /// </remarks>
        public static string Forward(double lat, double lon, int prec)
        {
            if (Abs(lat) > 90)
                throw new GeographicException($"Latitude {lat}d not in [-90d, 90d]");
            if (double.IsNaN(lat) || double.IsNaN(lon))
            {
                return "INVALID";
            }

            lon = AngNormalize(lon); // lon in [-180,180)
            if (lat == 90) lat *= (1 - DBL_EPSILON / 2);
            prec = Max(-1, Min(maxprec_, prec));

            if (prec == 1) ++prec;      // Disallow prec = 1
                                        // The C++ standard mandates 64 bits for long long.  But
                                        // check, to make sure.

            var m = 60000000000L;
            long
              x = (long)Floor(lon * m) - lonorig_ * m,
              y = (long)Floor(lat * m) - latorig_ * m;
            int ilon = (int)(x / m); int ilat = (int)(y / m);
            Span<char> georef1 = stackalloc char[maxlen_];
            georef1[0] = lontile_[ilon / tile_];
            georef1[1] = lattile_[ilat / tile_];
            if (prec >= 0)
            {
                georef1[2] = degrees_[ilon % tile_];
                georef1[3] = degrees_[ilat % tile_];
                if (prec > 0)
                {
                    x -= m * ilon; y -= m * ilat;
                    long d = (long)Pow(base_, maxprec_ - prec);
                    x /= d; y /= d;
                    for (int c = prec; c-- > 0;)
                    {
                        georef1[baselen_ + c] = digits_[(int)(x % base_)]; x /= base_;
                        georef1[baselen_ + c + prec] = digits_[(int)(y % base_)]; y /= base_;
                    }
                }
            }

            return georef1.Slice(0, baselen_ + 2 * prec).ToString();
        }

        /// <summary>
        /// Convert from <see cref="Georef"/> to geographic coordinates.
        /// </summary>
        /// <param name="georef">the <see cref="Georef"/> string.</param>
        /// <param name="centerp">
        /// if <see langword="true"/> (the default) return the center <i>georef</i>, otherwise return the south-west corner.
        /// </param>
        /// <returns><i>lat</i>, latitude of point (degrees), <i>lon</i>, longitude of point (degrees) and <i>prec</i>, 
        /// precision of the input <see cref="Georef"/> string.</returns>
        /// <remarks>
        /// The case of the letters in <i>georef</i> is ignored. <i>prec</i> is in the range [−1, 11] and gives the precision of <i>georef</i> as follows:
        /// <list type="bullet">
        /// <item><i>prec</i> = −1 (min), 15°</item>
        /// <item><i>prec</i> = 0, 1°</item>
        /// <item><i>prec</i> = 1, converted to <i>prec</i> = 2</item>
        /// <item><i>prec</i> = 2, 1'</item>
        /// <item><i>prec</i> = 3, 0.1'</item>
        /// <item><i>prec</i> = 4, 0.01'</item>
        /// <item><i>prec</i> = 5, 0.001'</item>
        /// <item>…</item>
        /// <item><i>prec</i> = 11 (max), 10^−9'</item>
        /// </list>
        /// <para>
        /// If the first 3 characters of <i>georef</i> are "INV", then <i>lat</i> and <i>lon</i> are set to <see cref="double.NaN"/> and <i>prec</i> is unchanged.
        /// </para>
        /// </remarks>
        public static (double lat, double lon, int prec) Reverse(ReadOnlySpan<char> georef, bool centerp = true)
        {
            var len = georef.Length;
            if (len >= 3 && georef.StartsWith("INV".AsSpan(), StringComparison.OrdinalIgnoreCase))
            {
                return (double.NaN, double.NaN, 0);
            }

            if (len < baselen_ - 2)
                throw new GeographicException("Georef must start with at least 2 letters:" + georef.ToString());

            int prec1 = (2 + len - baselen_) / 2 - 1;
            int k;
            k = lontile_.IndexOf(georef[0]);
            if (k < 0)
                throw new GeographicException("Bad longitude tile letter in georef: " + georef.ToString());
            double lon1 = k + lonorig_ / tile_;
            k = lattile_.IndexOf(georef[1]);
            if (k < 0)
                throw new GeographicException("Bad latitude tile letter in georef: " + georef.ToString());
            double lat1 = k + latorig_ / tile_;
            double unit = 1;

            if (len > 2)
            {
                unit *= tile_;
                k = degrees_.IndexOf(georef[2]);
                if (k < 0)
                    throw new GeographicException("Bad longitude degree letter in georef: " + georef.ToString());
                lon1 = lon1 * tile_ + k;
                if (len < 4)
                    throw new GeographicException("Missing latitude degree letter in georef: " + georef.ToString());
                k = degrees_.IndexOf(georef[3]);
                if (k < 0)
                    throw new GeographicException("Bad latitude degree letter in georef: " + georef.ToString());
                lat1 = lat1 * tile_ + k;
                if (prec1 > 0)
                {
                    if (georef.Slice(baselen_).FindFirstNotOf(digits_) != -1)
                        throw new GeographicException("Non digits in trailing portion of georef: " + georef.Slice(baselen_).ToString());
                    if (len % 2 != 0)
                        throw new GeographicException("Georef must end with an even number of digits: " + georef.Slice(baselen_).ToString());
                    if (prec1 == 1)
                        throw new GeographicException("Georef needs at least 4 digits for minutes: " + georef.Slice(baselen_).ToString());
                    if (prec1 > maxprec_)
                        throw new GeographicException($"More than {2 * maxprec_} digits in georef: " + georef.Slice(baselen_).ToString());
                    for (int i = 0; i < prec1; ++i)
                    {
                        int m = i != 0 ? base_ : 6;
                        unit *= m;
                        int
                          x = digits_.IndexOf(georef[baselen_ + i]),
                          y = digits_.IndexOf(georef[baselen_ + i + prec1]);
                        if (!(i != 0 || (x < m && y < m)))
                            throw new GeographicException("Minutes terms in georef must be less than 60 "
                                                + georef.Slice(baselen_).ToString());
                        lon1 = m * lon1 + x;
                        lat1 = m * lat1 + y;
                    }
                }
            }
            if (centerp)
            {
                unit *= 2; lat1 = 2 * lat1 + 1; lon1 = 2 * lon1 + 1;
            }
            return ((tile_ * lat1) / unit, (tile_ * lon1) / unit, prec1);
        }

        /// <summary>
        /// Gets the angular resolution of a <see cref="Georef"/>.
        /// </summary>
        /// <param name="prec">the precision of the <see cref="Georef"/>.</param>
        /// <returns>the latitude-longitude resolution (degrees).</returns>
        /// <remarks>
        /// Internally, <i>prec</i> is first put in the range [−1, 11].
        /// </remarks>
        public static double Resolution(int prec)
        {
            if (prec < 1)
                return prec < 0 ? 15 : 1;
            else
            {
                // Treat prec = 1 as 2.
                prec = Max(2, Min(maxprec_, prec));
                // Need extra real because, since C++11, pow(float, int) returns double
                return 1 / (60 * Pow(base_, prec - 2));
            }
        }

        /// <summary>
        /// Gets the <see cref="Georef"/> precision required to meet a given geographic resolution.
        /// </summary>
        /// <param name="res">the minimum of resolution in latitude and longitude (degrees).</param>
        /// <returns><see cref="Georef"/> precision</returns>
        /// <remarks>
        /// The returned length is in the range [0, 11].
        /// </remarks>
        public static double Precision(double res)
        {
            res = Abs(res);
            for (int prec = 0; prec < maxprec_; ++prec)
            {
                if (prec == 1)
                    continue;
                if (Resolution(prec) <= res)
                    return prec;
            }
            return maxprec_;
        }
    }
}
