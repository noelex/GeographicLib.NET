using System;
using System.Collections.Generic;
using System.Text;

using static System.Math;
using static GeographicLib.MathEx;

namespace GeographicLib.Geocodes
{
    /// <summary>
    /// Conversions for the Global Area Reference System (<see cref="GARS"/>).
    /// </summary>
    /// <remarks>
    /// The Global Area Reference System is described in
    /// <list type="bullet">
    /// <item><a href="https://en.wikipedia.org/wiki/Global_Area_Reference_System"></a></item>
    /// <item><a href="https://earth-info.nga.mil/index.php?dir=coordsys&amp;action=coordsys#tab_gars"></a></item>
    /// </list>
    /// <para>
    /// It provides a compact string representation of a geographic area (expressed as latitude and longitude).
    /// The classes <see cref="Georef"/> and <see cref="Geohash"/> implement similar compact representations.
    /// </para>
    /// </remarks>
    public static class GARS
    {
        private const string digits_ = "0123456789";
        private const string letters_ = "ABCDEFGHJKLMNPQRSTUVWXYZ";

        private const int
              lonorig_ = -hd,          // Origin for longitude
              latorig_ = -qd,           // Origin for latitude
              baselon_ = 10,            // Base for longitude tiles
              baselat_ = 24,            // Base for latitude tiles
              lonlen_ = 3,
              latlen_ = 2,
              baselen_ = lonlen_ + latlen_,
              mult1_ = 2,               // base precision = 1/2 degree
              mult2_ = 2,               // 6th char gives 2x more precision
              mult3_ = 3,               // 7th char gives 3x more precision
              m_ = mult1_ * mult2_ * mult3_,
              maxprec_ = 2,
              maxlen_ = baselen_ + maxprec_;

        /// <summary>
        /// Convert from geographic coordinates to <see cref="GARS"/>.
        /// </summary>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <param name="prec">the precision of the resulting <see cref="GARS"/>.</param>
        /// <returns>the <see cref="GARS"/> string.</returns>
        /// <remarks>
        /// <i>prec</i> specifies the precision of gars as follows:
        /// <list type="bullet">
        /// <item><i>prec</i> = 0 (min), 30' precision, e.g., 006AG;</item>
        /// <item><i>prec</i> = 1, 15' precision, e.g., 006AG3;</item>
        /// <item><i>prec</i> = 2 (max), 5' precision, e.g., 006AG39.</item>
        /// </list>
        /// <para>
        /// If <i>lat</i> or <i>lon</i> is <see cref="double.NaN"/>, then "INVALID" is returned.
        /// </para>
        /// </remarks>
        public static string Forward(double lat, double lon, int prec)
        {
            if (Abs(lat) > qd)
                throw new GeographicException($"Latitude {lat}d not in [-{qd}d, {qd}d]");
            if (double.IsNaN(lat) || double.IsNaN(lon))
            {
                return "INVALID";
            }

            lon = AngNormalize(lon);
            if (lon == hd) lon = -hd; // lon now in [-180,180)
            if (lat == qd) lat *= (1 - DBL_EPSILON / 2);

            prec = Max(0, Min(maxprec_, prec));
            int
              x = (int)Floor(lon * m_) - lonorig_ * m_,
              y = (int)Floor(lat * m_) - latorig_ * m_,
              ilon = x * mult1_ / m_,
              ilat = y * mult1_ / m_;
            x -= ilon * m_ / mult1_; y -= ilat * m_ / mult1_;
            Span<char> gars1=stackalloc char[maxlen_];
            ++ilon;
            for (int c = lonlen_; c--!=0;)
            {
                gars1[c] = digits_[ilon % baselon_]; ilon /= baselon_;
            }
            for (int c = latlen_; c-- != 0;)
            {
                gars1[lonlen_ + c] = letters_[ilat % baselat_]; ilat /= baselat_;
            }
            if (prec > 0)
            {
                ilon = x / mult3_; ilat = y / mult3_;
                gars1[baselen_] = digits_[mult2_ * (mult2_ - 1 - ilat) + ilon + 1];
                if (prec > 1)
                {
                    ilon = x % mult3_; ilat = y % mult3_;
                    gars1[baselen_ + 1] = digits_[mult3_ * (mult3_ - 1 - ilat) + ilon + 1];
                }
            }

            return gars1.Slice(0, baselen_ + prec).ToString();
        }

        /// <summary>
        /// Convert from <see cref="GARS"/> to geographic coordinates.
        /// </summary>
        /// <param name="gars">a <see cref="GARS"/> string.</param>
        /// <param name="centerp">if <see langword="true"/> (the default) return the center of the gars, otherwise return the south-west corner.</param>
        /// <returns>
        /// <i>lat</i>, latitude of point (degrees), <i>lon</i>, longitude of point (degrees) and <i>prec</i>, 
        /// precision of the input <see cref="GARS"/> string.
        /// </returns>
        /// <remarks>
        /// The case of the letters in gars is ignored. <i>prec</i> is in the range [0, 2] and gives the precision of gars as follows:
        /// <list type="bullet">
        /// <item><i>prec</i> = 0 (min), 30' precision, e.g., 006AG;</item>
        /// <item><i>prec</i> = 1, 15' precision, e.g., 006AG3;</item>
        /// <item><i>prec</i> = 2 (max), 5' precision, e.g., 006AG39.</item>
        /// </list>
        /// <para>
        /// If the first 3 characters of <i>gars</i> are "INV", then <i>lat</i> and <i>lon</i> are set to NaN and <i>prec</i> is unchanged.
        /// </para>
        /// </remarks>
        public static (double lat, double lon, double prec) Reverse(ReadOnlySpan<char> gars, bool centerp = true)
        {
            var len = gars.Length;
            if (len >= 3 && gars.StartsWith("INV".AsSpan(), StringComparison.OrdinalIgnoreCase))
            {
                return (double.NaN, double.NaN, 0);
            }

            if (len < baselen_)
                throw new GeographicException("GARS must have at least 5 characters: " + gars.ToString());
            if (len > maxlen_)
                throw new GeographicException("GARS can have at most 7 characters: " + gars.ToString());

            int prec1 = len - baselen_;
            int ilon = 0;
            for (int c = 0; c < lonlen_; ++c)
            {
                int k = digits_.IndexOf(gars[c]);
                if (k < 0)
                    throw new GeographicException("GARS must start with 3 digits: " + gars.ToString());
                ilon = ilon * baselon_ + k;
            }

            if (!(ilon >= 1 && ilon <= 2 * td))
                throw new GeographicException("Initial digits in GARS must lie in [1, 720]: " + gars.ToString());

            --ilon;
            int ilat = 0;
            for (int c = 0; c < latlen_; ++c)
            {
                int k = letters_.IndexOf(gars[lonlen_ + c]);
                if (k < 0)
                    throw new GeographicException("Illegal letters in GARS " + gars.Slice(3, 2).ToString());
                ilat = ilat * baselat_ + k;
            }
            if (!(ilat < td))
                throw new GeographicException("GARS letters must lie in [AA, QZ] " + gars.ToString());
            double
              unit = mult1_,
              lat1 = ilat + latorig_ * unit,
              lon1 = ilon + lonorig_ * unit;
            if (prec1 > 0)
            {
                int k = digits_.IndexOf(digits_, gars[baselen_]);
                if (!(k >= 1 && k <= mult2_ * mult2_))
                    throw new GeographicException("6th character in GARS must [1, 4] " + gars.ToString());
                --k;
                unit *= mult2_;
                lat1 = mult2_ * lat1 + (mult2_ - 1 - k / mult2_);
                lon1 = mult2_ * lon1 + (k % mult2_);
                if (prec1 > 1)
                {
                    k = digits_.IndexOf(digits_, gars[baselen_ + 1]);
                    if (!(k >= 1 /* && k <= mult3_ * mult3_ */))
                        throw new GeographicException("7th character in GARS must [1, 9] " + gars.ToString());
                    --k;
                    unit *= mult3_;
                    lat1 = mult3_ * lat1 + (mult3_ - 1 - k / mult3_);
                    lon1 = mult3_ * lon1 + (k % mult3_);
                }
            }
            if (centerp)
            {
                unit *= 2; lat1 = 2 * lat1 + 1; lon1 = 2 * lon1 + 1;
            }

            return (lat1 / unit, lon1 / unit, prec1);
        }

        /// <summary>
        /// Gets the angular resolution of a <see cref="GARS"/>.
        /// </summary>
        /// <param name="prec">the precision of the <see cref="GARS"/></param>
        /// <returns>the latitude-longitude resolution (degrees).</returns>
        /// <remarks>
        /// Internally, <i>prec</i> is first put in the range [0, 2].
        /// </remarks>
        public static double Resolution(int prec)=> 
            1 / (prec <= 0 ? mult1_ : (prec == 1 ? mult1_ * mult2_ : mult1_ * mult2_ * mult3_));

        /// <summary>
        /// Gets the <a href="https://geographiclib.sourceforge.io/html/classGeographicLib_1_1GARS.html">GARS</a> precision required to meet a given geographic resolution.
        /// </summary>
        /// <param name="res">the minimum of resolution in latitude and longitude (degrees).</param>
        /// <returns><a href="https://geographiclib.sourceforge.io/html/classGeographicLib_1_1GARS.html">GARS</a> precision.</returns>
        /// <remarks>
        /// The returned length is in the range [0, 2].
        /// </remarks>
        public static double Precision(double res)
        {
            res = Abs(res);
            for (int prec = 0; prec < maxprec_; ++prec)
                if (Resolution(prec) <= res)
                    return prec;
            return maxprec_;
        }
    }
}
