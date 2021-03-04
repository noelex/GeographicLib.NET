using GeographicLib.Projections;
using System;
using System.Collections.Generic;
using System.Text;

using static System.Math;

namespace GeographicLib.Geocodes
{
    /// <summary>
    /// Ordnance Survey grid system for Great Britain.
    /// </summary>
    /// <remarks>
    /// The class implements the coordinate system used by the Ordnance Survey for maps of Great Britain and conversions to the grid reference system.
    /// <para>
    /// <list type="bullet">
    /// <item><a href="http://www.ordnancesurvey.co.uk/docs/support/guide-coordinate-systems-great-britain.pdf">
    /// A guide to coordinate systems in Great Britain</a></item>
    /// <item><a href="http://www.ordnancesurvey.co.uk/docs/support/national-grid.pdf">Guide to the National Grid</a></item>
    /// </list>
    /// </para>
    /// <para>
    /// WARNING:
    /// </para>
    /// <para>
    /// The latitudes and longitudes for the Ordnance Survey grid system do not use the WGS84 datum.
    /// Do not use the values returned by this class in the <see cref="UTMUPS"/>, <see cref="MGRS"/>, or <see cref="Geoid"/> classes
    /// without first converting the datum (and vice versa).
    /// </para>
    /// </remarks>
    public static class OSGB
    {
        private const string letters_ = "ABCDEFGHJKLMNOPQRSTUVWXYZ";
        private const string digits_ = "0123456789";

        private const int
            base_ = 10,
            tile_ = 100000,
            tilelevel_ = 5,
            tilegrid_ = 5,
            tileoffx_ = 2 * tilegrid_,
            tileoffy_ = 1 * tilegrid_,
            minx_ = -tileoffx_ * tile_,
            miny_ = -tileoffy_ * tile_,
            maxx_ = (tilegrid_ * tilegrid_ - tileoffx_) * tile_,
            maxy_ = (tilegrid_ * tilegrid_ - tileoffy_) * tile_,
            // Maximum precision is um
            maxprec_ = 5 + 6;

        private static void CheckCoords(double x, double y)
        {
            // Limits are all multiples of 100km and are all closed on the lower end
            // and open on the upper end -- and this is reflected in the error
            // messages.  NaNs are let through.
            if (x < minx_ || x >= maxx_)
                throw new GeographicException($"Easting {(int)Floor(x / 1000)} km not in OSGB range [{ minx_ / 1000 } km, {maxx_ / 1000} km)");
            if (y < miny_ || y >= maxy_)
                throw new GeographicException($"Northing {(int)Floor(x / 1000)} km not in OSGB range [{ miny_ / 1000 } km, {maxy_ / 1000} km)");
        }

        /// <summary>
        /// Forward without returning the convergence and scale.
        /// </summary>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <returns></returns>
        public static (double x, double y) Forward(double lat, double lon) => Forward(lat, lon, out _, out _);

        /// <summary>
        /// Forward projection, from geographic to <see cref="OSGB"/> coordinates.
        /// </summary>
        /// <param name="lat"></param>
        /// <param name="lon"></param>
        /// <param name="gamma"></param>
        /// <param name="k"></param>
        /// <returns></returns>
        public static (double x, double y) Forward(double lat, double lon, out double gamma, out double k)
        {
            var (x, y) = OSGBTM.Forward(OriginLongitude, lat, lon, out gamma, out k);
            x += FalseEasting;
            y += NorthOffset;

            return (x, y);
        }

        /// <summary>
        /// Reverse without returning the convergence and scale.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static (double lat, double lon) Reverse(double x, double y) => Reverse(x, y, out _, out _);

        /// <summary>
        /// Reverse projection, from <see cref="OSGB"/> coordinates to geographic.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="gamma"></param>
        /// <param name="k"></param>
        /// <returns></returns>
        public static (double lat, double lon) Reverse(double x, double y, out double gamma, out double k)
        {
            x -= FalseEasting;
            y -= NorthOffset;
            var (lat, lon) = OSGBTM.Reverse(OriginLongitude, x, y, out gamma, out k);

            return (lat, lon);
        }

        /// <summary>
        /// Convert <see cref="OSGB"/> coordinates to a grid reference.
        /// </summary>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <param name="prec">precision relative to 100 km.</param>
        /// <returns>National Grid reference string.</returns>
        /// <remarks>
        /// <i>prec</i> specifies the precision of the grid reference string as follows:
        /// <list type="bullet">
        /// <item><i>prec</i> = 0 (min), 100km</item>
        /// <item><i>prec</i> = 1, 10km</item>
        /// <item><i>prec</i> = 2, 1km</item>
        /// <item><i>prec</i> = 3, 100m</item>
        /// <item><i>prec</i> = 4, 10m</item>
        /// <item><i>prec</i> = 5, 1m</item>
        /// <item><i>prec</i> = 6, 0.1m</item>
        /// <item><i>prec</i> = 11 (max), 1μm</item>
        /// </list>
        /// <para>
        /// The easting must be in the range [−1000 km, 1500 km) and the northing must be in the range [−500 km, 2000 km).
        /// These bounds are consistent with rules for the letter designations for the grid system.
        /// </para>
        /// <para>
        /// If <i>x</i> or <i>y</i> is <see cref="double.NaN"/>, the returned grid reference is "INVALID".
        /// </para>
        /// </remarks>
        public static string ToGridReference(double x, double y, int prec)
        {
            CheckCoords(x, y);

            if (!(prec >= 0 && prec <= maxprec_))
                throw new GeographicException($"OSGB precision {prec} not in [0, {maxprec_}]");

            if (double.IsNaN(x) || double.IsNaN(y))
            {
                return "INVALID";
            }

            Span<char> grid = stackalloc char[2 + 2 * maxprec_];

            int
              xh = (int)Floor(x / tile_),
              yh = (int)Floor(y / tile_);
            double
              xf = x - tile_ * xh,
              yf = y - tile_ * yh;

            xh += tileoffx_;
            yh += tileoffy_;

            int z = 0;

            grid[z++] = letters_[(tilegrid_ - (yh / tilegrid_) - 1)
                                * tilegrid_ + (xh / tilegrid_)];
            grid[z++] = letters_[(tilegrid_ - (yh % tilegrid_) - 1)
                                * tilegrid_ + (xh % tilegrid_)];

            // Need extra real because, since C++11, pow(float, int) returns double
            var mult = Pow(base_, Max(tilelevel_ - prec, 0));
            int
              ix = (int)Floor(xf / mult),
              iy = (int)Floor(yf / mult);
            for (int c = Min(prec, tilelevel_); c-- > 0;)
            {
                grid[z + c] = digits_[ix % base_];
                ix /= base_;
                grid[z + c + prec] = digits_[iy % base_];
                iy /= base_;
            }
            if (prec > tilelevel_)
            {
                xf -= Floor(xf / mult);
                yf -= Floor(yf / mult);
                mult = Pow(base_, prec - tilelevel_);
                ix = (int)Floor(xf * mult);
                iy = (int)Floor(yf * mult);
                for (int c = prec - tilelevel_; c-- > 0;)
                {
                    grid[z + c + tilelevel_] = digits_[ix % base_];
                    ix /= base_;
                    grid[z + c + tilelevel_ + prec] = digits_[iy % base_];
                    iy /= base_;
                }
            }
            int mlen = z + 2 * prec;

            return grid.Slice(0, mlen).ToString();
        }

        /// <summary>
        /// Convert grid reference to a <see cref="OSGB"/> coordinate.
        /// </summary>
        /// <param name="gridref">National Grid reference.</param>
        /// <param name="centerp">if <see langword="true"/> (default), return center of the grid square, else return SW (lower left) corner.</param>
        /// <returns>
        /// <i>x</i>, easting of point (meters), <i>y</i>, northing of point (meters) and <i>prec</i>, precision relative to 100 km.
        /// </returns>
        public static (double x, double y, int prec) FromGridReference(ReadOnlySpan<char> gridref, bool centerp = true)
        {
            int
              len = gridref.Length,
              p = 0;

            if (len >= 2 && gridref.StartsWith("IN".AsSpan(), StringComparison.OrdinalIgnoreCase))
            {
                return (double.NaN, double.NaN, -2); // For compatibility with MGRS::Reverse.
            }

            Span<char> grid = stackalloc char[2 + 2 * maxprec_];

            for (int i = 0; i < len; ++i)
            {
                if (!char.IsWhiteSpace(gridref[i]))
                {
                    if (p >= 2 + 2 * maxprec_)
                        throw new GeographicException("OSGB string " + gridref.ToString() + " too long");
                    grid[p++] = gridref[i];
                }
            }

            len = p;
            p = 0;
            if (len < 2)
                throw new GeographicException("OSGB string " + gridref.ToString() + " too short");
            if (len % 2 != 0)
                throw new GeographicException("OSGB string " + gridref.ToString() + " has odd number of characters");

            int
              xh = 0,
              yh = 0;

            while (p < 2)
            {
                var i = letters_.IndexOf(grid[p++]);
                if (i < 0)
                    throw new GeographicException("Illegal prefix character " + gridref.ToString());
                yh = yh * tilegrid_ + tilegrid_ - (i / tilegrid_) - 1;
                xh = xh * tilegrid_ + (i % tilegrid_);
            }

            xh -= tileoffx_;
            yh -= tileoffy_;

            int prec1 = (len - p) / 2;

            double
              unit = tile_,
              x1 = unit * xh,
              y1 = unit * yh;

            for (int i = 0; i < prec1; ++i)
            {
                unit /= base_;
                int
                  ix = digits_.IndexOf(grid[p + i]),
                  iy = digits_.IndexOf(grid[p + i + prec1]);
                if (ix < 0 || iy < 0)
                    throw new GeographicException("Encountered a non-digit in " + gridref.ToString());
                x1 += unit * ix;
                y1 += unit * iy;
            }

            if (centerp)
            {
                x1 += unit / 2;
                y1 += unit / 2;
            }

            return (x1, y1, prec1);
        }

        /// <summary>
        /// Gets a value representing the equatorial radius of the Airy 1830 ellipsoid (meters).
        /// </summary>
        /// <remarks>
        /// This is <c>20923713</c> ft converted to meters using the rule <c>1</c> ft = <c>10^(9.48401603−10)</c> m.
        /// The Airy 1830 value is returned because the <see cref="OSGB"/> projection is based on this ellipsoid.
        /// The conversion factor from feet to meters is the one used for the 1936 retriangulation of Britain;
        /// see Section A.1 (p. 37) of <i>A guide to coordinate systems in Great Britain</i>, v2.2 (Dec. 2013).
        /// </remarks>
        public static double EquatorialRadius { get; } = Pow(10, (48401603d - 100000000) / 100000000) * 20923713;

        /// <summary>
        /// Gets a value representing the inverse flattening of the Airy 1830 ellipsoid.
        /// </summary>
        /// <remarks>
        /// For the Airy 1830 ellipsoid, <i>a</i> = <c>20923713</c> ft and <i>b</i> = <c>20853810</c> ft;
        /// thus the flattening = <c>(20923713 − 20853810)/20923713</c> = <c>7767/2324857</c> = <c>1/299.32496459...</c>
        /// (The Airy 1830 value is returned because the <see cref="OSGB"/> projection is based on this ellipsoid.)
        /// </remarks>
        public static double Flattenting { get; } = (20923713d - 20853810) / 20923713;

        /// <summary>
        /// Gets a value representing the central scale for the <see cref="OSGB"/> projection (<c>0.9996012717...</c>).
        /// </summary>
        /// <remarks>
        /// C. J. Mugnier, Grids &amp; Datums, PE&amp;RS, Oct. 2003, states that this is defined as <c>10^(9.9998268−10)</c>.
        /// </remarks>
        public static double CentralScale { get; } = Pow(10, (9998268d - 10000000) / 10000000);

        /// <summary>
        /// Gets value representing latitude of the origin for the <see cref="OSGB"/> projection (<c>49</c> degrees).
        /// </summary>
        public static double OriginLatitude => 49;

        /// <summary>
        /// Gets value representing longitude of the origin for the <see cref="OSGB"/> projection (<c>-2</c> degrees).
        /// </summary>
        public static double OriginLongitude => -2;

        /// <summary>
        /// Gets value representing false northing of the <see cref="OSGB"/> projection (<c>-100000</c> meters).
        /// </summary>
        public static double FalseNorthing => -100000;

        /// <summary>
        /// Gets value representing false easting of the <see cref="OSGB"/> projection (<c>400000</c> meters).
        /// </summary>
        public static double FalseEasting => 400000;

        private static TransverseMercator OSGBTM { get; } = new TransverseMercator(EquatorialRadius, Flattenting, CentralScale);

        private static double NorthOffset { get; } = FalseNorthing - OSGBTM.Forward(0, OriginLatitude, 0).y;
    }
}
