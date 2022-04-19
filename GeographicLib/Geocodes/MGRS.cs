using System;
using System.Collections.Generic;
using System.Text;

using static System.Math;
using static GeographicLib.MathEx;

namespace GeographicLib.Geocodes
{
    /// <summary>
    /// Convert between UTM/UPS and MGRS.
    /// </summary>
    /// <remarks>
    /// <see cref="MGRS"/> is defined in Chapter 3 of
    /// <list type="bullet">
    /// <item>
    /// J. W. Hager, L. L. Fry, S. S. Jacks, D. R. Hill,
    /// <a href="https://web.archive.org/web/20161214054445/http://earth-info.nga.mil/GandG/publications/tm8358.1/pdf/TM8358_1.pdf">Datums, Ellipsoids, Grids, and Grid Reference Systems</a>,
    /// Defense Mapping Agency, Technical Manual TM8358.1 (1990).
    /// </item>
    /// </list>
    /// This document has been updated by the two NGA documents
    /// <list type="bullet">
    /// <item><a href="https://earth-info.nga.mil/php/download.php?file=coord-grids">
    /// Universal Grids and Grid Reference Systems</a>, NGA.STND.0037 (2014).</item>
    /// <item><a href="https://earth-info.nga.mil/php/download.php?file=coord-utmups">
    /// The Universal Grids and the Transverse Mercator and Polar Stereographic Map Projections</a>, NGA.SIG.0012 (2014).</item>
    /// </list>
    /// This implementation has the following properties:
    /// <list type="bullet">
    /// <item>The conversions are closed, i.e., output from Forward is legal input for Reverse and vice versa.
    /// Conversion in both directions preserve the UTM/UPS selection and the UTM zone.</item>
    /// <item>Forward followed by Reverse and vice versa is approximately the identity.
    /// (This is affected in predictable ways by errors in determining the latitude band and by loss of precision in the <see cref="MGRS"/> coordinates.)</item>
    /// <item>The trailing digits produced by Forward are consistent as the precision is varied.
    /// Specifically, the digits are obtained by operating on the easting with ⌊10^6 <i>x</i>⌋ 
    /// and extracting the required digits from the resulting number (and similarly for the northing).</item>
    /// <item>All <see cref="MGRS"/> coordinates truncate to legal 100 km blocks.
    /// All <see cref="MGRS"/> coordinates with a legal 100 km block prefix are legal (even though the latitude band letter may now belong to a neighboring band).</item>
    /// <item>The range of UTM/UPS coordinates allowed for conversion to <see cref="MGRS"/> coordinates 
    /// is the maximum consistent with staying within the letter ranges of the <see cref="MGRS"/> scheme.</item>
    /// <item>All the transformations are implemented as static methods in the <see cref="MGRS"/> class.</item>
    /// </list>
    /// The <a href="http://www.nga.mil">NGA</a> software package <a href="https://earth-info.nga.mil/index.php?dir=wgs84&amp;action=wgs84#tab_geotrans">geotrans</a>
    /// also provides conversions to and from <see cref="MGRS"/>. Version 3.0 (and earlier) suffers from some drawbacks:
    /// <list type="bullet">
    /// <item>Inconsistent rules are used to determine the whether a particular <see cref="MGRS"/> coordinate is legal. A more systematic approach is taken here.</item>
    /// <item>The underlying projections are not very accurately implemented.</item>
    /// </list>
    /// </remarks>
    public static class MGRS
    {
        private const string hemispheres_ = "SN";
        private static readonly string[] utmcols_ = new[] { "ABCDEFGH", "JKLMNPQR", "STUVWXYZ" };
        private const string utmrow_ = "ABCDEFGHJKLMNPQRSTUV";
        private static readonly string[] upscols_ = new[] { "JKLPQRSTUXYZ", "ABCFGHJKLPQR", "RSTUXYZ", "ABCFGHJ" };
        private static readonly string[] upsrows_ = new[] { "ABCDEFGHJKLMNPQRSTUVWXYZ", "ABCDEFGHJKLMNP" };
        private const string latband_ = "CDEFGHJKLMNPQRSTUVWX";
        private const string upsband_ = "ABYZ";
        private const string digits_ = "0123456789";

        // Entries are [band, x, y] either side of the band boundaries.  Units for
        // x, y are t = 100km.
        private static readonly short[] tab =
            new short[] {
              0, 5,  0,   0, 9,  0,     // south edge of band 0
              0, 5,  8,   0, 9,  8,     // north edge of band 0
              1, 5,  9,   1, 9,  9,     // south edge of band 1
              1, 5, 17,   1, 9, 17,     // north edge of band 1
              2, 5, 18,   2, 9, 18,     // etc.
              2, 5, 26,   2, 9, 26,
              3, 5, 27,   3, 9, 27,
              3, 5, 35,   3, 9, 35,
              4, 5, 36,   4, 9, 36,
              4, 5, 44,   4, 9, 44,
              5, 5, 45,   5, 9, 45,
              5, 5, 53,   5, 9, 53,
              6, 5, 54,   6, 9, 54,
              6, 5, 62,   6, 9, 62,
              7, 5, 63,   7, 9, 63,
              7, 5, 70,   7, 7, 70,   7, 7, 71,   7, 9, 71, // y = 71t crosses boundary
              8, 5, 71,   8, 6, 71,   8, 6, 72,   8, 9, 72, // between bands 7 and 8.
              8, 5, 79,   8, 8, 79,   8, 8, 80,   8, 9, 80, // y = 80t crosses boundary
              9, 5, 80,   9, 7, 80,   9, 7, 81,   9, 9, 81, // between bands 8 and 9.
              9, 5, 95,   9, 9, 95,     // north edge of band 9
            };

        private const int
            base_ = 10,
            // Top-level tiles are 10^5 m = 100 km on a side
            tilelevel_ = 5,
            // Period of UTM row letters
            utmrowperiod_ = 20,
            // Row letters are shifted by 5 for even zones
            utmevenrowshift_ = 5,
            // Maximum precision is um
            maxprec_ = 5 + 6,
            // For generating digits at maxprec
            mult_ = 1000000;

        internal const int
            tile_ = 100000,            // Size MGRS blocks
            minutmcol_ = 1,
            maxutmcol_ = 9,
            minutmSrow_ = 10,
            maxutmSrow_ = 100,         // Also used for UTM S false northing
            minutmNrow_ = 0,           // Also used for UTM N false northing
            maxutmNrow_ = 95,
            minupsSind_ = 8,           // These 4 ind's apply to easting and northing
            maxupsSind_ = 32,
            minupsNind_ = 13,
            maxupsNind_ = 27,
            upseasting_ = 20,          // Also used for UPS false northing
            utmeasting_ = 5,           // UTM false easting
                                       // Difference between S hemisphere northing and N hemisphere northing
            utmNshift_ = (maxutmSrow_ - minutmNrow_) * tile_;

        private static readonly int[] mineasting_ = new[] { minupsSind_, minupsNind_, minutmcol_, minutmcol_ };
        private static readonly int[] maxeasting_ = new[] { maxupsSind_, maxupsNind_, maxutmcol_, maxutmcol_ };
        private static readonly int[] minnorthing_ = new[] { minupsSind_, minupsNind_, minutmSrow_, minutmSrow_ - (maxutmSrow_ - minutmNrow_) };
        private static readonly int[] maxnorthing_ = new[] { maxupsSind_, maxupsNind_, maxutmNrow_ + (maxutmSrow_ - minutmNrow_), maxutmNrow_ };

        // The smallest angle s.t., 90 - angeps() < 90 (approx 50e-12 arcsec)
        // 7 = ceil(log_2(90))
        private static readonly double angeps = Ldexp(1, -(DBL_MANT_DIG - 7));

        // The smallest length s.t., 1.0e7 - eps() < 1.0e7 (approx 1.9 nm)
        // 25 = ceil(log_2(2e7)) -- use half circumference here because
        // northing 195e5 is a legal in the "southern" hemisphere.
        private static readonly double eps = Ldexp(1, -(DBL_MANT_DIG - 25));

        private static void CheckCoords(bool utmp, ref bool northp, ref double x, ref double y)
        {
            // Limits are all multiples of 100km and are all closed on the lower end
            // and open on the upper end -- and this is reflected in the error
            // messages.  However if a coordinate lies on the excluded upper end (e.g.,
            // after rounding), it is shifted down by eps.  This also folds UTM
            // northings to the correct N/S hemisphere.

            int
              ix = (int)Floor(x / tile_),
              iy = (int)Floor(y / tile_),
              ind = (utmp ? 2 : 0) + (northp ? 1 : 0);
            if (!(ix >= mineasting_[ind] && ix < maxeasting_[ind]))
            {
                if (ix == maxeasting_[ind] && x == maxeasting_[ind] * tile_)
                    x -= eps;
                else
                    throw new GeographicException($"Easting {(int)Floor(x / 1000)}km not in MGRS/{(utmp ? "UTM" : "UPS")} " +
                        $"range for {(northp ? "N" : "S")} hemisphere [{mineasting_[ind] * tile_ / 1000}km, " +
                        $"{maxeasting_[ind] * tile_ / 1000}km)");
            }
            if (!(iy >= minnorthing_[ind] && iy < maxnorthing_[ind]))
            {
                if (iy == maxnorthing_[ind] && y == maxnorthing_[ind] * tile_)
                    y -= eps;
                else
                    throw new GeographicException($"Northing {(int)Floor(y / 1000)}km not in MGRS/{(utmp ? "UTM" : "UPS")} " +
                        $"range for {(northp ? "N" : "S")} hemisphere " +
                        $"[{minnorthing_[ind] * tile_ / 1000}km, {maxnorthing_[ind] * tile_ / 1000}km)");
            }

            // Correct the UTM northing and hemisphere if necessary
            if (utmp)
            {
                if (northp && iy < minutmNrow_)
                {
                    northp = false;
                    y += utmNshift_;
                }
                else if (!northp && iy >= maxutmSrow_)
                {
                    if (y == maxutmSrow_ * tile_)
                        // If on equator retain S hemisphere
                        y -= eps;
                    else
                    {
                        northp = true;
                        y -= utmNshift_;
                    }
                }
            }
        }

        private static int UTMRow(int iband, int icol, int irow)
        {
            // Input is iband = band index in [-10, 10) (as returned by LatitudeBand),
            // icol = column index in [0,8) with origin of easting = 100km, and irow =
            // periodic row index in [0,20) with origin = equator.  Output is true row
            // index in [-90, 95).  Returns maxutmSrow_ = 100, if irow and iband are
            // incompatible.

            // Estimate center row number for latitude band
            // 90 deg = 100 tiles; 1 band = 8 deg = 100*8/90 tiles
            var c = 100 * (8 * iband + 4) / (double)qd;
            bool northp = iband >= 0;
            // These are safe bounds on the rows
            //  iband minrow maxrow
            //   -10    -90    -81
            //    -9    -80    -72
            //    -8    -71    -63
            //    -7    -63    -54
            //    -6    -54    -45
            //    -5    -45    -36
            //    -4    -36    -27
            //    -3    -27    -18
            //    -2    -18     -9
            //    -1     -9     -1
            //     0      0      8
            //     1      8     17
            //     2     17     26
            //     3     26     35
            //     4     35     44
            //     5     44     53
            //     6     53     62
            //     7     62     70
            //     8     71     79
            //     9     80     94
            int
              minrow = iband > -10 ?
              (int)Floor(c - 4.3 - 0.1 * (northp ? 1 : 0)) : -90,
              maxrow = iband < 9 ?
              (int)Floor(c + 4.4 - 0.1 * (northp ? 1 : 0)) : 94,
              baserow = (minrow + maxrow) / 2 - utmrowperiod_ / 2;
            // Offset irow by the multiple of utmrowperiod_ which brings it as close as
            // possible to the center of the latitude band, (minrow + maxrow) / 2.
            // (Add maxutmSrow_ = 5 * utmrowperiod_ to ensure operand is positive.)
            irow = (irow - baserow + maxutmSrow_) % utmrowperiod_ + baserow;
            if (!(irow >= minrow && irow <= maxrow))
            {
                // Outside the safe bounds, so need to check...
                // Northing = 71e5 and 80e5 intersect band boundaries
                //   y = 71e5 in scol = 2 (x = [3e5,4e5] and x = [6e5,7e5])
                //   y = 80e5 in scol = 1 (x = [2e5,3e5] and x = [7e5,8e5])
                // This holds for all the ellipsoids given in NGA.SIG.0012_2.0.0_UTMUPS.
                // The following deals with these special cases.
                int
                  // Fold [-10,-1] -> [9,0]
                  sband = iband >= 0 ? iband : -iband - 1,
                  // Fold [-90,-1] -> [89,0]
                  srow = irow >= 0 ? irow : -irow - 1,
                  // Fold [4,7] -> [3,0]
                  scol = icol < 4 ? icol : -icol + 7;
                // For example, the safe rows for band 8 are 71 - 79.  However row 70 is
                // allowed if scol = [2,3] and row 80 is allowed if scol = [0,1].
                if (!((srow == 70 && sband == 8 && scol >= 2) ||
                         (srow == 71 && sband == 7 && scol <= 2) ||
                         (srow == 79 && sband == 9 && scol >= 1) ||
                         (srow == 80 && sband == 8 && scol <= 1)))
                    irow = maxutmSrow_;
            }
            return irow;
        }

        /// <summary>
        /// Return latitude band number [-10, 10) for the given latitude (degrees).
        /// The bands are reckoned in include their southern edges.
        /// </summary>
        /// <param name="lat"></param>
        /// <returns></returns>
        internal static int LatitudeBand(double lat)
        {
            var ilat = (int)Floor(lat);
            return Max(-10, Min(9, (ilat + 80) / 8 - 10));
        }

        /// <summary>
        /// Return approximate latitude band number [-10, 10) for the given northing
        /// (meters).  With this rule, each 100km tile would have a unique band
        /// letter corresponding to the latitude at the center of the tile.  This
        /// function isn't currently used.
        /// </summary>
        /// <param name="y"></param>
        /// <returns></returns>
        static int ApproxLatitudeBand(double y)
        {
            // northing at tile center in units of tile = 100km
            var ya = Floor(Min(88, Abs(y / tile_))) + 0.5;
            // convert to lat (mult by 90/100) and then to band (divide by 8)
            // the +1 fine tunes the boundary between bands 3 and 4
            var b = (int)Floor(((ya * 9 + 1) / 10) / 8);
            // For the northern hemisphere we have
            // band rows  num
            // N 0   0:8    9
            // P 1   9:17   9
            // Q 2  18:26   9
            // R 3  27:34   8
            // S 4  35:43   9
            // T 5  44:52   9
            // U 6  53:61   9
            // V 7  62:70   9
            // W 8  71:79   9
            // X 9  80:94  15
            return y >= 0 ? b : -(b + 1);
        }

        /// <summary>
        /// Convert UTM or UPS coordinate to an MGRS coordinate.
        /// </summary>
        /// <param name="zone">UTM zone (zero means UPS).</param>
        /// <param name="northp">hemisphere (<see langword="true"/> means north, <see langword="false"/> means south).</param>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <param name="prec">precision relative to 100 km.</param>
        /// <returns>A <see cref="MGRS"/> string.</returns>
        /// <remarks>
        /// <i>prec</i> specifies the precision of the <see cref="MGRS"/> as follows:
        /// <list type="bullet">
        /// <item><i>prec</i> = -1 (min), only the grid zone is returned</item>
        /// <item><i>prec</i> = 0, 100km</item>
        /// <item><i>prec</i> = 1, 10km</item>
        /// <item><i>prec</i> = 2, 1km</item>
        /// <item><i>prec</i> = 3, 100m</item>
        /// <item><i>prec</i> = 4, 10m</item>
        /// <item><i>prec</i> = 5, 1m</item>
        /// <item><i>prec</i> = 6, 0.1m</item>
        /// <item>...</item>
        /// <item><i>prec</i> = 11 (max), 1μm</item>
        /// </list>
        /// <para>
        /// UTM eastings are allowed to be in the range [100 km, 900 km],
        /// northings are allowed to be in in [0 km, 9500 km] for the northern hemisphere and in [1000 km, 10000 km] for the southern hemisphere.
        /// (However UTM northings can be continued across the equator. So the actual limits on the northings are [−9000 km, 9500 km]
        /// for the "northern" hemisphere and [1000 km, 19500 km] for the "southern" hemisphere.)
        /// </para>
        /// <para>
        /// UPS eastings/northings are allowed to be in the range [1300 km, 2700 km] in the northern hemisphere
        /// and in [800 km, 3200 km] in the southern hemisphere.
        /// </para>
        /// <para>
        /// The ranges are 100 km more restrictive than for the conversion between geographic coordinates and UTM and UPS given by <see cref="UTMUPS"/>.
        /// These restrictions are dictated by the allowed letters in <see cref="MGRS"/> coordinates.
        /// The choice of 9500 km for the maximum northing for northern hemisphere and of 1000 km as the minimum northing for southern hemisphere
        /// provide at least 0.5 degree extension into standard UPS zones.
        /// The upper ends of the ranges for the UPS coordinates is dictated by requiring symmetry about the meridians <c>0E</c> and <c>90E</c>.
        /// </para>
        /// <para>
        /// All allowed UTM and UPS coordinates may now be converted to legal <see cref="MGRS"/> coordinates with the proviso that eastings and northings
        /// on the upper boundaries are silently reduced by about 4 nm (4 nanometers) to place them within the allowed range.
        /// (This includes reducing a southern hemisphere northing of 10000 km by 4 nm so that it is placed in latitude band <c>M</c>.)
        /// The UTM or UPS coordinates are truncated to requested precision to determine the <see cref="MGRS"/> coordinate.
        /// Thus in UTM zone <c>38N</c>, the square area with easting in [444 km, 445 km) and northing in [3688 km, 3689 km) maps to <see cref="MGRS"/>
        /// coordinate <c>38SMB4488</c> (at <i>prec</i> = 2, 1 km), Khulani Sq., Baghdad.
        /// </para>
        /// <para>
        /// The UTM/UPS selection and the UTM zone is preserved in the conversion to <see cref="MGRS"/> coordinate.
        /// Thus for <i>zone</i> > 0, the <see cref="MGRS"/> coordinate begins with the zone number followed by one of [<c>C</c>–<c>M</c>] for the southern hemisphere
        /// and [<c>N</c>–<c>X</c>] for the northern hemisphere. For <i>zone</i> = 0, the <see cref="MGRS"/> coordinates begins with one of [<c>AB</c>] for the 
        /// southern hemisphere and [<c>XY</c>] for the northern hemisphere.
        /// </para>
        /// <para>
        /// The conversion to the <see cref="MGRS"/> is exact for prec in [0, 5] except that a neighboring 
        /// latitude band letter may be given if the point is within 5nm of a band boundary. For <i>prec</i> in [6, 11], the conversion is accurate to roundoff.
        /// </para>
        /// <para>
        /// If <i>prec</i> = −1, then the "grid zone designation", e.g., <c>18T</c>, is returned.
        /// This consists of the UTM zone number (absent for UPS) and the first letter of the <see cref="MGRS"/> string which labels the latitude band
        /// for UTM and the hemisphere for UPS.
        /// </para>
        /// <para>
        /// If <i>x</i> or <i>y</i> is <see cref="double.NaN"/> or if zone is <see cref="ZoneSpec.Invalid"/>, the returned <see cref="MGRS"/> string is "INVALID".
        /// </para>
        /// </remarks>
        public static string Forward(int zone, bool northp, double x, double y, int prec)
        {
            double lat;
            if (zone > 0)
            {
                // Does a rough estimate for latitude determine the latitude band?
                var ys = northp ? y : y - utmNshift_;
                // A cheap calculation of the latitude which results in an "allowed"
                // latitude band would be
                //   lat = ApproxLatitudeBand(ys) * 8 + 4;
                //
                // Here we do a more careful job using the band letter corresponding to
                // the actual latitude.
                ys /= tile_;
                if (Abs(ys) < 1)
                    lat = 0.9 * ys;         // accurate enough estimate near equator
                else
                {
                    double
                      // The poleward bound is a fit from above of lat(x,y)
                      // for x = 500km and y = [0km, 950km]
                      latp = 0.901 * ys + (ys > 0 ? 1 : -1) * 0.135,
                      // The equatorward bound is a fit from below of lat(x,y)
                      // for x = 900km and y = [0km, 950km]
                      late = 0.902 * ys * (1 - 1.85e-6 * ys * ys);
                    if (LatitudeBand(latp) == LatitudeBand(late))
                        lat = latp;
                    else
                        // bounds straddle a band boundary so need to compute lat accurately
                        (lat, _) = UTMUPS.Reverse(zone, northp, x, y);
                }
            }
            else
                // Latitude isn't needed for UPS specs or for INVALID
                lat = 0;

            return Forward(zone, northp, x, y, lat, prec);
        }

        /// <summary>
        /// Convert UTM or UPS coordinate to an <see cref="MGRS"/> coordinate when the latitude is known.
        /// </summary>
        /// <param name="zone">UTM zone (zero means UPS).</param>
        /// <param name="northp">hemisphere (true means north, false means south).</param>
        /// <param name="x">easting of point (meters).</param>
        /// <param name="y">northing of point (meters).</param>
        /// <param name="lat">latitude (degrees).</param>
        /// <param name="prec">precision relative to 100 km.</param>
        /// <returns>A <see cref="MGRS"/> string.</returns>
        /// <remarks>
        /// The latitude is ignored for <i>zone</i> = 0 (UPS);
        /// otherwise the latitude is used to determine the latitude band and this is checked for consistency using the same tests as <see cref="Reverse"/>.
        /// </remarks>
        public static string Forward(int zone, bool northp, double x, double y, double lat, int prec)
        {
            if (zone == (int)ZoneSpec.Invalid || double.IsNaN(x) || double.IsNaN(y) || double.IsNaN(lat))
            {
                return "INVALID";
            }

            bool utmp = zone != 0;
            CheckCoords(utmp, ref northp, ref x, ref y);

            if (!(zone >= (int)ZoneSpec.MinZone && zone <= (int)ZoneSpec.MaxZone))
                throw new GeographicException($"Zone {zone} not in [0,60]");
            if (!(prec >= -1 && prec <= maxprec_))
                throw new GeographicException($"MGRS precision {prec} not in [-1, {maxprec_} ]");

            // Fixed char array for accumulating string.  Allow space for zone, 3 block
            // letters, easting + northing.  Don't need to allow for terminating null.
            Span<char> mgrs1 = stackalloc char[2 + 3 + 2 * maxprec_];
            int
              zone1 = zone - 1,
              z = utmp ? 2 : 0,
              mlen = z + 3 + 2 * prec;

            if (utmp)
            {
                mgrs1[0] = digits_[zone / base_];
                mgrs1[1] = digits_[zone % base_];
                // This isn't necessary...!  Keep y non-neg
                // if (!northp) y -= maxutmSrow_ * tile_;
            }

            long
             ix = (long)Floor(x * mult_),
             iy = (long)Floor(y * mult_),
             m = (long)mult_ * tile_;
            int xh = (int)(ix / m), yh = (int)(iy / m);

            if (utmp)
            {
                int
                  // Correct fuzziness in latitude near equator
                  iband = Abs(lat) < angeps ? (northp ? 0 : -1) : LatitudeBand(lat),
                  icol = xh - minutmcol_,
                  irow = UTMRow(iband, icol, yh % utmrowperiod_);

                if (irow != yh - (northp ? minutmNrow_ : maxutmSrow_))
                    throw new GeographicException($"Latitude {lat} is inconsistent with UTM coordinates");

                mgrs1[z++] = latband_[10 + iband];
                mgrs1[z++] = utmcols_[zone1 % 3][icol];
                mgrs1[z++] = utmrow_[(yh + (((zone1 & 1) == 1) ? utmevenrowshift_ : 0))
                                   % utmrowperiod_];
            }
            else
            {
                bool eastp = xh >= upseasting_;
                int iband = (northp ? 2 : 0) + (eastp ? 1 : 0);
                mgrs1[z++] = upsband_[iband];
                mgrs1[z++] = upscols_[iband][xh - (eastp ? upseasting_ :
                                                   (northp ? minupsNind_ :
                                                    minupsSind_))];
                mgrs1[z++] = upsrows_[northp ? 1 : 0][yh - (northp ? minupsNind_ : minupsSind_)];
            }
            if (prec > 0)
            {
                ix -= m * xh; iy -= m * yh;
                var d = (long)Pow(base_, maxprec_ - prec);
                ix /= d; iy /= d;
                for (int c = prec; c-- > 0;)
                {
                    mgrs1[z + c] = digits_[(int)(ix % base_)]; ix /= base_;
                    mgrs1[z + c + prec] = digits_[(int)(iy % base_)]; iy /= base_;
                }
            }

            return mgrs1.Slice(0, mlen).ToString();
        }

        /// <summary>
        /// Convert a <see cref="MGRS"/> coordinate to UTM or UPS coordinates.
        /// </summary>
        /// <param name="mgrs"><see cref="MGRS"/> string.</param>
        /// <param name="centerp">if <see langword="true"/> (default), return center of the <see cref="MGRS"/> square, else return SW (lower left) corner.</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>zone</i>, UTM zone (zero means UPS).</item>
        /// <item><i>northp</i>, hemisphere (<see langword="true"/> means north, <see langword="false"/> means south).</item>
        /// <item><i>x</i>, easting of point (meters).</item>
        /// <item><i>y</i>, northing of point (meters).</item>
        /// <item><i>prec</i>, precision relative to 100 km.</item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// All conversions from <see cref="MGRS"/> to UTM/UPS are permitted provided the <see cref="MGRS"/> coordinate is a possible result
        /// of a conversion in the other direction. (The leading 0 may be dropped from an input <see cref="MGRS"/>  coordinate for UTM zones 1–9.)
        /// In addition, <see cref="MGRS"/>  coordinates with a neighboring latitude band letter are permitted provided that some portion of the 100 km block
        /// is within the given latitude band. Thus
        /// <list type="bullet">
        /// <item><c>38VLS</c> and <c>38WLS</c> are allowed (latitude <c>64N</c> intersects the square <c>38[VW]LS</c>);
        /// but <c>38VMS</c> is not permitted (all of <c>38WMS</c> is north of <c>64N</c>)</item>
        /// <item><c>38MPE</c> and <c>38NPF</c> are permitted (they straddle the equator);
        /// but <c>38NPE</c> and <c>38MPF</c> are not permitted (the equator does not intersect either block).</item>
        /// <item>Similarly <c>ZAB</c> and <c>YZB</c> are permitted (they straddle the prime meridian);
        /// but <c>YAB</c> and <c>ZZB</c> are not (the prime meridian does not intersect either block).</item>
        /// </list>
        /// <para>
        /// The UTM/UPS selection and the UTM zone is preserved in the conversion from <see cref="MGRS"/> coordinate.
        /// The conversion is exact for <i>prec</i> in [0, 5].
        /// With <i>centerp</i> = <see langword="true"/>, the conversion from <see cref="MGRS"/> to geographic and back is stable.
        /// This is not assured if <i>centerp</i> = <see langword="false"/>.
        /// </para>
        /// <para>
        /// If a "grid zone designation" (for example, <c>18T</c> or <c>A</c>) is given,
        /// then some suitable (but essentially arbitrary) point within that grid zone is returned.
        /// The main utility of the conversion is to allow zone and  <i>northp</i> to be determined.
        /// In this case, the centerp parameter is ignored and <i>prec</i> is set to <c>−1</c>.
        /// </para>
        /// <para>
        /// If the first 3 characters of <paramref name="mgrs"/> are "INV", 
        /// then <i>x</i> and <i>y</i> are set to <see cref="double.NaN"/>, <i>zone</i> is set to <see cref="ZoneSpec.Invalid"/>,
        /// and <i>prec</i> is set to <c>−2</c>.
        /// </para>
        /// </remarks>
        public static (int zone, bool northp, double x, double y, int prec) Reverse(ReadOnlySpan<char> mgrs, bool centerp = true)
        {
            int
              p = 0,
              len = mgrs.Length;

            int prec;
            double x, y;

            if (len >= 3 && mgrs.StartsWith("INV".AsSpan(), StringComparison.OrdinalIgnoreCase))
            {
                return ((int)ZoneSpec.Invalid, false, double.NaN, double.NaN, -2);
            }

            int zone1 = 0;
            while (p < len)
            {
                int i = digits_.IndexOf(mgrs[p]);
                if (i < 0)
                    break;
                zone1 = 10 * zone1 + i;
                ++p;
            }

            if (p > 0 && !(zone1 >= (int)ZoneSpec.MinUTMZone && zone1 <= (int)ZoneSpec.MaxUTMZone))
                throw new GeographicException($"Zone {zone1} not in [1,60]");
            if (p > 2)
                throw new GeographicException("More than 2 digits at start of MGRS: " + mgrs.Slice(0, p).ToString());
            if (len - p < 1)
                throw new GeographicException("MGRS string too short:" + mgrs.ToString());

            bool utmp = zone1 != (int)ZoneSpec.UPS;
            int zonem1 = zone1 - 1;
            var band = utmp ? latband_ : upsband_;
            int iband = band.IndexOf(mgrs[p++], StringComparison.InvariantCultureIgnoreCase);

            if (iband < 0)
                throw new GeographicException($"Band letter {mgrs[p - 1]} not in {(utmp ? "UTM" : "UPS")} set {band}");

            bool northp1 = iband >= (utmp ? 10 : 2);
            if (p == len)
            {             // Grid zone only (ignore centerp)
                          // Approx length of a degree of meridian arc in units of tile.
                var deg = (double)utmNshift_ / (qd * tile_);

                if (utmp)
                {
                    // Pick central meridian except for 31V
                    x = ((zone1 == 31 && iband == 17) ? 4 : 5) * tile_;
                    // Pick center of 8deg latitude bands
                    y = Floor(8 * (iband - 9.5) * deg + 0.5) * tile_ + (northp1 ? 0 : utmNshift_);
                }
                else
                {
                    // Pick point at lat 86N or 86S
                    x = ((((iband & 1) == 1) ? 1 : -1) * Floor(4 * deg + 0.5) + upseasting_) * tile_;
                    // Pick point at lon 90E or 90W.
                    y = upseasting_ * tile_;
                }
                prec = -1;

                return (zone1, northp1, x, y, prec);
            }
            else if (len - p < 2)
                throw new GeographicException("Missing row letter in: " + mgrs.ToString());

            var col = utmp ? utmcols_[zonem1 % 3] : upscols_[iband];
            var row = utmp ? utmrow_ : upsrows_[northp1 ? 1 : 0];
            int icol = col.IndexOf(mgrs[p++], StringComparison.InvariantCultureIgnoreCase);

            if (icol < 0)
                throw new GeographicException($"Column letter {mgrs[p - 1]} not in " +
                    $"{(utmp ? "zone " + mgrs.Slice(0, p - 2).ToString() : "UPS band " + mgrs[p - 2])} set {col}");

            int irow = row.IndexOf(mgrs[p++], StringComparison.InvariantCultureIgnoreCase);
            if (irow < 0)
                throw new GeographicException($"Row letter {mgrs[p - 1]} not in " +
                    $"{(utmp ? "UTM" : ("UPS " + hemispheres_[northp1 ? 1 : 0]))} set {row}");

            if (utmp)
            {
                if ((zonem1 & 1) is 1)
                    irow = (irow + utmrowperiod_ - utmevenrowshift_) % utmrowperiod_;
                iband -= 10;
                irow = UTMRow(iband, icol, irow);

                if (irow == maxutmSrow_)
                    throw new GeographicException($"Block {mgrs.Slice(p - 2, 2).ToString()} not in zone/band {mgrs.Slice(0, p - 2).ToString()}");

                irow = northp1 ? irow : irow + 100;
                icol = icol + minutmcol_;
            }
            else
            {
                bool eastp = (iband & 1) is 1;
                icol += eastp ? upseasting_ : (northp1 ? minupsNind_ : minupsSind_);
                irow += northp1 ? minupsNind_ : minupsSind_;
            }

            int prec1 = (len - p) / 2;
            double
              unit = 1,
              x1 = icol,
              y1 = irow;

            for (int i = 0; i < prec1; ++i)
            {
                unit *= base_;
                int
                  ix = digits_.IndexOf(mgrs[p + i]),
                  iy = digits_.IndexOf(mgrs[p + i + prec1]);
                if (ix < 0 || iy < 0)
                    throw new GeographicException("Encountered a non-digit in " + mgrs.Slice(p).ToString());
                x1 = base_ * x1 + ix;
                y1 = base_ * y1 + iy;
            }

            if (((len - p) % 2) != 0)
            {
                if (digits_.IndexOf(mgrs[len - 1]) < 0)
                    throw new GeographicException("Encountered a non-digit in " + mgrs.Slice(p).ToString());
                else
                    throw new GeographicException("Not an even number of digits in " + mgrs.Slice(p).ToString());
            }

            if (prec1 > maxprec_)
                throw new GeographicException($"More than {2 * maxprec_} digits in: {mgrs.Slice(p).ToString()}");

            if (centerp)
            {
                unit *= 2; x1 = 2 * x1 + 1; y1 = 2 * y1 + 1;
            }

            return (zone1, northp1, (tile_ * x1) / unit, (tile_ * y1) / unit, prec1);
        }

        /// <summary>
        /// Gets a value representing the equatorial radius of the WGS84 ellipsoid (meters).
        /// </summary>
        /// <remarks>
        /// (The WGS84 value is returned because the UTM and UPS projections are based on this ellipsoid.)
        /// </remarks>
        public static double EquatorialRadius => UTMUPS.EquatorialRadius;

        /// <summary>
        /// Gets a value representing the flattening of the WGS84 ellipsoid.
        /// </summary>
        /// <remarks>
        /// (The WGS84 value is returned because the UTM and UPS projections are based on this ellipsoid.)
        /// </remarks>
        public static double Flattening => UTMUPS.Flattening;

        /// <summary>
        /// Perform some checks on the <see cref="UTMUPS"/> coordinates on this ellipsoid.
        /// Throw an error if any of the assumptions made in the <see cref="MGRS"/> class is not true.
        /// This check needs to be carried out if the ellipsoid parameters (or the UTM/UPS scales) are ever changed.
        /// </summary>
        public static void Check()
        {
            double t = tile_;
            var (_, lon) = UTMUPS.Reverse(31, true, 1 * t, 0 * t);
            if (!(lon < 0))
                throw new GeographicException("MGRS::Check: equator coverage failure");
            var (lat, _) = UTMUPS.Reverse(31, true, 1 * t, 95 * t);
            if (!(lat > 84))
                throw new GeographicException("MGRS::Check: UTM doesn't reach latitude = 84");
            (lat, _) = UTMUPS.Reverse(31, false, 1 * t, 10 * t);
            if (!(lat < -80))
                throw new GeographicException("MGRS::Check: UTM doesn't reach latitude = -80");
            var (_, _, x, _) = UTMUPS.Forward(56, 3, 32);
            if (!(x > 1 * t))
                throw new GeographicException("MGRS::Check: Norway exception creates a gap");
            (_, _, x, _) = UTMUPS.Forward(72, 21, 35);
            if (!(x > 1 * t))
                throw new GeographicException("MGRS::Check: Svalbard exception creates a gap");
            (lat, _) = UTMUPS.Reverse(0, true, 20 * t, 13 * t);
            if (!(lat < 84))
                throw new GeographicException("MGRS::Check: North UPS doesn't reach latitude = 84");
            (lat, _) = UTMUPS.Reverse(0, false, 20 * t, 8 * t);
            if (!(lat > -80))
                throw new GeographicException("MGRS::Check: South UPS doesn't reach latitude = -80");

            var bandchecks = tab.Length / 3;
            for (int i = 0; i < bandchecks; ++i)
            {
                (lat, _) = UTMUPS.Reverse(38, true, tab[3 * i + 1] * t, tab[3 * i + 2] * t);
                if (!(LatitudeBand(lat) == tab[3 * i + 0]))
                    throw new GeographicException(
                        $"MGRS::Check: Band error, b = {tab[3 * i + 0]}, x = {tab[3 * i + 1]}00km, y = {tab[3 * i + 2]}00km");
            }
        }
    }
}
