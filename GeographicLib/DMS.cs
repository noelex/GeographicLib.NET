using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using static System.Math;
using static GeographicLib.MathEx;

namespace GeographicLib
{
    /// <summary>
    /// Convert between degrees and the DMS representation.
    /// </summary>
    /// <remarks>
    /// Parse a string representing degree, minutes, and seconds and return the angle in degrees and format an angle in degrees as degree,
    /// minutes, and seconds.In addition, handle NANs and infinities on input and output.
    /// </remarks>
    public static class DMS
    {
        private const string hemispheres_ = "SNWE";
        private const string signs_ = "-+";
        private const string digits_ = "0123456789";
        private const string dmsindicators_ = "D'\":";

        private static readonly string[] components_ = new[] { "degrees", "minutes", "seconds" };

        private static (double angle, HemisphereIndicator indicator) InternalDecode(string dmsa)
        {
            string errormsg = null;
            var ind = HemisphereIndicator.None;
            do
            {                       // Executed once (provides the ability to break)
                int sign = 1;
                int
                  beg = 0,
                  end = dmsa.Length;
                var ind1 = HemisphereIndicator.None;
                int k = -1;
                if (end > beg && (k = hemispheres_.IndexOf(dmsa[beg], StringComparison.OrdinalIgnoreCase)) >= 0)
                {
                    ind1 = (k / 2) != 0 ? HemisphereIndicator.Longitude : HemisphereIndicator.Latitude;
                    sign = k % 2 != 0 ? 1 : -1;
                    ++beg;
                }
                if (end > beg && (k = hemispheres_.IndexOf(dmsa[end - 1], StringComparison.OrdinalIgnoreCase)) >= 0)
                {
                    if (k >= 0)
                    {
                        if (ind1 != HemisphereIndicator.None)
                        {
                            if (char.ToUpper(dmsa[beg - 1]) == char.ToUpper(dmsa[end - 1]))
                                errormsg = $"Repeated hemisphere indicators {dmsa[beg - 1]} in {dmsa.Substring(beg - 1, end - beg + 1)}";
                            else
                                errormsg = $"Contradictory hemisphere indicators " +
                                    $"{dmsa[beg - 1]} and {dmsa[end - 1]} in {dmsa.Substring(beg - 1, end - beg + 1)}";
                            break;
                        }
                        ind1 = (k / 2) != 0 ? HemisphereIndicator.Longitude : HemisphereIndicator.Latitude;
                        sign = k % 2 != 0 ? 1 : -1;
                        --end;
                    }
                }
                if (end > beg && (k = signs_.IndexOf(dmsa[beg], StringComparison.OrdinalIgnoreCase)) >= 0)
                {
                    if (k >= 0)
                    {
                        sign *= k != 0 ? 1 : -1;
                        ++beg;
                    }
                }
                if (end == beg)
                {
                    errormsg = "Empty or incomplete DMS string " + dmsa;
                    break;
                }
                Span<double> ipieces = stackalloc[] { 0d, 0, 0 };
                Span<double> fpieces = stackalloc[] { 0d, 0, 0 };
                int npiece = 0;
                double icurrent = 0;
                double fcurrent = 0;
                int ncurrent = 0, p = beg;
                bool pointseen = false;
                int digcount = 0, intcount = 0;
                while (p < end)
                {
                    char x = dmsa[p++];
                    if ((k = digits_.IndexOf(x, StringComparison.OrdinalIgnoreCase)) >= 0)
                    {
                        ++ncurrent;
                        if (digcount > 0)
                            ++digcount;         // Count of decimal digits
                        else
                        {
                            icurrent = 10 * icurrent + k;
                            ++intcount;
                        }
                    }
                    else if (x == '.')
                    {
                        if (pointseen)
                        {
                            errormsg = $"Multiple decimal points in {dmsa.Substring(beg, end - beg)}";
                            break;
                        }
                        pointseen = true;
                        digcount = 1;
                    }
                    else if ((k = dmsindicators_.IndexOf(x, StringComparison.OrdinalIgnoreCase)) >= 0)
                    {
                        if (k >= 3)
                        {
                            if (p == end)
                            {
                                errormsg = "Illegal for : to appear at the end of " +
                                  dmsa.Substring(beg, end - beg);
                                break;
                            }
                            k = npiece;
                        }
                        if (k == npiece - 1)
                        {
                            errormsg = $"Repeated {components_[k]} component in {dmsa.Substring(beg, end - beg)}";
                            break;
                        }
                        else if (k < npiece)
                        {
                            errormsg = $"{components_[k]} component follows {components_[npiece - 1]} component in {dmsa.Substring(beg, end - beg)}";
                            break;
                        }
                        if (ncurrent == 0)
                        {
                            errormsg = $"Missing numbers in {components_[k]} component of {dmsa.Substring(beg, end - beg)}";
                            break;
                        }
                        if (digcount > 0)
                        {
                            fcurrent = double.Parse(dmsa.Substring(p - intcount - digcount - 1,
                                        intcount + digcount));
                            icurrent = 0;
                        }
                        ipieces[k] = icurrent;
                        fpieces[k] = icurrent + fcurrent;
                        if (p < end)
                        {
                            npiece = k + 1;
                            icurrent = fcurrent = 0;
                            ncurrent = digcount = intcount = 0;
                        }
                    }
                    else if (signs_.IndexOf(x, StringComparison.OrdinalIgnoreCase) >= 0)
                    {
                        errormsg = $"Internal sign in DMS string {dmsa.Substring(beg, end - beg)}";
                        break;
                    }
                    else
                    {
                        errormsg = $"Illegal character {x} in DMS string {dmsa.Substring(beg, end - beg)}";
                        break;
                    }
                }
                if (!string.IsNullOrEmpty(errormsg))
                    break;
                if (dmsindicators_.IndexOf(dmsa[p - 1], StringComparison.OrdinalIgnoreCase) < 0)
                {
                    if (npiece >= 3)
                    {
                        errormsg = $"Extra text following seconds in DMS string {dmsa.Substring(beg, end - beg)}";
                        break;
                    }
                    if (ncurrent == 0)
                    {
                        errormsg = $"Missing numbers in trailing component of {dmsa.Substring(beg, end - beg)}";
                        break;
                    }
                    if (digcount > 0)
                    {
                        fcurrent = double.Parse(dmsa.Substring(p - intcount - digcount,
                                      intcount + digcount));
                        icurrent = 0;
                    }
                    ipieces[npiece] = icurrent;
                    fpieces[npiece] = icurrent + fcurrent;
                }
                if (pointseen && digcount == 0)
                {
                    errormsg = "Decimal point in non-terminal component of {dmsa.Substring(beg, end - beg)}";
                    break;
                }
                // Note that we accept 59.999999... even though it rounds to 60.
                if (ipieces[1] >= 60 || fpieces[1] > 60)
                {
                    errormsg = $"Minutes {fpieces[1]} not in range [0, 60)";
                    break;
                }
                if (ipieces[2] >= 60 || fpieces[2] > 60)
                {
                    errormsg = $"Seconds {fpieces[2]} not in range [0, 60)";
                    break;
                }
                // Assume check on range of result is made by calling routine (which
                // might be able to offer a better diagnostic).
                ind = ind1;
                return (sign *
                  (fpieces[2] != 0 ?
                    (60 * (60 * fpieces[0] + fpieces[1]) + fpieces[2]) / 3600 :
                    (fpieces[1] != 0 ?
                      (60 * fpieces[0] + fpieces[1]) / 60 : fpieces[0])), ind);
            } while (false);
            var val = dmsa.NumMatch();
            if (val == 0)
                throw new GeographicException(errormsg);
            else
                ind = HemisphereIndicator.None;
            return (val, ind);
        }

        /// <summary>
        /// Convert a string in <see cref="DMS"/> to an angle.
        /// </summary>
        /// <param name="dms">string input</param>
        /// <returns>
        /// Angle in degrees and a <see cref="HemisphereIndicator"/> value indicating the presence of a hemisphere indicator.
        /// </returns>
        /// <remarks>
        /// Degrees, minutes, and seconds are indicated by the characters d, ' (single quote), " (double quote), and these components may only be
        /// given in this order. Any (but not all) components may be omitted and other symbols (e.g., the ° symbol for degrees and the unicode prime
        /// and double prime symbols for minutes and seconds) may be substituted; two single quotes can be used instead of ". The last component 
        /// indicator may be omitted and is assumed to be the next smallest unit (thus <c>33d10</c> is interpreted as <c>33d10'</c>).
        /// The final component may be a decimal fraction but the non-final components must be integers.
        /// Instead of using d, ', and " to indicate degrees, minutes, and seconds, : (colon) may be used to separate these components (numbers must
        /// appear before and after each colon); thus <c>50d30'10.3"</c> may be written as <c>50:30:10.3</c>, <c>5.5'</c> may be written <c>0:5.5</c>,
        /// and so on. The integer parts of the minutes and seconds components must be less than 60. A single leading sign is permitted.
        /// A hemisphere designator (N, E, W, S) may be added to the beginning or end of the string.
        /// The result is multiplied by the implied sign of the hemisphere designator (negative for S and W).
        /// In addition ind is set to <see cref="HemisphereIndicator.Latitude"/> if N or S is present, to <see cref="HemisphereIndicator.Longitude"/>
        /// if E or W is present, and to <see cref="HemisphereIndicator.None"/> otherwise. Throws an error on a malformed string.
        /// No check is performed on the range of the result. Examples of legal and illegal strings are
        /// <list type="bullet">
        /// <item>
        /// <i>LEGAL</i> (all the entries on each line are equivalent)
        /// <list type="bullet">
        /// <item><c>-20.51125, 20d30'40.5"S, -20°30'40.5, -20d30.675, N-20d30'40.5", -20:30:40.5</c></item>
        /// <item><c>4d0'9, 4d9", 4d9'', 4:0:9, 004:00:09, 4.0025, 4.0025d, 4d0.15, 04:.15</c></item>
        /// <item><c>4:59.99999999999999, 4:60.0, 4:59:59.9999999999999, 4:59:60.0, 5</c></item>
        /// <item>
        /// <i>ILLEGAL</i> (the exception thrown explains the problem)
        /// <list type="bullet">
        /// <item><c>4d5"4', 4::5, 4:5:, :4:5, 4d4.5'4", -N20.5, 1.8e2d, 4:60, 4:59:60</c></item>
        /// </list>
        /// </item>
        /// </list>
        /// </item>
        /// </list>
        /// The decoding operation can also perform addition and subtraction operations.
        /// If the string includes internal signs (i.e., not at the beginning nor immediately after an initial hemisphere designator),
        /// then the string is split immediately before such signs and each piece is decoded according to the above rules and the results added;
        /// thus <c>S3-2.5+4.1N</c> is parsed as the sum of <c>S3</c>, <c>-2.5</c>, <c>+4.1N</c>. Any piece can include a hemisphere designator;
        /// however, if multiple designators are given, they must compatible; e.g., you cannot mix N and E. In addition, the designator can appear
        /// at the beginning or end of the first piece, but must be at the end of all subsequent pieces (a hemisphere designator is not allowed after
        /// the initial sign). Examples of legal and illegal combinations are
        /// <list type="bullet">
        /// <item>
        /// <i>LEGAL</i> (these are all equivalent): <c>070:00:45, 70:01:15W+0:0.5, 70:01:15W-0:0:30W, W70:01:15+0:0:30E</c>
        /// </item>
        /// <item>
        /// <i>ILLEGAL</i> (the exception thrown explains the problem): <c>70:01:15W+0:0:15N, W70:01:15+W0:0:15</c>
        /// </item>
        /// </list>
        /// <para>
        /// <b>WARNING</b>: The "exponential" notation is not recognized. Thus 7.0E1 is illegal, while 7.0E+1 is parsed as (7.0E) + (+1), yielding the same result as 8.0E.
        /// </para>
        /// </remarks>
        public static (double angle, HemisphereIndicator ind) Decode(string dms)
        {
            // Here's a table of the allowed characters

            // S unicode   dec  UTF-8      descripton

            // DEGREE
            // d U+0064    100  64         d
            // D U+0044     68  44         D
            // ° U+00b0    176  c2 b0      degree symbol
            // º U+00ba    186  c2 ba      alt symbol
            // ⁰ U+2070   8304  e2 81 b0   sup zero
            // ˚ U+02da    730  cb 9a      ring above
            // ∘ U+2218   8728  e2 88 98   compose function
            // * U+002a     42  2a         GRiD symbol for degrees

            // MINUTES
            // ' U+0027     39  27         apostrophe
            // ` U+0060     96  60         grave accent
            // ′ U+2032   8242  e2 80 b2   prime
            // ‵ U+2035   8245  e2 80 b5   back prime
            // ´ U+00b4    180  c2 b4      acute accent
            // ‘ U+2018   8216  e2 80 98   left single quote (also ext ASCII 0x91)
            // ’ U+2019   8217  e2 80 99   right single quote (also ext ASCII 0x92)
            // ‛ U+201b   8219  e2 80 9b   reversed-9 single quote
            // ʹ U+02b9    697  ca b9      modifier letter prime
            // ˊ U+02ca    714  cb 8a      modifier letter acute accent
            // ˋ U+02cb    715  cb 8b      modifier letter grave accent

            // SECONDS
            // " U+0022     34  22         quotation mark
            // ″ U+2033   8243  e2 80 b3   double prime
            // ‶ U+2036   8246  e2 80 b6   reversed double prime
            // ˝ U+02dd    733  cb 9d      double acute accent
            // “ U+201c   8220  e2 80 9c   left double quote (also ext ASCII 0x93)
            // ” U+201d   8221  e2 80 9d   right double quote (also ext ASCII 0x94)
            // ‟ U+201f   8223  e2 80 9f   reversed-9 double quote
            // ʺ U+02ba    698  ca ba      modifier letter double prime

            // PLUS
            // + U+002b     43  2b         plus sign
            // ➕ U+2795  10133  e2 9e 95   heavy plus
            //   U+2064   8292  e2 81 a4   invisible plus |⁤|

            // MINUS
            // - U+002d     45  2d         hyphen
            // ‐ U+2010   8208  e2 80 90   dash
            // ‑ U+2011   8209  e2 80 91   non-breaking hyphen
            // – U+2013   8211  e2 80 93   en dash (also ext ASCII 0x96)
            // — U+2014   8212  e2 80 94   em dash (also ext ASCII 0x97)
            // − U+2212   8722  e2 88 92   minus sign
            // ➖ U+2796  10134  e2 9e 96   heavy minus

            // IGNORED
            //   U+00a0    160  c2 a0      non-breaking space
            //   U+2007   8199  e2 80 87   figure space | |
            //   U+2009   8201  e2 80 89   thin space   | |
            //   U+200a   8202  e2 80 8a   hair space   | |
            //   U+200b   8203  e2 80 8b   invisible space |​|
            //   U+202f   8239  e2 80 af   narrow space | |
            //   U+2063   8291  e2 81 a3   invisible separator |⁣|
            // « U+00ab    171  c2 ab      left guillemot (for cgi-bin)
            // » U+00bb    187  c2 bb      right guillemot (for cgi-bin)
            var dmsa = dms
                .Replace('°', 'd') // U+00b0 degree symbol
                .Replace('º', 'd') // U+00ba alt symbol
                .Replace('⁰', 'd') // U+2070 sup zero
                .Replace('˚', 'd') // U+02da ring above
                .Replace('∘', 'd') // U+2218 compose function

                .Replace('′', '\'') // U+2032 prime
                .Replace('‵', '\'') // U+2035 back prime
                .Replace('´', '\'') // U+00b4 acute accent
                .Replace('‘', '\'') // U+2018 left single quote
                .Replace('’', '\'') // U+2019 right single quote
                .Replace('‛', '\'') // U+201b reversed-9 single quote
                .Replace('ʹ', '\'') // U+02b9 modifier letter prime
                .Replace('ˊ', '\'') // U+02ca modifier letter acute accent
                .Replace('ˋ', '\'') // U+02cb modifier letter grave accent

                .Replace('″', '\"') // U+2033 double prime
                .Replace('‶', '\"') // U+2036 reversed double prime
                .Replace('˝', '\"') // U+02dd double acute accent
                .Replace('“', '\"') // U+201c left double quote
                .Replace('”', '\"') // U+201d right double quote
                .Replace('‟', '\"') // U+201f reversed-9 double quote
                .Replace('ʺ', '\"') // U+02ba modifier letter double prime

                .Replace('➕', '+') // U+2795 heavy plus
                .Replace('⁤', '+') // U+2064 invisible plus

                .Replace('‐', '-') // U+2010 dash
                .Replace('‑', '-') // U+2011 non-breaking hyphen
                .Replace('–', '-') // U+2013 en dash
                .Replace('—', '-') // U+2014 em dash
                .Replace('−', '-') // U+2212 minus sign
                .Replace('➖', '-') // U+2796 heavy minus

                .Replace(' ', '\0') // U+00a0 non-breaking space
                .Replace(' ', '\0') // U+2007 figure space
                .Replace(' ', '\0') // U+2007 thin space
                .Replace(' ', '\0') // U+200a hair space
                .Replace('​', '\0') // U+200b invisible space
                .Replace(' ', '\0') // U+202f narrow space
                .Replace('⁣', '\0') // U+2063 invisible separator

                .Replace('\xb0', 'd') // 0xb0 bare degree symbol
                .Replace('\xba', 'd') // 0xba bare alt symbol
                .Replace('*', 'd') // GRiD symbol for degree
                .Replace('`', '\'') // grave accent
                .Replace('\xb4', '\'') // 0xb4 bare acute accent
                                       // Don't implement these alternatives; they are only relevant for cgi-bin
                                       // replace(dmsa, "\x91",      '\'') // 0x91 ext ASCII left single quote
                                       // replace(dmsa, "\x92",      '\'') // 0x92 ext ASCII right single quote
                                       // replace(dmsa, "\x93",      '"' ) // 0x93 ext ASCII left double quote
                                       // replace(dmsa, "\x94",      '"' ) // 0x94 ext ASCII right double quote
                                       // replace(dmsa, "\x96",      '-' ) // 0x96 ext ASCII en dash
                                       // replace(dmsa, "\x97",      '-' ) // 0x97 ext ASCII em dash
                .Replace('\xa0', '\0') // 0xa0 bare non-breaking space
                .Replace("''", "\"") // '' -> "
                .Trim();

            // The trimmed string in [beg, end)
            double v = 0;
            int i = 0;
            var ind1 = HemisphereIndicator.None;

            // p is pointer to the next piece that needs decoding
            for (int p = 0, pb; p < dmsa.Length; p = pb, ++i)
            {
                var pa = p;
                // Skip over initial hemisphere letter (for i == 0)
                if (i == 0 && hemispheres_.IndexOf(dmsa[pa], StringComparison.OrdinalIgnoreCase) >= 0)
                    ++pa;
                // Skip over initial sign (checking for it if i == 0)
                if (i > 0 || (pa < dmsa.Length && signs_.IndexOf(dmsa[pa], StringComparison.OrdinalIgnoreCase) >= 0))
                    ++pa;
                // Find next sign
                pb = (int)Min((uint)dmsa.IndexOfAny(signs_.ToCharArray(), pa), dmsa.Length);
                var (x, ind2) = InternalDecode(dmsa.Substring(p, pb - p));
                v += x;
                if (ind1 == HemisphereIndicator.None)
                    ind1 = ind2;
                else if (!(ind2 == HemisphereIndicator.None || ind1 == ind2))
                    throw new GeographicException("Incompatible hemisphere specifier in " +
                                        dmsa.Substring(0, pb));
            }

            if (i == 0)
                throw new GeographicException("Empty or incomplete DMS string " + dmsa);
            return (v, ind1);
        }

        /// <summary>
        /// Convert <see cref="DMS"/> to an angle.
        /// </summary>
        /// <param name="d">degrees</param>
        /// <param name="m">arc minutes</param>
        /// <param name="s">arc seconds</param>
        /// <returns>Angle in degrees.</returns>
        /// <remarks>
        /// This does not propagate the sign on <paramref name="d"/> to the other components, so <c>-3d20'</c> would need to be represented as
        /// <c>-DMS.Decode(3.0, 20.0)</c> or <c>DMS.Decode(-3.0, -20.0)</c>.
        /// </remarks>
        public static double Decode(double d, double m = 0, double s = 0)
            => d + (m + s / 60) / 60;

        /// <summary>
        /// Convert a pair of strings to latitude and longitude.
        /// </summary>
        /// <param name="stra">first string.</param>
        /// <param name="strb">second string.</param>
        /// <param name="longfirst">
        /// if <see langword="true"/> assume longitude is given before latitude in the absence of hemisphere designators (default <see langword="false"/>).
        /// </param>
        /// <returns>latitude and longitude in degrees.</returns>
        /// <remarks>
        /// By default, the <i>lat</i> (resp., <i>lon</i>) is assigned to the results of decoding <paramref name="stra"/> (resp., <paramref name="strb"/>).
        /// However this is overridden if either <paramref name="stra"/> or <paramref name="strb"/> contain a latitude or longitude hemisphere designator
        /// (N, S, E, W).
        /// </remarks>
        public static (double lat, double lon) Decode(string stra, string strb, bool longfirst = false)
        {
            var (a, ia) = Decode(stra);
            var (b, ib) = Decode(strb);

            if (ia == HemisphereIndicator.None && ib == HemisphereIndicator.None)
            {
                // Default to lat, long unless longfirst
                ia = longfirst ? HemisphereIndicator.Longitude : HemisphereIndicator.Latitude;
                ib = longfirst ? HemisphereIndicator.Latitude : HemisphereIndicator.Longitude;
            }
            else if (ia == HemisphereIndicator.None)
                ia = (HemisphereIndicator)((int)HemisphereIndicator.Latitude + HemisphereIndicator.Longitude - ib);
            else if (ib == HemisphereIndicator.None)
                ib = (HemisphereIndicator)((int)HemisphereIndicator.Latitude + HemisphereIndicator.Longitude - ia);
            if (ia == ib)
                throw new GeographicException($"Both {stra.ToString()} and {strb.ToString()} interpreted as" +
                    $" {(ia == HemisphereIndicator.Latitude ? "latitudes" : "longitudes")}");
            double
              lat1 = ia == HemisphereIndicator.Latitude ? a : b,
              lon1 = ia == HemisphereIndicator.Latitude ? b : a;
            if (Abs(lat1) > 90)
                throw new GeographicException($"Latitude {lat1}d not in [-90d, 90d]");

            return (lat1, lon1);
        }

        /// <summary>
        /// Convert a string to an angle in degrees.
        /// </summary>
        /// <param name="angstr">input string.</param>
        /// <returns>Angle in degrees.</returns>
        /// <remarks>
        /// No hemisphere designator is allowed and no check is done on the range of the result.
        /// </remarks>
        public static double DecodeAngle(string angstr)
        {
            var (ang, ind) = Decode(angstr);
            if (ind != HemisphereIndicator.None)
                throw new GeographicException($"Arc angle {angstr} includes a hemisphere, N/E/W/S");
            return ang;
        }

        /// <summary>
        /// Convert a string to an azimuth in degrees.
        /// </summary>
        /// <param name="azistr">input string.</param>
        /// <returns>azimuth (degrees) reduced to the range [−180°, 180°].</returns>
        /// <remarks>
        /// A hemisphere designator E/W can be used; the result is multiplied by −1 if W is present.
        /// </remarks>
        public static double DecodeAzimuth(string azistr)
        {
            var (azi, ind) = Decode(azistr);
            if (ind == HemisphereIndicator.Latitude)
                throw new GeographicException("Azimuth " + azistr
                                    + " has a latitude hemisphere, N/S");
            return AngNormalize(azi);
        }

        /// <summary>
        /// Convert angle into a <see cref="DMS"/> string (using d, ', and ") selecting the trailing component based on the precision.
        /// </summary>
        /// <param name="angle">input angle (degrees)</param>
        /// <param name="prec">the precision relative to 1 degree.</param>
        /// <param name="ind"><see cref="HemisphereIndicator"/> value indicating additional formatting.</param>
        /// <param name="dmssep">if not <c>0</c>, use as the <see cref="DMS"/> separator character (instead of d, ', " delimiters).</param>
        /// <returns>formatted string</returns>
        /// <remarks>
        /// <paramref name="prec"/> indicates the precision relative to 1 degree, e.g., <paramref name="prec"/> = 3 gives a result accurate to 0.1' and
        /// <paramref name="prec"/> = 4 gives a result accurate to 1". <paramref name="ind"/> is interpreted as in 
        /// <see cref="Encode(double, TrailingUnit, int, HemisphereIndicator, char)"/> with the additional facility that <see cref="HemisphereIndicator.Number"/>
        /// represents angle as a number in fixed format with precision <paramref name="prec"/>.
        /// </remarks>
        public static string Encode(double angle, int prec, HemisphereIndicator ind = HemisphereIndicator.None, char dmssep = '\0')
            => ind == HemisphereIndicator.Number ? angle.ToFixedString(prec) :
                Encode(angle,
                       prec < 2 ? TrailingUnit.Degree : (prec < 4 ? TrailingUnit.Minute : TrailingUnit.Second),
                       prec < 2 ? prec : (prec < 4 ? prec - 2 : prec - 4),
                       ind, dmssep);

        /// <summary>
        /// Convert angle (in degrees) into a DMS string (using d, ', and ").
        /// </summary>
        /// <param name="angle">input angle (degrees)</param>
        /// <param name="trailing"><see cref="TrailingUnit"/> value indicating the trailing units of the string (this component is given as a decimal number if necessary).</param>
        /// <param name="prec">the number of digits after the decimal point for the trailing component.</param>
        /// <param name="ind"><see cref="HemisphereIndicator"/> value indicating additional formatting.</param>
        /// <param name="dmssep">if not <c>0</c>, use as the <see cref="DMS"/> separator character (instead of d, ', " delimiters).</param>
        /// <returns>formatted string</returns>
        /// <remarks>
        /// The interpretation of <paramref name="ind"/> is as follows:
        /// <list type="bullet">
        /// <item><paramref name="ind"/> == <see cref="HemisphereIndicator.None"/>,
        /// signed result no leading zeros on degrees except in the units place, e.g., <c>-8d03'</c>.</item>
        /// <item><paramref name="ind"/> == <see cref="HemisphereIndicator.Latitude"/>,
        /// trailing N or S hemisphere designator, no sign, pad degrees to 2 digits, e.g., <c>08d03'S</c>.
        /// </item>
        /// <item><paramref name="ind"/> == <see cref="HemisphereIndicator.Longitude"/>,
        /// trailing E or W hemisphere designator, no sign, pad degrees to 3 digits, e.g., <c>008d03'W</c>.
        /// </item>
        /// <item><paramref name="ind"/> == <see cref="HemisphereIndicator.Azimuth"/>,
        /// convert to the range [0, 360°), no sign, pad degrees to 3 digits, e.g., <c>351d57'</c>.
        /// </item>
        /// </list>
        /// The integer parts of the minutes and seconds components are always given with 2 digits.
        /// </remarks>
        public static string Encode(double angle, TrailingUnit trailing, int prec, HemisphereIndicator ind = HemisphereIndicator.None, char dmssep = '\0')
        {
            // Assume check on range of input angle has been made by calling
            // routine (which might be able to offer a better diagnostic).
            if (!IsFinite(angle))
                return angle < 0 ? "-inf" :
                  (angle > 0 ? "inf" : "nan");

            // 15 - 2 * trailing = ceiling(log10(2^53/90/60^trailing)).
            // This suffices to give full real precision for numbers in [-90,90]
            prec = Min(15 + 0 - 2 * (int)trailing, prec);
            var scale = 1d;
            for (var i = 0; i < (int)trailing; ++i)
                scale *= 60;
            for (var i = 0; i < prec; ++i)
                scale *= 10;
            if (ind == HemisphereIndicator.Azimuth)
                angle -= Floor(angle / 360) * 360;
            int sign = angle < 0 ? -1 : 1;
            angle *= sign;

            // Break off integer part to preserve precision in manipulation of
            // fractional part.
            double
              idegree = Floor(angle),
              fdegree = (angle - idegree) * scale + 0.5;
            {
                // Implement the "round ties to even" rule
                var f = Floor(fdegree);
                fdegree = (f == fdegree && (f % 2) == 1) ? f - 1 : f;
            }
            fdegree /= scale;
            if (fdegree >= 1)
            {
                idegree += 1;
                fdegree -= 1;
            }
            Span<double> pieces = stackalloc[] { fdegree, 0, 0 };
            for (var i = 1; i <= (int)trailing; ++i)
            {
                double
                  ip = Floor(pieces[i - 1]),
                  fp = pieces[i - 1] - ip;
                pieces[i] = fp * 60;
                pieces[i - 1] = ip;
            }
            pieces[0] += idegree;
            var s = new StringBuilder();
            if (ind == HemisphereIndicator.None && sign < 0)
                s.Append('-');
            var w = 0;
            switch (trailing)
            {
                case TrailingUnit.Degree:
                    if (ind != HemisphereIndicator.None)
                    {
                        w = 1 + Min((int)ind, 2) + prec + (prec != 0 ? 1 : 0);
                    }
                    s.Append(pieces[0].ToFixedString(prec).PadLeft(w, '0'));
                    // Don't include degree designator (d) if it is the trailing component.
                    break;
                default:
                    if (ind != HemisphereIndicator.None)
                    {
                        w = 1 + Min((int)ind, 2);
                    }
                    s.Append(pieces[0].ToString().PadLeft(w, '0'))
                     .Append(dmssep != 0 ? dmssep : char.ToLower(dmsindicators_[0]));
                    switch (trailing)
                    {
                        case TrailingUnit.Minute:
                            w = 2 + prec + (prec != 0 ? 1 : 0);
                            s.Append(pieces[1].ToFixedString(prec).PadLeft(w, '0'));
                            if (dmssep != 0)
                                s.Append(char.ToLower(dmsindicators_[1]));
                            break;
                        case TrailingUnit.Second:
                            s.AppendFormat("{0:D2}", (int)pieces[1])
                             .Append(dmssep != 0 ? dmssep : char.ToLower(dmsindicators_[1]))
                             .Append(pieces[2].ToFixedString(prec).PadLeft(2 + prec + (prec != 0 ? 1 : 0), '0'));
                            if (dmssep == 0)
                                s.Append(char.ToLower(dmsindicators_[2]));
                            break;
                        default:
                            break;
                    }
                    break;
            }
            if (ind != HemisphereIndicator.None && ind != HemisphereIndicator.Azimuth)
                s.Append(hemispheres_[(ind == HemisphereIndicator.Latitude ? 0 : 2) + (sign < 0 ? 0 : 1)]);
            return s.ToString();
        }

        /// <summary>
        /// Split angle into degrees and minutes.
        /// </summary>
        /// <param name="ang">angle (degrees)</param>
        /// <returns>degrees and arc minutes.</returns>
        public static (double degrees, double minutes) EncodeDM(double ang)
        {
            var d = (int)ang;
            return (d, 60 * (ang - d));
        }

        /// <summary>
        /// Split angle into degrees and minutes and seconds.
        /// </summary>
        /// <param name="ang">angle (degrees)</param>
        /// <returns>degrees and arc minutes and arc seconds</returns>
        public static (double degrees, double minutes, double seconds) Encode(double ang)
        {
            var d = (int)ang;
            ang = 60 * (ang - d);
            var m = (int)ang;

            return (d, m, 60 * (ang - m));
        }
    }
}
