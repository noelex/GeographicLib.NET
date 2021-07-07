using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

using static System.Math;
using static GeographicLib.MathEx;
using static GeographicLib.Macros;

namespace GeographicLib
{
    /// <summary>
    /// Looking up the height of the geoid above the ellipsoid.
    /// </summary>
    /// <remarks>
    /// This class evaluates the height of one of the standard geoids, EGM84, EGM96, or EGM2008 by bilinear or cubic interpolation
    /// into a rectangular grid of data.These geoid models are documented in
    /// <list type="bullet">
    /// <item>
    /// EGM84: <a href="https://earth-info.nga.mil/index.php?dir=wgs84&amp;action=wgs84#tab_egm84"></a>
    /// </item>
    /// <item>
    /// EGM96: <a href="https://earth-info.nga.mil/index.php?dir=wgs84&amp;action=wgs84#tab_egm96"></a>
    /// </item>
    /// <item>
    /// EGM2008: <a href="https://earth-info.nga.mil/index.php?dir=wgs84&amp;action=wgs84#tab_egm2008"></a>
    /// </item>
    /// </list>
    /// The geoids are defined in terms of spherical harmonics.
    /// However in order to provide a quick and flexible method of evaluating the geoid heights,
    /// this class evaluates the height by interpolation into a grid of precomputed values.
    /// <para>
    /// The height of the geoid above the ellipsoid, <i>N</i>, is sometimes called the geoid undulation.
    /// It can be used to convert a height above the ellipsoid, <i>h</i>, to the corresponding height above the geoid
    /// (the orthometric height, roughly the height above mean sea level), <i>H</i>, using the relations
    /// </para>
    /// <para>
    /// <i>h</i> = <i>N</i> + <i>H</i>;   <i>H</i> = −<i>N</i> + <i>h</i>.
    /// </para>
    /// <para>
    /// See <a href="https://geographiclib.sourceforge.io/html/geoid.html">Geoid height</a> for details of how to install the data sets, the data format, estimates of the interpolation errors, and how to use caching.
    /// </para>
    /// <para>
    /// This class is typically <i>not</i> thread safe in that a single instantiation cannot be safely used by multiple threads because of the way
    /// the object reads the data set and because it maintains a single-cell cache. If multiple threads need to calculate geoid heights they should
    /// all construct thread-local instantiations. Alternatively, set the optional threadsafe parameter to true in the constructor. This causes the
    /// constructor to read all the data into memory and to turn off the single-cell caching which results in a <see cref="Geoid"/> object which is thread safe.
    /// </para>
    /// </remarks>
    public sealed class Geoid : IEllipsoid, IDisposable
    {
        private static readonly uint
            pixel_size_ = 2,
            pixel_max_ = 0xffffu;

        private const int stencilsize_ = 12;
        private const int nterms_ = ((3 + 1) * (3 + 2)) / 2; // for a cubic fit

        private const int c0_ = 240;
        private const int c0n_ = 372;
        private const int c0s_ = 372;

        private static readonly Memory<int> c3_ = new[]
        {
          9, -18, -88,    0,  96,   90,   0,   0, -60, -20,
         -9,  18,   8,    0, -96,   30,   0,   0,  60, -20,
          9, -88, -18,   90,  96,    0, -20, -60,   0,   0,
        186, -42, -42, -150, -96, -150,  60,  60,  60,  60,
         54, 162, -78,   30, -24,  -90, -60,  60, -60,  60,
         -9, -32,  18,   30,  24,    0,  20, -60,   0,   0,
         -9,   8,  18,   30, -96,    0, -20,  60,   0,   0,
         54, -78, 162,  -90, -24,   30,  60, -60,  60, -60,
        -54,  78,  78,   90, 144,   90, -60, -60, -60, -60,
          9,  -8, -18,  -30, -24,    0,  20,  60,   0,   0,
         -9,  18, -32,    0,  24,   30,   0,   0, -60,  20,
          9, -18,  -8,    0, -24,  -30,   0,   0,  60,  20,
        };
        private static readonly ReadOnlyMemory<int> c3n_ = new[]
        {
            0, 0, -131, 0,  138,  144, 0,   0, -102, -31,
            0, 0,    7, 0, -138,   42, 0,   0,  102, -31,
            62, 0,  -31, 0,    0,  -62, 0,   0,    0,  31,
        124, 0,  -62, 0,    0, -124, 0,   0,    0,  62,
        124, 0,  -62, 0,    0, -124, 0,   0,    0,  62,
            62, 0,  -31, 0,    0,  -62, 0,   0,    0,  31,
            0, 0,   45, 0, -183,   -9, 0,  93,   18,   0,
            0, 0,  216, 0,   33,   87, 0, -93,   12, -93,
            0, 0,  156, 0,  153,   99, 0, -93,  -12, -93,
            0, 0,  -45, 0,   -3,    9, 0,  93,  -18,   0,
            0, 0,  -55, 0,   48,   42, 0,   0,  -84,  31,
            0, 0,   -7, 0,  -48,  -42, 0,   0,   84,  31,
        };
        private static readonly ReadOnlyMemory<int> c3s_ = new[]{
             18,  -36, -122,   0,  120,  135, 0,   0,  -84, -31,
            -18,   36,   -2,   0, -120,   51, 0,   0,   84, -31,
             36, -165,  -27,  93,  147,   -9, 0, -93,   18,   0,
            210,   45, -111, -93,  -57, -192, 0,  93,   12,  93,
            162,  141,  -75, -93, -129, -180, 0,  93,  -12,  93,
            -36,  -21,   27,  93,   39,    9, 0, -93,  -18,   0,
              0,    0,   62,   0,    0,   31, 0,   0,    0, -31,
              0,    0,  124,   0,    0,   62, 0,   0,    0, -62,
              0,    0,  124,   0,    0,   62, 0,   0,    0, -62,
              0,    0,   62,   0,    0,   31, 0,   0,    0, -31,
            -18,   36,  -64,   0,   66,   51, 0,   0, -102,  31,
             18,  -36,    2,   0,  -66,  -51, 0,   0,  102,  31,
          };

        private readonly string _name, _dir, _filename;
        private readonly bool _cubic;
        private readonly double _eps;
        private readonly Stream _file;
        private readonly double _rlonres, _rlatres;
        private readonly string _description;
        private readonly DateTime? _datetime;
        private readonly double _offset, _scale, _maxerror, _rmserror;
        private readonly int _width, _height;
        private readonly long _datastart, _swidth;
        private readonly bool _threadsafe;

        private List<Memory<ushort>> _data=new List<Memory<ushort>>();
        private bool _cache;
        private int _xoffset, _yoffset, _xsize, _ysize;
        private int _ix, _iy;
        private double _v00, _v01, _v10, _v11;
        private Memory<double> _t = new double[nterms_];

        static Geoid()
        {
            const string GEOGRAPHICLIB_GEOID_DEFAULT_NAME = "egm96-5";

            string GetGeoidPath()
            {
                var path = Environment.GetEnvironmentVariable("GEOGRAPHICLIB_GEOID_PATH");
                if (!string.IsNullOrWhiteSpace(path))
                    return path;
                path = Environment.GetEnvironmentVariable("GEOGRAPHICLIB_DATA");

                return Path.Combine((!string.IsNullOrWhiteSpace(path) ? path : GEOGRAPHICLIB_DATA), "geoids");
            }

            string GetGeoidName()
            {
                var name = Environment.GetEnvironmentVariable("GEOGRAPHICLIB_GEOID_NAME");
                return !string.IsNullOrWhiteSpace(name) ? name : GEOGRAPHICLIB_GEOID_DEFAULT_NAME;
            }

            DefaultGeoidPath = GetGeoidPath();
            DefaultGeoidName = GetGeoidName();
        }

        /// <summary>
        /// Construct a geoid.
        /// </summary>
        /// <param name="name">the name of the geoid.</param>
        /// <param name="path">directory for data file.</param>
        /// <param name="cubic">interpolation method; <see langword="false"/> means bilinear, <see langword="true"/> (the default) means cubic.</param>
        /// <param name="threadsafe">if <see langword="true"/>, construct a thread safe object. The default is <see langword="false"/></param>
        /// <remarks>
        /// The data file is formed by appending ".pgm" to the name.
        /// If <paramref name="path"/> is specified (and is non-empty), then the file is loaded from directory, path.
        /// Otherwise the path is given by <see cref="DefaultGeoidPath"/>.
        /// If the <paramref name="threadsafe"/> parameter is <see langword="true"/>, the data set is read into memory, the data file is closed,
        /// and single-cell caching is turned off; this results in a <see cref="Geoid"/> object which is thread safe.
        /// </remarks>
        public Geoid(string name, string path = "", bool cubic = true, bool threadsafe = false)
        {
            _name = name;
            _dir = path;
            _cubic = cubic;
            //_a = Constants.WGS84_a;
            //_e2 = (2 - Constants.WGS84_f) * Constants.WGS84_f;
            //_degree = Degree;
            _eps = Sqrt(DBL_EPSILON);
            _threadsafe = false;

            Debug.Assert(pixel_size_ == Macros.GEOGRAPHICLIB_GEOID_PGM_PIXEL_WIDTH, "pixel_t has the wrong size");

            if (string.IsNullOrWhiteSpace(_dir))
            {
                _dir = DefaultGeoidPath;
            }

            _filename = Path.Combine(_dir, Path.ChangeExtension(_name, pixel_size_ != 4 ? "pgm" : "pgm4"));
            _file = File.OpenRead(_filename);

            using (var sr = new StreamReader(_file, Encoding.UTF8, true, bufferSize: 1, leaveOpen: true))
            {
                var s = sr.ReadLine();

                if (s != "P5")
                {
                    throw new GeographicException("File not in PGM format: " + _filename);
                }

                _offset = double.MaxValue;
                _scale = 0;
                _maxerror = _rmserror = -1;
                _description = "NONE";
                _datetime = null;

                while ((s = sr.ReadLine()) != null)
                {
                    if (s.StartsWith("#"))
                    {
                        var match = Regex.Match(s, @"^#\s+([A-Za-z]+)\s+(.+)$");
                        if (!match.Success)
                        {
                            continue;
                        }

                        var key = match.Groups[1].Value;
                        if (key == "Description")
                        {
                            _description = match.Groups[2].Value.Trim();
                        }
                        else if (key == "DateTime")
                        {
                            _datetime = System.DateTime.Parse( match.Groups[2].Value.Trim());
                        }
                        else if (key == "Offset")
                        {
                            if (!double.TryParse(match.Groups[2].Value.Trim(), out _offset))
                            {
                                throw new GeographicException("Error reading offset: " + _filename);
                            }
                        }
                        else if (key == "Scale")
                        {
                            if (!double.TryParse(match.Groups[2].Value.Trim(), out _scale))
                            {
                                throw new GeographicException("Error reading scale: " + _filename);
                            }
                        }
                        else if (key == (_cubic ? "MaxCubicError" : "MaxBilinearError"))
                        {
                            // It's not an error if the error can't be read
                            double.TryParse(match.Groups[2].Value.Trim(), out _maxerror);
                        }
                        else if (key == (_cubic ? "RMSCubicError" : "RMSBilinearError"))
                        {
                            // It's not an error if the error can't be read
                            double.TryParse(match.Groups[2].Value.Trim(), out _rmserror);
                        }
                    }
                    else
                    {
                        var items = s.Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                        if (items.Length != 2 || !int.TryParse(items[0], out _width) || !int.TryParse(items[1], out _height))
                        {
                            throw new GeographicException("Error reading raster size: " + _filename);
                        }
                        break;
                    }
                }


                if (!uint.TryParse(s = sr.ReadLine(), out var maxval))
                    throw new GeographicException("Error reading maxval: " + _filename);
                if (maxval != pixel_max_)
                    throw new GeographicException("Incorrect value of maxval: " + _filename);

                // HACK: Get start position of binary data.
                sr.BaseStream.Seek(0, SeekOrigin.Begin);
                sr.DiscardBufferedData();
                var buff = new char[1024];
                var sp = buff.AsSpan();
                sr.ReadBlock(buff, 0, buff.Length);
                var end = sp.IndexOf((s + '\n').AsSpan()) + s.Length + 1;

                // Add 1 for whitespace after maxval
                _datastart = Encoding.UTF8.GetByteCount(sp.Slice(0, end).ToArray()); // +1 ?
                _swidth = _width;
            }

            if (_offset == double.MaxValue)
                throw new GeographicException("Offset not set: " + _filename);
            if (_scale == 0)
                throw new GeographicException("Scale not set " + _filename);
            if (_scale < 0)
                throw new GeographicException("Scale must be positive " + _filename);
            if (_height < 2 || _width < 2)
                // Coarsest grid spacing is 180deg.
                throw new GeographicException("Raster size too small " + _filename);
            if ((_width & 1) != 0)
                // This is so that longitude grids can be extended thru the poles.
                throw new GeographicException("Raster width is odd " + _filename);
            if ((_height & 1) == 0)
                // This is so that latitude grid includes the equator.
                throw new GeographicException("Raster height is even " + _filename);

            _file.Seek(0, SeekOrigin.End);
            if (_datastart + pixel_size_ * _swidth * _height != _file.Position)
                // Possibly this test should be "<" because the file contains, e.g., a
                // second image.  However, for now we are more strict.
                throw new GeographicException("File has the wrong length " + _filename);
            _rlonres = _width / 360.0;
            _rlatres = (_height - 1) / 180.0;
            _cache = false;
            _ix = _width;
            _iy = _height;

            // Ensure that file errors throw exceptions
            if (threadsafe)
            {
                CacheAll();
                _file.Close();
                _threadsafe = true;
            }
        }

        /// <summary>
        /// Set up a cache.
        /// </summary>
        /// <param name="south">latitude (degrees) of the south edge of the cached area.</param>
        /// <param name="west">longitude (degrees) of the west edge of the cached area.</param>
        /// <param name="north">latitude (degrees) of the north edge of the cached area.</param>
        /// <param name="east">longitude (degrees) of the east edge of the cached area.</param>
        /// <remarks>
        /// Cache the data for the specified "rectangular" area bounded by the parallels <paramref name="south"/> and <paramref name="north"/>
        /// and the meridians <paramref name="west"/> and <paramref name="east"/>.
        /// <paramref name="east"/> is always interpreted as being east of <paramref name="west"/>,
        /// if necessary by adding 360° to its value. <paramref name="south"/> and <paramref name="north"/> should be in the range [−90°, 90°].
        /// </remarks>
        public void CacheArea(double south, double west, double north, double east)
        {
            if (_threadsafe)
                throw new GeographicException("Attempt to change cache of threadsafe Geoid");
            if (south > north)
            {
                CacheClear();
                return;
            }
            south = LatFix(south);
            north = LatFix(north);
            west = AngNormalize(west); // west in [-180, 180)
            east = AngNormalize(east);
            if (east <= west)
                east += 360;              // east - west in (0, 360]
            int
              iw = (int)Floor(west * _rlonres),
              ie = (int)Floor(east * _rlonres),
              in_ = (int)Floor(-north * _rlatres) + (_height - 1) / 2,
              is_ = (int)Floor(-south * _rlatres) + (_height - 1) / 2;

            in_ = Max(0, Min(_height - 2, in_));
            is_ = Max(0, Min(_height - 2, is_));
            is_ += 1;
            ie += 1;

            if (_cubic)
            {
                in_ -= 1;
                is_ += 1;
                iw -= 1;
                ie += 1;
            }
            if (ie - iw >= _width - 1)
            {
                // Include entire longitude range
                iw = 0;
                ie = _width - 1;
            }
            else
            {
                ie += iw < 0 ? _width : (iw >= _width ? -_width : 0);
                iw += iw < 0 ? _width : (iw >= _width ? -_width : 0);
            }

            _xsize = ie - iw + 1;
            _ysize = is_ - in_ + 1;
            _xoffset = iw;
            _yoffset = in_;


            _data.Clear();
            for (int iy = _ysize; iy-- > 0;)
            {
                _data.Add(new ushort[_xsize]);
            }

            try
            {
                for (int iy = in_; iy <= is_; ++iy)
                {
                    int iy1 = iy, iw1 = iw;
                    if (iy < 0 || iy >= _height)
                    {
                        // Allow points "beyond" the poles to support interpolation
                        iy1 = iy1 < 0 ? -iy1 : 2 * (_height - 1) - iy1;
                        iw1 += _width / 2;
                        if (iw1 >= _width)
                            iw1 -= _width;
                    }
                    int xs1 = Min(_width - iw1, _xsize);
                    filepos(iw1, iy1);

                    Utility.ReadArray(_file, _data[iy - in_].Span.Slice(0, xs1), true);
                    if (xs1 < _xsize)
                    {
                        // Wrap around longitude = 0
                        filepos(0, iy1);
                        Utility.ReadArray(_file, _data[iy - in_].Span.Slice(xs1, _xsize - xs1), true);
                    }
                }
                _cache = true;
            }
            catch (Exception e)
            {
                CacheClear();
                throw new GeographicException("Error filling cache: " + e.Message);
            }
        }

        /// <summary>
        /// Cache all the data.
        /// </summary>
        /// <remarks>
        /// On most computers, this is fast for data sets with grid resolution of 5' or coarser.
        /// For a 1' grid, the required RAM is 450MB; a 2.5' grid needs 72MB; and a 5' grid needs 18MB.
        /// </remarks>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void CacheAll() => CacheArea(-90, 0, 90, 360);

        /// <summary>
        /// Clear the cache. This never throws an error. (This does nothing with a thread safe <see cref="Geoid"/>.)
        /// </summary>
        public void CacheClear()
        {
            if (!_threadsafe)
            {
                _cache = false;
                try
                {
                    _data.Clear();
                }
                catch
                {
                }
            }
        }

        /// <summary>
        /// Compute the geoid height at a point.
        /// </summary>
        /// <param name="lat">latitude of the point (degrees).</param>
        /// <param name="lon">longitude of the point (degrees).</param>
        /// <returns>the height of the geoid above the ellipsoid (meters).</returns>
        /// <remarks>The latitude should be in [−90°, 90°].</remarks>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double Evaluate(double lat, double lon)
            => height(lat, lon);

        /// <summary>
        /// Compute the geoid height at a point.
        /// </summary>
        /// <param name="coords">coordinate of the point.</param>
        /// <returns>the height of the geoid above the ellipsoid (meters).</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double Evaluate(GeoCoords coords) => Evaluate(coords.Latitude, coords.Longitude);

        /// <summary>
        /// Convert a height above the geoid to a height above the ellipsoid and vice versa.
        /// </summary>
        /// <param name="lat">latitude of the point (degrees).</param>
        /// <param name="lon">longitude of the point (degrees).</param>
        /// <param name="h">height of the point (degrees).</param>
        /// <param name="d">a <see cref="ConvertFlag"/> specifying the direction of the conversion;
        /// <see cref="ConvertFlag.GeoidToEllipsoid"/> means convert a height above the geoid to a height above the ellipsoid;
        /// <see cref="ConvertFlag.EllipsoidToGeoid"/> means convert a height above the ellipsoid to a height above the geoid.</param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double ConvertHeight(double lat, double lon, double h, ConvertFlag d)
            => h + (int)d * height(lat, lon);

        /// <summary>
        /// Gets a value representing decription of the geoid, if available, in the data file; if absent, return "NONE".
        /// </summary>
        public string Description => _description;

        /// <summary>
        /// Gets a value representing date of the geoid, if available, in the data file; if absent, return <see langword="null"/>.
        /// </summary>
        public DateTime? DateTime => _datetime;

        /// <summary>
        /// Gets a value representing the name used to load the geoid data (from the first argument of the constructor).
        /// </summary>
        public string GeoidNmae => _name;

        /// <summary>
        /// Gets a value representing the full file name used to load the geoid data.
        /// </summary>
        public string GeoidFile => _filename;

        /// <summary>
        /// Gets a value representing the directory used to load the geoid data.
        /// </summary>
        public string GeoidDirectory => _dir;

        /// <summary>
        /// Gets a value representing the interpolation method of the geoid data ("cubic" or "bilinear").
        /// </summary>
        public string Interpolation => _cubic ? "cubic" : "bilinear";

        /// <summary>
        /// Gets a value representing estimate of the maximum interpolation and quantization error (meters).
        /// </summary>
        /// <remarks>
        /// This relies on the value being stored in the data file. If the value is absent, return −1.
        /// </remarks>
        public double MaxError => _maxerror;

        /// <summary>
        /// Gets a value representing estimate of the RMS interpolation and quantization error (meters).
        /// </summary>
        /// <remarks>
        /// This relies on the value being stored in the data file. If the value is absent, return −1.
        /// </remarks>
        public double RMSError => _rmserror;

        /// <summary>
        /// Gets a value representing offset (meters) of the geoid data.
        /// </summary>
        /// <remarks>
        /// This in used in converting from the pixel values in the data file to geoid heights.
        /// </remarks>
        public double Offset => _offset;

        /// <summary>
        /// Gets a value representing scale (meters) of the geoid data.
        /// </summary>
        /// <remarks>
        /// This in used in converting from the pixel values in the data file to geoid heights.
        /// </remarks>
        public double Scale => _scale;

        /// <summary>
        /// Gets a value representing whether current <see cref="Geoid"/> instance is thread safe,
        /// </summary>
        public bool IsThreadSafe => _threadsafe;

        /// <summary>
        /// Gets a value representing whether data cache is enabled.
        /// </summary>
        public bool IsCacheEnabled => _cache;

        /// <summary>
        /// Gets a value representing west edge of the cached area; the cache includes this edge.
        /// </summary>
        public double CacheWest
            => _cache ? ((_xoffset + (_xsize == _width ? 0 : (_cubic ? 1 : 0))
                        + _width / 2) % _width - _width / 2) / _rlonres :
               0;

        /// <summary>
        /// Gets a value representing east edge of the cached area; the cache excludes this edge.
        /// </summary>
        public double CacheEast
            => _cache ?
                    CacheWest +
                    (_xsize - (_xsize == _width ? 0 : 1 + 2 * (_cubic ? 1 : 0))) / _rlonres :
                    0;

        /// <summary>
        /// Gets a value representing north edge of the cached area; the cache includes this edge.
        /// </summary>
        public double CacheNorth
            => _cache ? 90 - (_yoffset + (_cubic ? 1 : 0)) / _rlatres : 0;

        /// <summary>
        /// Gets a value representing south edge of the cached area; the cache excludes this edge unless it's the south pole.
        /// </summary>
        public double CacheSouth
            => _cache ? 90 - (_yoffset + _ysize - 1 - (_cubic ? 1 : 0)) / _rlatres : 0;

        /// <inheritdoc/>
        public double EquatorialRadius => Constants.WGS84_a;

        /// <inheritdoc/>
        public double Flattening => Constants.WGS84_f;

        /// <summary>
        /// Gets a value representing the default path for geoid data files.
        /// </summary>
        /// <remarks>
        /// This is the value of the environment variable GEOGRAPHICLIB_GEOID_PATH, if set;
        /// otherwise, it is $GEOGRAPHICLIB_DATA/geoids if the environment variable GEOGRAPHICLIB_DATA is set;
        /// otherwise, it is a compile-time default (/usr/local/share/GeographicLib/geoids on non-Windows systems and C:/ProgramData/GeographicLib/geoids on Windows systems).
        /// </remarks>
        public static string DefaultGeoidPath { get; }

        /// <summary>
        /// Gets a value representing the default name for the geoid.
        /// </summary>
        /// <remarks>
        /// This is the value of the environment variable GEOGRAPHICLIB_GEOID_NAME, if set; otherwise, it is "egm96-5".
        /// The <see cref="Geoid"/> class does not use this function;
        /// it is just provided as a convenience for a calling program when constructing a <see cref="Geoid"/> object.
        /// </remarks>
        public static string DefaultGeoidName { get; }

        private void filepos(int ix, int iy)
        {
            _file.Seek(
                 (long)(_datastart +
                   pixel_size_ * ((uint)iy * _swidth + (uint)ix)), SeekOrigin.Begin);
        }

        private double rawval(int ix, int iy)
        {
            if (ix < 0)
                ix += _width;
            else if (ix >= _width)
                ix -= _width;
            if (_cache && iy >= _yoffset && iy < _yoffset + _ysize &&
                ((ix >= _xoffset && ix < _xoffset + _xsize) ||
                 (ix + _width >= _xoffset && ix + _width < _xoffset + _xsize)))
            {
                return (double)_data[iy - _yoffset].Span
                            [ix >= _xoffset ? ix - _xoffset : ix + _width - _xoffset];
            }
            else
            {
                if (iy < 0 || iy >= _height)
                {
                    iy = iy < 0 ? -iy : 2 * (_height - 1) - iy;
                    ix += (ix < _width / 2 ? 1 : -1) * _width / 2;
                }
                try
                {
                    filepos(ix, iy);
                    // initial values to suppress warnings in case get fails
                    byte a = (byte)_file.ReadByte(), b = (byte)_file.ReadByte();

                    var r = ((a) << 8) | (b);
                    if (pixel_size_ == 4)
                    {
                        a = (byte)_file.ReadByte(); b = (byte)_file.ReadByte();
                        r = (r << 16) | ((a) << 8) | (b);
                    }
                    return r;
                }
                catch (Exception e)
                {
                    // throw GeographicErr("Error reading " + _filename + ": "
                    //                      + e.what());
                    // triggers complaints about the "binary '+'" under Visual Studio.
                    // So use '+=' instead.
                    var err = "Error reading ";
                    err += _filename;
                    err += ": ";
                    err += e.Message;
                    throw new GeographicException(err);
                }
            }
        }

        private double height(double lat, double lon)
        {
            lat = LatFix(lat);
            if (double.IsNaN(lat) || double.IsNaN(lon))
            {
                return double.NaN;
            }
            lon = AngNormalize(lon);
            double
              fx = lon * _rlonres,
              fy = -lat * _rlatres;
            int
              ix = (int)Floor(fx),
              iy = Min((_height - 1) / 2 - 1, (int)Floor(fy));
            fx -= ix;
            fy -= iy;
            iy += (_height - 1) / 2;
            ix += ix < 0 ? _width : (ix >= _width ? -_width : 0);
            double v00 = 0, v01 = 0, v10 = 0, v11 = 0;
            Span<double> t = stackalloc double[nterms_];

            if (_threadsafe || !(ix == _ix && iy == _iy))
            {
                if (!_cubic)
                {
                    v00 = rawval(ix, iy);
                    v01 = rawval(ix + 1, iy);
                    v10 = rawval(ix, iy + 1);
                    v11 = rawval(ix + 1, iy + 1);
                }
                else
                {
                    Span<double> v = stackalloc double[stencilsize_];

                    int k = 0;
                    v[k++] = rawval(ix, iy - 1);
                    v[k++] = rawval(ix + 1, iy - 1);
                    v[k++] = rawval(ix - 1, iy);
                    v[k++] = rawval(ix, iy);
                    v[k++] = rawval(ix + 1, iy);
                    v[k++] = rawval(ix + 2, iy);
                    v[k++] = rawval(ix - 1, iy + 1);
                    v[k++] = rawval(ix, iy + 1);
                    v[k++] = rawval(ix + 1, iy + 1);
                    v[k++] = rawval(ix + 2, iy + 1);
                    v[k++] = rawval(ix, iy + 2);
                    v[k++] = rawval(ix + 1, iy + 2);

                    var c3x = iy == 0 ? c3n_.Span : (iy == _height - 2 ? c3s_.Span : c3_.Span);
                    int c0x = iy == 0 ? c0n_ : (iy == _height - 2 ? c0s_ : c0_);
                    for (var i = 0; i < nterms_; ++i)
                    {
                        t[i] = 0;
                        for (var j = 0; j < stencilsize_; ++j)
                            t[i] += v[j] * c3x[nterms_ * j + i];
                        t[i] /= c0x;
                    }
                }
            }
            else
            { // same cell; used cached coefficients
                if (!_cubic)
                {
                    v00 = _v00;
                    v01 = _v01;
                    v10 = _v10;
                    v11 = _v11;
                }
                else
                    _t.Span.CopyTo(t);
            }
            if (!_cubic)
            {
                double
                  a = (1 - fx) * v00 + fx * v01,
                  b = (1 - fx) * v10 + fx * v11,
                  c = (1 - fy) * a + fy * b,
                  h = _offset + _scale * c;
                if (!_threadsafe)
                {
                    _ix = ix;
                    _iy = iy;
                    _v00 = v00;
                    _v01 = v01;
                    _v10 = v10;
                    _v11 = v11;
                }
                return h;
            }
            else
            {
                double h = t[0] + fx * (t[1] + fx * (t[3] + fx * t[6])) +
                  fy * (t[2] + fx * (t[4] + fx * t[7]) +
                       fy * (t[5] + fx * t[8] + fy * t[9]));
                h = _offset + _scale * h;
                if (!_threadsafe)
                {
                    _ix = ix;
                    _iy = iy;
                    t.CopyTo(_t.Span);
                }
                return h;
            }
        }

        /// <summary>
        /// Release all underlying resource used by current <see cref="Geoid"/> instance.
        /// </summary>
        public void Dispose()
        {
            _file.Dispose();
            _data.Clear();
        }
    }
}
