﻿using GeographicLib.SphericalHarmonics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.InteropServices;
using System.Text;
using static GeographicLib.Macros;
using static GeographicLib.MathEx;
using static System.Math;

namespace GeographicLib
{
    /// <summary>
    /// Model of the earth's magnetic field. 
    /// </summary>
    /// <remarks>
    /// Evaluate the earth's magnetic field according to a model. 
    /// At present only internal magnetic fields are handled.
    /// These are due to the earth's code and crust; these vary slowly (over many years).
    /// Excluded are the effects of currents in the ionosphere and magnetosphere which have daily and annual variations.
    /// <para>See <a href="https://geographiclib.sourceforge.io/html/magnetic.html">Magnetic models</a> 
    /// for details of how to install the magnetic models and the data format.</para>
    /// <para>
    /// See
    /// <list type="bullet">
    /// <item>
    /// General information:
    /// <list type="bullet">
    /// <item><a href="http://geomag.org/models/index.html"></a></item>
    /// </list>
    /// </item>
    /// <item>
    /// WMM2010:
    /// <list type="bullet">
    /// <item><a href="https://ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml"></a></item>
    /// <item><a href="https://ngdc.noaa.gov/geomag/WMM/data/WMM2010/WMM2010COF.zip"></a></item>
    /// </list>
    /// </item>
    /// <item>
    /// WMM2015 (deprecated):
    /// <list type="bullet">
    /// <item><a href="https://ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml"></a></item>
    /// <item><a href="https://ngdc.noaa.gov/geomag/WMM/data/WMM2015/WMM2015COF.zip"></a></item>
    /// </list>
    /// </item>
    /// <item>
    /// WMM2015V2:
    /// <list type="bullet">
    /// <item><a href="https://ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml"></a></item>
    /// <item><a href="https://ngdc.noaa.gov/geomag/WMM/data/WMM2015/WMM2015v2COF.zip"></a></item>
    /// </list>
    /// </item>
    /// <item>
    /// WMM2020:
    /// <list type="bullet">
    /// <item><a href="https://ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml"></a></item>
    /// <item><a href="https://ngdc.noaa.gov/geomag/WMM/data/WMM2020/WMM2020COF.zip"></a></item>
    /// </list>
    /// </item>
    /// <item>
    /// IGRF11:
    /// <list type="bullet">
    /// <item><a href="https://ngdc.noaa.gov/IAGA/vmod/igrf.html"></a></item>
    /// <item><a href="https://ngdc.noaa.gov/IAGA/vmod/igrf11coeffs.txt"></a></item>
    /// <item><a href="https://ngdc.noaa.gov/IAGA/vmod/geomag70_linux.tar.gz"></a></item>
    /// </list>
    /// </item>
    /// <item>
    /// EMM2010:
    /// <list type="bullet">
    /// <item><a href="https://ngdc.noaa.gov/geomag/EMM/index.html"></a></item>
    /// <item><a href="https://ngdc.noaa.gov/geomag/EMM/data/geomag/EMM2010_Sph_Windows_Linux.zip"></a></item>
    /// </list>
    /// </item>
    /// <item>
    /// EMM2015:
    /// <list type="bullet">
    /// <item><a href="https://ngdc.noaa.gov/geomag/EMM/index.html"></a></item>
    /// <item><a href="https://www.ngdc.noaa.gov/geomag/EMM/data/geomag/EMM2015_Sph_Linux.zip"></a></item>
    /// </list>
    /// </item>
    /// <item>
    /// EMM2017:
    /// <list type="bullet">
    /// <item><a href="https://ngdc.noaa.gov/geomag/EMM/index.html"></a></item>
    /// <item><a href="https://www.ngdc.noaa.gov/geomag/EMM/data/geomag/EMM2017_Sph_Linux.zip"></a></item>
    /// </list>
    /// </item>
    /// </list>
    /// </para>
    /// </remarks>
    public class MagneticModel : IEllipsoid
    {
        private const int idlength_ = 8;

        private readonly DateTime? _date;
        private readonly string _name, _dir, _description, _filename, _id;
        private readonly double _t0, _dt0, _tmin, _tmax, _a, _hmin, _hmax;
        private readonly int _Nmodels, _Nconstants, _nmx, _mmx;
        private readonly Normalization _norm;
        private readonly Geocentric _earth;
        private readonly List<SphericalHarmonic> _harm = new List<SphericalHarmonic>();

        static MagneticModel()
        {
            string GetDefaultPath()
            {
                var path = Environment.GetEnvironmentVariable("GEOGRAPHICLIB_MAGNETIC_PATH");
                if (!string.IsNullOrEmpty(path))
                    return path;

                var datapath = Environment.GetEnvironmentVariable("GEOGRAPHICLIB_DATA");
                return Path.Combine(!string.IsNullOrEmpty(datapath) ? datapath : GEOGRAPHICLIB_DATA, "magnetic");
            }

            DefaultMagneticPath = GetDefaultPath();
        }

        /// <summary>
        /// Construct a magnetic model.
        /// </summary>
        /// <param name="name">the name of the model</param>
        /// <param name="path">(optional) directory for data file.</param>
        /// <param name="earth">(optional) <see cref="Geocentric"/> object for converting coordinates; default <see cref="Geocentric.WGS84"/>.</param>
        /// <param name="Nmax">(optional) if non-negative, truncate the degree of the model this value.</param>
        /// <param name="Mmax">(optional) if non-negative, truncate the order of the model this value.</param>
        /// <remarks>
        /// A filename is formed by appending ".wmm" (World Magnetic Model) to the name.
        /// If path is specified (and is non-empty), then the file is loaded from directory, <paramref name="path"/>.
        /// Otherwise the path is given by the <see cref="DefaultMagneticPath"/>.
        /// <para>
        /// This file contains the metadata which specifies the properties of the model.
        /// The coefficients for the spherical harmonic sums are obtained from a file obtained by appending ".cof" to metadata file (so the filename ends in ".wwm.cof").
        /// </para>
        /// <para>
        /// The model is not tied to a particular ellipsoidal model of the earth.
        /// The final earth argument to the constructor specifies an ellipsoid to allow geodetic coordinates to the 
        /// transformed into the spherical coordinates used in the spherical harmonic sum.
        /// </para>
        /// <para>
        /// If <paramref name="Nmax"/> ≥ 0 and <paramref name="Mmax"/> &lt; 0, then <paramref name="Mmax"/> is set to <paramref name="Nmax"/>.
        /// After the model is loaded, the maximum degree and order of the model can be found by the <see cref="Degree"/> and <see cref="Order"/> methods.
        /// </para>
        /// </remarks>
        public MagneticModel(string name, string path = "", Geocentric earth = null, int Nmax = -1, int Mmax = -1)
            : this(earth)
        {
            bool truncate = Nmax >= 0 || Mmax >= 0;
            if (truncate)
            {
                if (Nmax >= 0 && Mmax < 0) Mmax = Nmax;
                if (Nmax < 0) Nmax = int.MaxValue;
                if (Mmax < 0) Mmax = int.MaxValue;
            }

            _dir = string.IsNullOrEmpty(path) ? DefaultMagneticPath : path;
            _filename = _dir + "/" + name + ".wmm";
            _name = name;

            using (var stream = File.OpenRead(_filename))
            {
                ReadMetadata(
                    _filename,
                    stream, ref _id, ref _name, ref _description, ref _date,
                    ref _a, ref _t0, ref _dt0, ref _tmin, ref _tmax, ref _hmin, ref _hmax,
                    ref _Nmodels, ref _Nconstants, ref _norm);
            }

            string coeff = _filename + ".cof";
            using (var stream = File.OpenRead(coeff))
            {
                ReadCoefficients(coeff, stream, truncate, Nmax, Mmax, ref _nmx, ref _mmx);
            }
        }

        /// <summary>
        /// Construct a magnetic model from the given <paramref name="metadataStream"/> and <paramref name="coefficientsStream"/>.
        /// </summary>
        /// <param name="metadataStream">A <see cref="Stream"/> which contains the metadata of the magnetic model.</param>
        /// <param name="coefficientsStream">A <see cref="Stream"/> which contains the coefficients of the magnetic model.</param>
        /// <param name="earth">(optional) <see cref="Geocentric"/> object for converting coordinates; default <see cref="Geocentric.WGS84"/>.</param>
        /// <param name="Nmax">(optional) if non-negative, truncate the degree of the model this value.</param>
        /// <param name="Mmax">(optional) if non-negative, truncate the order of the model this value.</param>
        /// <param name="leaveOpen">
        /// <see langword="true"/> to leave the streams open after the constructor returns; otherwise, <see langword="false"/>.
        /// </param>
        /// <remarks>
        /// <para>
        /// The model is not tied to a particular ellipsoidal model of the earth.
        /// The final earth argument to the constructor specifies an ellipsoid to allow geodetic coordinates to the 
        /// transformed into the spherical coordinates used in the spherical harmonic sum.
        /// </para>
        /// <para>
        /// If <paramref name="Nmax"/> ≥ 0 and <paramref name="Mmax"/> &lt; 0, then <paramref name="Mmax"/> is set to <paramref name="Nmax"/>.
        /// After the model is loaded, the maximum degree and order of the model can be found by the <see cref="Degree"/> and <see cref="Order"/> methods.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"/>
        public MagneticModel(
            Stream metadataStream,
            Stream coefficientsStream,
            Geocentric earth = null,
            int Nmax = -1,
            int Mmax = -1,
            bool leaveOpen = false)
            : this(earth)
        {
            if (metadataStream == null)
            {
                throw new ArgumentNullException(nameof(metadataStream));
            }

            if (coefficientsStream == null)
            {
                throw new ArgumentNullException(nameof(coefficientsStream));
            }

            _dir = null;
            _filename = null;
            _name = null;

            bool truncate = Nmax >= 0 || Mmax >= 0;
            if (truncate)
            {
                if (Nmax >= 0 && Mmax < 0) Mmax = Nmax;
                if (Nmax < 0) Nmax = int.MaxValue;
                if (Mmax < 0) Mmax = int.MaxValue;
            }

            try
            {
                ReadMetadata(
                    null,
                    metadataStream,
                    ref _id, ref _name, ref _description, ref _date,
                    ref _a, ref _t0, ref _dt0, ref _tmin, ref _tmax, ref _hmin, ref _hmax,
                    ref _Nmodels, ref _Nconstants, ref _norm);

                ReadCoefficients(
                    null,
                    coefficientsStream,
                    truncate, Nmax, Mmax, ref _nmx, ref _mmx);
            }
            finally
            {
                if (!leaveOpen)
                {
                    metadataStream.Dispose();
                    coefficientsStream.Dispose();
                }
            }
        }

        /// <summary>
        /// Construct a magnetic model from the given <paramref name="metadataBytes"/> and <paramref name="coefficientsBytes"/>.
        /// </summary>
        /// <param name="metadataBytes">A <see cref="Byte"/> array which contains the metadata of the magnetic model.</param>
        /// <param name="coefficientsBytes">A <see cref="Byte"/> array which contains the coefficients of the magnetic model.</param>
        /// <param name="earth">(optional) <see cref="Geocentric"/> object for converting coordinates; default <see cref="Geocentric.WGS84"/>.</param>
        /// <param name="Nmax">(optional) if non-negative, truncate the degree of the model this value.</param>
        /// <param name="Mmax">(optional) if non-negative, truncate the order of the model this value.</param>
        /// <remarks>
        /// <para>
        /// The model is not tied to a particular ellipsoidal model of the earth.
        /// The final earth argument to the constructor specifies an ellipsoid to allow geodetic coordinates to the 
        /// transformed into the spherical coordinates used in the spherical harmonic sum.
        /// </para>
        /// <para>
        /// If <paramref name="Nmax"/> ≥ 0 and <paramref name="Mmax"/> &lt; 0, then <paramref name="Mmax"/> is set to <paramref name="Nmax"/>.
        /// After the model is loaded, the maximum degree and order of the model can be found by the <see cref="Degree"/> and <see cref="Order"/> methods.
        /// </para>
        /// </remarks>
        public MagneticModel(
            byte[] metadataBytes,
            byte[] coefficientsBytes,
            Geocentric earth = null,
            int Nmax = -1,
            int Mmax = -1)
            : this(
                 new MemoryStream(metadataBytes),
                 new MemoryStream(coefficientsBytes),
                 earth, Nmax, Mmax, leaveOpen: false)
        {

        }

        private MagneticModel(Geocentric earth)
        {
            _description = "NONE";
            _t0 = double.NaN;
            _dt0 = 1;
            _tmin = double.NaN;
            _tmax = double.NaN;
            _a = double.NaN;
            _hmin = double.NaN;
            _hmax = double.NaN;
            _Nmodels = 1;
            _Nconstants = 0;
            _nmx = -1;
            _mmx = -1;
            _norm = Normalization.Schmidt;
            _earth = earth ?? Geocentric.WGS84;
        }

        private void ReadMetadata(
            string name,
            Stream metadataStream, ref string _id, ref string _name, ref string _description, ref DateTime? _date,
            ref double _a, ref double _t0, ref double _dt0, ref double _tmin, ref double _tmax,
            ref double _hmin, ref double _hmax, ref int _Nmodels, ref int _Nconstants, ref Normalization _norm)
        {
            name = name == null ? "metadata stream" : name;
            using (var metastr = new StreamReader(metadataStream, Encoding.UTF8, true, bufferSize: 1024, leaveOpen: true))
            {
                var line = metastr.ReadLine();
                if (!(line.Length >= 6 && line.StartsWith("WMMF-")))
                    throw new GeographicException(name + " does not contain WMMF-n signature");

                var parts = line.TrimEnd().Split(new[] { "WMMF-" }, StringSplitOptions.RemoveEmptyEntries);
                if (parts.Length == 0 || (parts[0] != "1" && parts[0] != "2"))
                    throw new GeographicException("Unknown version in " + name + ": " + parts[0]);

                while ((line = metastr.ReadLine()) != null)
                {
                    if (!Utility.ParseLine(line.AsSpan(), out var key, out var val))
                        continue;
                    // Process key words
                    switch (key)
                    {
                        case "Name":
                            _name = val; break;
                        case "Description":
                            _description = val; break;
                        case "ReleaseDate":
                            _date = val.ParseDateTime(); break;
                        case "Radius":
                            _a = val.ParseDouble(); break;
                        case "Type" when val.ToLower() != "linear":
                            throw new GeographicException("Only linear models are supported");
                        case "Epoch":
                            _t0 = val.ParseDouble(); break;
                        case "DeltaEpoch":
                            _dt0 = val.ParseDouble(); break;
                        case "NumModels":
                            _Nmodels = val.ParseInt32(); break;
                        case "NumConstants":
                            _Nconstants = val.ParseInt32(); break;
                        case "MinTime":
                            _tmin = val.ParseDouble(); break;
                        case "MaxTime":
                            _tmax = val.ParseDouble(); break;
                        case "MinHeight":
                            _hmin = val.ParseDouble(); break;
                        case "MaxHeight":
                            _hmax = val.ParseDouble(); break;
                        case "ID":
                            _id = val; break;
                        case "Normalization":
                            switch (val.ToLower())
                            {
                                case "full": _norm = Normalization.Full; break;
                                case "schmidt": _norm = Normalization.Schmidt; break;
                                default: throw new GeographicException("Unknown normalization " + val);
                            };
                            break;
                        case "ByteOrder":
                            switch (val.ToLower())
                            {
                                case "little": break;
                                case "big": throw new GeographicException("Only little-endian ordering is supported");
                                default: throw new GeographicException("Unknown byte ordering " + val);
                            };
                            break;
                            // else unrecognized keywords are skipped
                    }
                }
                // Check values
                if (!(IsFinite(_a) && _a > 0))
                    throw new GeographicException("Reference radius must be positive");
                if (!(_t0 > 0))
                    throw new GeographicException("Epoch time not defined");
                if (_tmin >= _tmax)
                    throw new GeographicException("Min time exceeds max time");
                if (_hmin >= _hmax)
                    throw new GeographicException("Min height exceeds max height");
                if (_id.Length != idlength_)
                    throw new GeographicException("Invalid ID");
                if (_Nmodels < 1)
                    throw new GeographicException("NumModels must be positive");
                if (!(_Nconstants == 0 || _Nconstants == 1))
                    throw new GeographicException("NumConstants must be 0 or 1");
                if (!(_dt0 > 0))
                {
                    if (_Nmodels > 1)
                        throw new GeographicException("DeltaEpoch must be positive");
                    else
                        _dt0 = 1;
                }
            }
        }

        private void ReadCoefficients(
            string name,
            Stream coefficientsStream,
            bool truncate,
            int Nmax,
            int Mmax,
            ref int _nmx,
            ref int _mmx)
        {
            name = name == null ? "coefficients stream" : name;

            Span<byte> id = stackalloc byte[idlength_];
            if (coefficientsStream.Read(id) != idlength_)
                throw new GeographicException("No header in " + name);
            if (MemoryMarshal.Cast<char, byte>(_id.AsSpan()).SequenceEqual(id))
                throw new GeographicException($"ID mismatch: {_id} vs {Encoding.ASCII.GetString(id.ToArray())}");

            for (int i = 0; i < _Nmodels + 1 + _Nconstants; ++i)
            {
                int N = 0, M = 0;
                if (truncate) { N = Nmax; M = Mmax; }

                var c = SphericalEngine.Coeff.FromStream(coefficientsStream, ref N, ref M, truncate);

                if (!(M < 0 || c.Cv(0) == 0))
                    throw new GeographicException("A degree 0 term is not permitted");

                var sh = new SphericalHarmonic(c, _a, _norm);

                _harm.Add(sh);
                _nmx = Max(_nmx, sh.Coefficients.Nmx);
                _mmx = Max(_mmx, sh.Coefficients.Mmx);
            }

            if (coefficientsStream.Position != coefficientsStream.Length)
                throw new GeographicException("Extra data in " + name);
        }

        private (double Bx, double By, double Bz) Field(double t, double lat, double lon, double h, bool diffp,
                       out double Bxt, out double Byt, out double Bzt)
        {
            Bxt = Byt = Bzt = default;

            Span<double> M = stackalloc double[Geocentric.dim2_];
            var (X, Y, Z) = _earth.IntForward(lat, lon, h, M);
            // Components in geocentric basis
            // initial values to suppress warning
            var (BX, BY, BZ) = FieldGeocentric(t, X, Y, Z, out var BXt, out var BYt, out var BZt);
            if (diffp)
                Geocentric.Unrotate(M, BXt, BYt, BZt, out Bxt, out Byt, out Bzt);
            Geocentric.Unrotate(M, BX, BY, BZ, out var Bx, out var By, out var Bz);

            return (Bx, By, Bz);
        }

        /// <summary>
        /// Compute the magnetic field in geocentric coordinate.
        /// </summary>
        /// <param name="t">the time (fractional years).</param>
        /// <param name="X"><i>X</i> component of geocentric coordinate (meters).</param>
        /// <param name="Y"><i>Y</i> component of geocentric coordinate (meters).</param>
        /// <param name="Z"><i>Z</i> component of geocentric coordinate (meters).</param>
        /// <param name="BXt">the rate of change of <i>BX</i> (nT/yr).</param>
        /// <param name="BYt">the rate of change of <i>BY</i> (nT/yr).</param>
        /// <param name="BZt">the rate of change of <i>BZ</i> (nT/yr).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>BX</i>, the <i>X</i> component of the magnetic field (nanotesla).</item>
        /// <item><i>BY</i>, the <i>Y</i> component of the magnetic field (nanotesla).</item>
        /// <item><i>BZ</i>, the <i>Z</i> component of the magnetic field (nanotesla).</item>
        /// </list>
        /// </returns>
        public (double BX, double BY, double BZ) FieldGeocentric(double t, double X, double Y, double Z,
                              out double BXt, out double BYt, out double BZt)
        {
            t -= _t0;
            int n = Max(Min((int)Floor(t / _dt0), _Nmodels - 1), 0);
            bool interpolate = n + 1 < _Nmodels;
            t -= n * _dt0;
            // Components in geocentric basis
            // initial values to suppress warning
            double BXc = 0, BYc = 0, BZc = 0;
            _harm[n].Evaluate(X, Y, Z, out var BX, out var BY, out var BZ);
            _harm[n + 1].Evaluate(X, Y, Z, out BXt, out BYt, out BZt);
            if (_Nconstants != 0)
                _harm[_Nmodels + 1].Evaluate(X, Y, Z, out BXc, out BYc, out BZc);
            if (interpolate)
            {
                // Convert to a time derivative
                BXt = (BXt - BX) / _dt0;
                BYt = (BYt - BY) / _dt0;
                BZt = (BZt - BZ) / _dt0;
            }
            BX += t * BXt + BXc;
            BY += t * BYt + BYc;
            BZ += t * BZt + BZc;

            BXt = BXt * -_a;
            BYt = BYt * -_a;
            BZt = BZt * -_a;

            BX *= -_a;
            BY *= -_a;
            BZ *= -_a;

            return (BX, BY, BZ);
        }

        /// <summary>
        /// Evaluate the components of the geomagnetic field.
        /// </summary>
        /// <param name="t">the time (fractional years).</param>
        /// <param name="lat">latitude of the point (degrees).</param>
        /// <param name="lon">longitude of the point (degrees).</param>
        /// <param name="h">the height of the point above the ellipsoid (meters).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>Bx</i>, the easterly component of the magnetic field (nanotesla).</item>
        /// <item><i>By</i>, the northerly component of the magnetic field (nanotesla).</item>
        /// <item><i>Bz</i>, the vertical (up) component of the magnetic field (nanotesla).</item>
        /// </list>
        /// </returns>
        public (double Bx, double By, double Bz) Evaluate(double t, double lat, double lon, double h)
            => Field(t, lat, lon, h, false, out _, out _, out _);

        /// <summary>
        /// Evaluate the components of the geomagnetic field and their time derivatives.
        /// </summary>
        /// <param name="t">the time (fractional years).</param>
        /// <param name="lat">latitude of the point (degrees).</param>
        /// <param name="lon">longitude of the point (degrees).</param>
        /// <param name="h">the height of the point above the ellipsoid (meters).</param>
        /// <param name="Bxt">the rate of change of <i>Bx</i> (nT/yr).</param>
        /// <param name="Byt">the rate of change of <i>By</i> (nT/yr).</param>
        /// <param name="Bzt">the rate of change of <i>Bz</i> (nT/yr).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>Bx</i>, the easterly component of the magnetic field (nanotesla).</item>
        /// <item><i>By</i>, the northerly component of the magnetic field (nanotesla).</item>
        /// <item><i>Bz</i>, the vertical (up) component of the magnetic field (nanotesla).</item>
        /// </list>
        /// </returns>
        /// <remarks>
        ///  Use <see cref="Utility.FractionalYear(string)"/> to convert a date of the form yyyy-mm or yyyy-mm-dd into a fractional year.
        /// </remarks>
        public (double Bx, double By, double Bz) Evaluate(double t, double lat, double lon, double h,
                             out double Bxt, out double Byt, out double Bzt)
            => Field(t, lat, lon, h, true, out Bxt, out Byt, out Bzt);

        /// <summary>
        /// Create a <see cref="MagneticCircle"/> object to allow the geomagnetic field at many points with constant
        /// <i>lat</i>, <i>h</i>, and <i>t</i> and varying lon to be computed efficiently.
        /// </summary>
        /// <param name="h">the height of the point above the ellipsoid (meters).</param>
        /// <param name="lat">latitude of the point (degrees).</param>
        /// <param name="t">the time (fractional years).</param>
        /// <returns>
        /// a <see cref="MagneticCircle"/> object whose <see cref="MagneticCircle.Evaluate(double, out double, out double, out double)"/> member function
        /// computes the field at particular values of <i>lon</i>.
        /// </returns>
        /// <remarks>
        /// If the field at several points on a circle of latitude need to be calculated then creating a <see cref="MagneticCircle"/>
        /// and using its member functions will be substantially faster, especially for high-degree models.
        /// </remarks>
        public MagneticCircle Circle(double t, double lat, double h)
        {
            var t1 = t - _t0;
            int n = Max(Min((int)Floor(t1 / _dt0), _Nmodels - 1), 0);
            bool interpolate = n + 1 < _Nmodels;
            t1 -= n * _dt0;
            Span<double> M = stackalloc double[Geocentric.dim2_];
            var (X, Y, Z) = _earth.IntForward(lat, 0, h, M);
            // Y = 0, cphi = M[7], sphi = M[8];

            return (_Nconstants == 0 ?
                   new MagneticCircle(_a, _earth.Flattening, lat, h, t,
                                   M[7], M[8], t1, _dt0, interpolate,
                                   _harm[n].Circle(X, Z, true),
                                   _harm[n + 1].Circle(X, Z, true)) :
                   new MagneticCircle(_a, _earth.Flattening, lat, h, t,
                                   M[7], M[8], t1, _dt0, interpolate,
                                   _harm[n].Circle(X, Z, true),
                                   _harm[n + 1].Circle(X, Z, true),
                                   _harm[_Nmodels + 1].Circle(X, Z, true)));
        }

        /// <summary>
        /// Compute various quantities dependent on the magnetic field.
        /// </summary>
        /// <param name="Bx">the <i>x</i> (easterly) component of the magnetic field (nanotesla).</param>
        /// <param name="By">the <i>y</i> (northerly) component of the magnetic field (nanotesla).</param>
        /// <param name="Bz">the <i>z</i> (vertical, up positive) component of the magnetic field (nanotesla).</param>
        /// <param name="Bxt">the rate of change of <i>Bx</i> (nT/yr).</param>
        /// <param name="Byt">the rate of change of <i>By</i> (nT/yr).</param>
        /// <param name="Bzt">the rate of change of <i>Bz</i> (nT/yr).</param>
        /// <param name="Ht">the rate of change of <i>H</i> (nT/yr).</param>
        /// <param name="Ft">the rate of change of <i>F</i> (nT/yr).</param>
        /// <param name="Dt">the rate of change of <i>D</i> (degrees/yr).</param>
        /// <param name="It">the rate of change of <i>I</i> (degrees/yr).</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>H</i>, the horizontal magnetic field (nT).</item>
        /// <item><i>F</i>, the total magnetic field (nT).</item>
        /// <item><i>D</i>, the declination of the field (degrees east of north).</item>
        /// <item><i>I</i>, the inclination of the field (degrees down from horizontal).</item>
        /// </list>
        /// </returns>
        public static (double H, double F, double D, double I) FieldComponents(double Bx, double By, double Bz,
                                double Bxt, double Byt, double Bzt,
                                out double Ht, out double Ft, out double Dt, out double It)
        {
            var H = Hypot(Bx, By);
            Ht = H != 0 ? (Bx * Bxt + By * Byt) / H : Hypot(Bxt, Byt);
            var D = H != 0 ? Atan2d(Bx, By) : Atan2d(Bxt, Byt);
            Dt = (H != 0 ? (By * Bxt - Bx * Byt) / Sq(H) : 0) / MathEx.Degree;
            var F = Hypot(H, Bz);
            Ft = F != 0 ? (H * Ht + Bz * Bzt) / F : Hypot(Ht, Bzt);
            var I = F != 0 ? Atan2d(-Bz, H) : Atan2d(-Bzt, Ht);
            It = (F != 0 ? (Bz * Ht - H * Bzt) / Sq(F) : 0) / MathEx.Degree;

            return (H, F, D, I);
        }

        /// <summary>
        /// Compute various quantities dependent on the magnetic field and its rate of change.
        /// </summary>
        /// <param name="Bx">the <i>x</i> (easterly) component of the magnetic field (nanotesla).</param>
        /// <param name="By">the <i>y</i> (northerly) component of the magnetic field (nanotesla).</param>
        /// <param name="Bz">the <i>z</i> (vertical, up positive) component of the magnetic field (nanotesla).</param>
        /// <param name="Bxt">the rate of change of <i>Bx</i> (nT/yr).</param>
        /// <param name="Byt">the rate of change of <i>By</i> (nT/yr).</param>
        /// <param name="Bzt">the rate of change of <i>Bz</i> (nT/yr).</param>
        ///  <returns>
        /// <list type="bullet">
        /// <item><i>H</i>, the horizontal magnetic field (nT).</item>
        /// <item><i>F</i>, the total magnetic field (nT).</item>
        /// <item><i>D</i>, the declination of the field (degrees east of north).</item>
        /// <item><i>I</i>, the inclination of the field (degrees down from horizontal).</item>
        /// </list>
        /// </returns>
        public static (double H, double F, double D, double I) FieldComponents(double Bx, double By, double Bz,
                                double Bxt, double Byt, double Bzt)
            => FieldComponents(Bx, By, Bz, Bxt, Byt, Bzt, out _, out _, out _, out _);

        /// <summary>
        /// This is the value of the environment variable GEOGRAPHICLIB_MAGNETIC_PATH, if set;
        /// otherwise, it is $GEOGRAPHICLIB_DATA/magnetic if the environment variable GEOGRAPHICLIB_DATA is set;
        /// otherwise, it is a compile-time default (/usr/local/share/GeographicLib/magnetic on non-Windows systems 
        /// and C:/ProgramData/GeographicLib/magnetic on Windows systems).
        /// </summary>
        public static string DefaultMagneticPath { get; }

        /// <summary>
        /// This is the value of the environment variable GEOGRAPHICLIB_MAGNETIC_NAME, if set; otherwise, it is "wmm2020".
        /// The <see cref="MagneticModel"/> class does not use this function;
        /// it is just provided as a convenience for a calling program when constructing a <see cref="MagneticModel"/> object.
        /// </summary>
        public static string DefaultMagneticName { get; } = Environment.GetEnvironmentVariable("GEOGRAPHICLIB_MAGNETIC_NAME") ?? "wmm2020";

        /// <summary>
        /// Gets a value representing the equatorial radius (<i>a</i>) of the ellipsoid.
        /// </summary>
        public double EquatorialRadius => _earth.EquatorialRadius;

        /// <summary>
        /// Gets a value representing the flattening (<i>f</i>) of the ellipsoid.
        /// </summary>
        public double Flattening => _earth.Flattening;

        /// <summary>
        /// Gets a value representing the maximum degree of the components of the model.
        /// </summary>
        public int Degree => _nmx;

        /// <summary>
        /// Gets a value representing the maximum order of the components of the model.
        /// </summary>
        public int Order => _mmx;

        /// <summary>
        /// Gets a value representing the minimum height above the ellipsoid (in meters) for which this <see cref="MagneticModel"/> should be used.
        /// </summary>
        /// <remarks>
        /// Because the model will typically provide useful results slightly outside the range of allowed heights,
        /// no check of t argument is made by <see cref="Evaluate(double, double, double, double, out double, out double, out double)"/> or 
        /// <see cref="Circle(double, double, double)"/>.
        /// </remarks>
        public double MinHeight => _hmin;

        /// <summary>
        /// Gets a value representing the maximum height above the ellipsoid (in meters) for which this <see cref="MagneticModel"/> should be used.
        /// </summary>
        /// <remarks>
        /// Because the model will typically provide useful results slightly outside the range of allowed heights,
        /// no check of t argument is made by <see cref="Evaluate(double, double, double, double, out double, out double, out double)"/> or 
        /// <see cref="Circle(double, double, double)"/>.
        /// </remarks>
        public double MaxHeight => _hmax;

        /// <summary>
        /// Gets a value representing the minimum time (in years) for which this <see cref="MagneticModel"/> should be used.
        /// </summary>
        /// <remarks>
        /// Because the model will typically provide useful results slightly outside the range of allowed heights,
        /// no check of t argument is made by <see cref="Evaluate(double, double, double, double, out double, out double, out double)"/> or 
        /// <see cref="Circle(double, double, double)"/>.
        /// </remarks>
        public double MinTime => _tmin;

        /// <summary>
        /// Gets a value representing the maximum time (in years) for which this <see cref="MagneticModel"/> should be used.
        /// </summary>
        /// <remarks>
        /// Because the model will typically provide useful results slightly outside the range of allowed heights,
        /// no check of t argument is made by <see cref="Evaluate(double, double, double, double, out double, out double, out double)"/> or 
        /// <see cref="Circle(double, double, double)"/>.
        /// </remarks>
        public double MaxTime => _tmax;

        /// <summary>
        /// Gets a value representing the description of the magnetic model, if available, from the Description file in the data file;
        /// if absent, return "NONE".
        /// </summary>
        public string Description => _description;

        /// <summary>
        /// Gets a value representing release date of the model, if available, from the ReleaseDate field in the data file;
        /// if absent, return <see langword="null"/>.
        /// </summary>
        public DateTime? DateTime => _date;

        /// <summary>
        /// Gets a value representing the full file name used to load the magnetic model.
        /// </summary>
        /// <remarks>
        /// This property returns <see langword="null"/> if the <see cref="MagneticModel"/> object
        /// is constructed from <see cref="Stream"/> or <see cref="Byte"/> array.
        /// </remarks>
        public string MagneticFile => _filename;

        /// <summary>
        /// Gets a value representing the "name" used to load the magnetic model
        /// (from the first argument of the constructor, but this may be overridden by the model file).
        /// </summary>
        /// <remarks>
        /// This property returns <see langword="null"/> if the <see cref="MagneticModel"/> object
        /// is constructed from <see cref="Stream"/> or <see cref="Byte"/> array, and the model file doesn't
        /// define <c>Name</c> attribute.
        /// </remarks>
        public string MagneticModelName => _name;

        /// <summary>
        /// Gets a value representing the directory used to load the magnetic model.
        /// </summary>
        /// <remarks>
        /// This property returns <see langword="null"/> if the <see cref="MagneticModel"/> object
        /// is constructed from <see cref="Stream"/> or <see cref="Byte"/> array.
        /// </remarks>
        public string MagneticModelDirectory => _dir;
    }
}
