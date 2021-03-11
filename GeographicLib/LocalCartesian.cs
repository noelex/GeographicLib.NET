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
    /// Local cartesian coordinates.
    /// </summary>
    /// <remarks>
    /// Convert between geodetic coordinates latitude = <i>lat</i>, longitude = <i>lon</i>, height = <i>h</i>
    /// (measured vertically from the surface of the ellipsoid) to local cartesian coordinates (<i>x</i>, <i>y</i>, <i>z</i>).
    /// The origin of local cartesian coordinate system is at <i>lat</i> = <i>lat0</i>, <i>lon</i> = <i>lon0</i>, <i>h</i> = <i>h0</i>.
    /// The <i>z</i> axis is normal to the ellipsoid; the <i>y</i> axis points due north. The plane <i>z</i> = - <i>h0</i> is tangent to the ellipsoid.
    /// <para>
    /// The conversions all take place via geocentric coordinates using a <see cref="Geocentric"/> object (by default <see cref="Geocentric.WGS84"/>).
    /// </para>
    /// </remarks>
    public class LocalCartesian : IEllipsoid
    {
        private const int
            dim_ = 3,
            dim2_ = dim_ * dim_;

        private readonly Geocentric _earth;

        private double _lat0, _lon0, _h0;
        private double _x0, _y0, _z0;

        private readonly Memory<double> _r = new double[dim2_];

        /// <summary>
        /// Initialize a new <see cref="LocalCartesian"/> instance with specified origin.
        /// </summary>
        /// <param name="lat0">latitude at origin (degrees).</param>
        /// <param name="lon0">longitude at origin (degrees).</param>
        /// <param name="h0">height above ellipsoid at origin (meters); default 0.</param>
        /// <param name="earth"><see cref="Geocentric"/> object for the transformation; default <see cref="Geocentric.WGS84"/>.</param>
        /// <remarks>
        /// <paramref name="lat0"/> should be in the range [−90°, 90°].
        /// </remarks>
        public LocalCartesian(double lat0, double lon0, double h0 = 0, Geocentric earth = null)
        {
            _earth = earth ?? Geocentric.WGS84;
            Reset(lat0, lon0, h0);
        }

        /// <summary>
        /// Initialize a new <see cref="LocalCartesian"/> instance.
        /// </summary>
        /// <param name="earth"><see cref="Geocentric"/> object for the transformation; default <see cref="Geocentric.WGS84"/>.</param>
        public LocalCartesian(Geocentric earth = null) : this(0, 0, earth: earth) { }

        private (double x, double y, double z) IntForward(double lat, double lon, double h, Span<double> M)
        {
            var (xc, yc, zc) = _earth.IntForward(lat, lon, h, M);
            xc -= _x0; yc -= _y0; zc -= _z0;
            var x = _r.Span[0] * xc + _r.Span[3] * yc + _r.Span[6] * zc;
            var y = _r.Span[1] * xc + _r.Span[4] * yc + _r.Span[7] * zc;
            var z = _r.Span[2] * xc + _r.Span[5] * yc + _r.Span[8] * zc;
            if (M != default)
                MatrixMultiply(M);

            return (x, y, z);
        }

        private (double lat, double lon, double h) IntReverse(double x, double y, double z, Span<double> M)
        {
            double
                xc = _x0 + _r.Span[0] * x + _r.Span[1] * y + _r.Span[2] * z,
                yc = _y0 + _r.Span[3] * x + _r.Span[4] * y + _r.Span[5] * z,
                zc = _z0 + _r.Span[6] * x + _r.Span[7] * y + _r.Span[8] * z;

            var (lat, lon, h) = _earth.IntReverse(xc, yc, zc, M);
            if (M!=default)
                MatrixMultiply(M);

            return (lat, lon, h);
        }

        private void MatrixMultiply(Span<double> M)
        {
            Span<double> t = stackalloc double[dim2_];

            M.CopyTo(t);

            for (var i = 0; i < dim2_; ++i)
            {
                int row = i / dim_, col = i % dim_;
                M[i] = _r.Span[row] * t[col] + _r.Span[row + 3] * t[col + 3] + _r.Span[row + 6] * t[col + 6];
            }
        }

        /// <summary>
        /// Reset the origin.
        /// </summary>
        /// <param name="lat0">latitude at origin (degrees).</param>
        /// <param name="lon0">longitude at origin (degrees).</param>
        /// <param name="h0">height above ellipsoid at origin (meters); default 0.</param>
        /// <remarks>
        /// <paramref name="lat0"/> should be in the range [−90°, 90°].
        /// </remarks>
        public void Reset(double lat0, double lon0, double h0 = 0)
        {
            _lat0 = LatFix(lat0);
            _lon0 = AngNormalize(lon0);
            _h0 = h0;
            (_x0, _y0, _z0) = _earth.Forward(_lat0, _lon0, _h0);

            SinCosd(_lat0, out var sphi, out var cphi);
            SinCosd(_lon0, out var slam, out var clam);

            Geocentric.Rotation(sphi, cphi, slam, clam, _r.Span);
        }

        /// <summary>
        /// Convert from geodetic to local cartesian coordinates and return rotation matrix.
        /// </summary>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <param name="h">height of point above the ellipsoid (meters).</param>
        /// <param name="M">if the length of the vector is 9, fill with the rotation matrix in row-major order.</param>
        /// <returns>
        /// local cartesian coordinate (meters).
        /// </returns>
        /// <remarks>
        /// <paramref name="lat"/> should be in the range [−90°, 90°].
        /// <para>
        /// Let <i>v</i> be a unit vector located at (<paramref name="lat"/>, <paramref name="lon"/>, <paramref name="h"/>).
        /// We can express <i>v</i> as column vectors in one of two ways
        /// <list type="bullet">
        /// <item>
        /// in east, north, up coordinates (where the components are relative to a local coordinate system at 
        /// (<paramref name="lat"/>, <paramref name="lon"/>, <paramref name="h"/>)); call this representation <i>v1</i>.
        /// </item>
        /// <item>
        /// in <i>x</i>, <i>y</i>, <i>z</i> coordinates (where the components are relative to the local coordinate system at
        /// (<i>lat0</i>, <i>lon0</i>, <i>h0</i>)); call this representation <i>v0</i>.
        /// </item>
        /// </list>
        /// Then we have <i>v0</i> = <paramref name="M"/> ⋅ <i>v1</i>.
        /// </para>
        /// </remarks>
        public (double x, double y, double z) Forward(double lat, double lon, double h, Span<double> M = default)
            => M.Length == dim2_ ? IntForward(lat, lon, h, M) : IntForward(lat, lon, h, default);

        /// <summary>
        /// Convert from local cartesian to geodetic coordinates.
        /// </summary>
        /// <param name="x"><i>x</i> component of local cartesian coordinate (meters).</param>
        /// <param name="y"><i>y</i> component of local cartesian coordinate (meters).</param>
        /// <param name="z"><i>z</i> component of local cartesian coordinate (meters).</param>
        /// <param name="M">if the length of the vector is 9, fill with the rotation matrix in row-major order.</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>lat</i>, latitude of point (degrees).</item>
        /// <item><i>lon</i>, longitude of point (degrees).</item>
        /// <item><i>h</i>, height of point above the ellipsoid (meters).</item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// In general, there are multiple solutions and the result which minimizes |<i>h</i> |is returned, i.e., (<i>lat</i>, <i>lon</i>) corresponds
        /// to the closest point on the ellipsoid. The value of lon returned is in the range [−180°, 180°].
        /// <para>
        /// Let <i>v</i> be a unit vector located at (<i>lat</i>, <i>lon</i>, <i>h</i>).
        /// We can express <i>v</i> as column vectors in one of two ways
        /// <list type="bullet">
        /// <item>
        /// in east, north, up coordinates (where the components are relative to a local coordinate system at 
        /// (<i>lat</i>, <i>lon</i>, <i>h</i>)); call this representation <i>v1</i>.
        /// </item>
        /// <item>
        /// in <i>x</i>, <i>y</i>, <i>z</i> coordinates (where the components are relative to the local coordinate system at
        /// (<i>lat0</i>, <i>lon0</i>, <i>h0</i>)); call this representation <i>v0</i>.
        /// </item>
        /// </list>
        /// Then we have <i>v1</i> = <paramref name="M"/>t ⋅ <i>v0</i>, where <paramref name="M"/>t is the transpose of <paramref name="M"/>.
        /// </para>
        /// </remarks>
        public (double lat, double lon, double h) Reverse(double x, double y, double z, Span<double> M = default)
            => M.Length == dim2_ ? IntReverse(x, y, z, M) : IntReverse(x, y, z, default);

        /// <summary>
        /// Gets a values representing latitude of the origin in degrees.
        /// </summary>
        public double LatitudeOrigin => _lat0;

        /// <summary>
        /// Gets a values representing longitude of the origin in degrees.
        /// </summary>
        public double LongitudeOrigin => _lon0;

        /// <summary>
        /// Gets a values representing height of the origin in meters.
        /// </summary>
        public double HeightOrigin => _h0;

        /// <inheritdoc/>
        public double EquatorialRadius => _earth.EquatorialRadius;

        /// <inheritdoc/>
        public double Flattening => _earth.Flattening;
    }
}
