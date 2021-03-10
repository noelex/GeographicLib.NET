using System;
using System.Collections.Generic;
using System.Text;

using static System.Math;
using static GeographicLib.MathEx;

namespace GeographicLib
{
    /// <summary>
    /// Provides conversion between geodetic and geocentric coordinates.
    /// </summary>
    /// <remarks>
    /// Convert between geodetic coordinates latitude = <i>lat</i>, longitude = <i>lon</i>, height = <i>h</i>
    /// (measured vertically from the surface of the ellipsoid)  to geocentric coordinates (<i>X</i>, <i>Y</i>, <i>Z</i>). 
    /// The origin of geocentric coordinates is at the center of the earth. 
    /// The <i>Z</i> axis goes thru the north pole, <i>lat</i> = 90°. 
    /// The <i>X</i> axis goes thru <i>lat</i> = 0, <i>lon</i> = 0. 
    /// Geocentric coordinates are also known as earth centered, earth fixed (ECEF) coordinates.
    /// <para>
    /// The conversion from geographic to geocentric coordinates is straightforward. For the reverse transformation we use
    /// <list type="bullet">
    /// <item>
    /// H. Vermeille, 
    /// <a href="https://doi.org/10.1007/s00190-002-0273-6">Direct transformation from geocentric coordinates to geodetic coordinates</a>,
    /// J. Geodesy 76, 451–454 (2002).
    /// </item>
    /// </list>
    /// </para>
    /// <para>
    /// Several changes have been made to ensure that the method returns accurate results for all finite inputs (even if <i>h</i> is infinite).
    /// The changes are described in Appendix B of
    /// <list type="bullet">
    /// <item>
    /// C. F. F. Karney,
    /// <a href="https://arxiv.org/abs/1102.1215v1">Geodesics on an ellipsoid of revolution</a>,
    /// Feb. 2011; preprint <a href="https://arxiv.org/abs/1102.1215v1">arxiv:1102.1215v1</a>.
    /// </item>
    /// </list>
    /// </para>
    /// <para>
    /// Vermeille similarly updated his method in
    /// <list type="bullet">
    /// <item>
    /// H. Vermeille,
    /// <a href="https://doi.org/10.1007/s00190-010-0419-x">An analytical method to transform geocentric into geodetic coordinates</a>,
    /// J. Geodesy 85, 105–117 (2011).
    /// </item>
    /// </list>
    /// </para>
    /// <para>
    /// See <a href="https://geographiclib.sourceforge.io/html/geocentric.html">Geocentric coordinates</a> for more information.
    /// </para>
    /// <para>
    /// The errors in these routines are close to round-off.
    /// Specifically, for points within 5000 km of the surface of the ellipsoid (either inside or outside the ellipsoid),
    /// the error is bounded by 7 nm (7 nanometers) for the WGS84 ellipsoid.
    /// See <a href="https://geographiclib.sourceforge.io/html/geocentric.html">Geocentric coordinates</a> for further information on the errors.
    /// </para>
    /// </remarks>
    public class Geocentric : IEllipsoid
    {
        private const int dim_ = 3;

        internal const int dim2_ = dim_ * dim_;

        private readonly double _a, _f, _e2, _e2m, _e2a, _e4a, _maxrad;

        /// <summary>
        /// Initialize a new <see cref="Geocentric"/> instance with specified equatorial radius and flattening of an ellipsoid.
        /// </summary>
        /// <param name="a">equatorial radius (meters).</param>
        /// <param name="f">
        /// flattening of ellipsoid. Setting <paramref name="f"/> = 0 gives a sphere.
        /// Negative <paramref name="f"/> gives a prolate ellipsoid.
        /// </param>
        public Geocentric(double a, double f)
        {
            _a = a;
            _f = f;
            _e2 = _f * (2 - _f);
            _e2m = Sq(1 - _f);
            _e2a = Abs(_e2);
            _e4a = Sq(_e2);
            _maxrad = 2 * _a / DBL_EPSILON;

            if (!(IsFinite(_a) && _a > 0))
                throw new GeographicException("Equatorial radius is not positive");
            if (!(IsFinite(_f) && _f < 1))
                throw new GeographicException("Polar semi-axis is not positive");
        }

        /// <summary>
        /// Initialize a new <see cref="Geocentric"/> instance with specified <see cref="IEllipsoid"/> instance.
        /// </summary>
        /// <param name="ellipsoid">The <see cref="IEllipsoid"/> to use.</param>
        public Geocentric(IEllipsoid ellipsoid) : this(ellipsoid.EquatorialRadius, ellipsoid.Flattening) { }

        /// <summary>
        /// Convert from geodetic to geocentric coordinates and return rotation matrix.
        /// </summary>
        /// <param name="lat">latitude of point (degrees).</param>
        /// <param name="lon">longitude of point (degrees).</param>
        /// <param name="h">height of point above the ellipsoid (meters).</param>
        /// <param name="M">if the length of the vector is 9, fill with the rotation matrix in row-major order.</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>X</i>, <i>x</i> component of geocentric coordinate (meters).</item>
        /// <item><i>Y</i>, <i>y</i> component of geocentric coordinate (meters).</item>
        /// <item><i>Z</i>, <i>z</i> component of geocentric coordinate (meters).</item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// Let <i>v</i> be a unit vector located at (<i>lat</i>, <i>lon</i>, <i>h</i>).
        /// We can express <i>v</i> as column vectors in one of two ways
        /// <list type="bullet">
        /// <item>
        /// in east, north, up coordinates (where the components are relative to a local coordinate system at (<i>lat</i>, <i>lon</i>, <i>h</i>));
        /// call this representation <i>v1</i>.
        /// </item>
        /// <item>
        /// in geocentric <i>X</i>, <i>Y</i>, <i>Z</i> coordinates; call this representation <i>v0</i>.
        /// </item>
        /// </list>
        /// <para>
        /// Then we have <i>v0</i> = <i>M</i> ⋅ <i>v1</i>.
        /// </para>
        /// </remarks>
        public (double X, double Y, double Z) Forward(double lat, double lon, double h, Span<double> M = default)
            => IntForward(lat, lon, h, M.Length == dim2_ ? M : default);

        /// <summary>
        /// Convert from geocentric to geodetic to coordinates.
        /// </summary>
        /// <param name="X"><i>x</i> component of geocentric coordinate (meters).</param>
        /// <param name="Y"><i>y</i> component of geocentric coordinate (meters).</param>
        /// <param name="Z"><i>z</i> component of geocentric coordinate (meters).</param>
        /// <param name="M"> the length of the vector is 9, fill with the rotation matrix in row-major order.</param>
        /// <returns>
        /// <list type="bullet">
        /// <item><i>X</i>, <i>x</i> component of geocentric coordinate (meters).</item>
        /// <item><i>Y</i>, <i>y</i> component of geocentric coordinate (meters).</item>
        /// <item><i>Z</i>, <i>z</i> component of geocentric coordinate (meters).</item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// Let <i>v</i> be a unit vector located at (<i>lat</i>, <i>lon</i>, <i>h</i>).
        /// We can express <i>v</i> as column vectors in one of two ways
        /// <list type="bullet">
        /// <item>
        /// in east, north, up coordinates (where the components are relative to a local coordinate system at (<i>lat</i>, <i>lon</i>, <i>h</i>));
        /// call this representation <i>v1</i>.
        /// </item>
        /// <item>
        /// in geocentric <i>X</i>, <i>Y</i>, <i>Z</i> coordinates; call this representation <i>v0</i>.
        /// </item>
        /// </list>
        /// <para>
        /// Then we have <i>v1</i> = <i>M</i>^T ⋅ <i>v0</i>, where <i>M</i>^T is the transpose of <i>M</i>.
        /// </para>
        /// <para>
        /// In general, there are multiple solutions and the result which minimizes |<i>h</i>|is returned, 
        /// i.e., (<i>lat</i>, <i>lon</i>) corresponds to the closest point on the ellipsoid.
        /// If there are still multiple solutions with different latitudes (applies only if <i>Z</i> = 0), 
        /// then the solution with <i>lat</i> > 0 is returned.
        /// If there are still multiple solutions with different longitudes (applies only if <i>X</i> = <i>Y</i> = 0)
        /// then <i>lon</i> = <c>0</c> is returned. The value of h returned satisfies <i>h</i> ≥ − <i>a</i> (1 − <i>e</i>^2) / sqrt(1 − <i>e</i>^2 sin^2<i>lat</i>).
        /// The value of lon returned is in the range [−180°, 180°].
        /// </para>
        /// </remarks>
        public (double lat, double lon, double h) Reverse(double X, double Y, double Z, Span<double> M=default)
            => IntReverse(X, Y, Z, M.Length == dim2_ ? M : default);

        /// <summary>
        /// Perform [X,Y,Z]^t = M.[x,y,z]^t (typically local cartesian to geocentric)
        /// </summary>
        /// <param name="M"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        /// <param name="X"></param>
        /// <param name="Y"></param>
        /// <param name="Z"></param>
        internal static void Rotate(ReadOnlySpan<double> M,
                       double x, double y, double z,
                       out double X, out double Y, out double Z)
        {
            X = M[0] * x + M[1] * y + M[2] * z;
            Y = M[3] * x + M[4] * y + M[5] * z;
            Z = M[6] * x + M[7] * y + M[8] * z;
        }

        /// <summary>
        /// Perform [x,y,z]^t = M^t.[X,Y,Z]^t (typically geocentric to local cartesian)
        /// </summary>
        /// <param name="M"></param>
        /// <param name="X"></param>
        /// <param name="Y"></param>
        /// <param name="Z"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        internal static void Unrotate(ReadOnlySpan<double> M,
                               double X, double Y, double Z,
                               out double x, out double y, out double z)
        {
            x = M[0] * X + M[3] * Y + M[6] * Z;
            y = M[1] * X + M[4] * Y + M[7] * Z;
            z = M[2] * X + M[5] * Y + M[8] * Z;
        }

        /// <summary>
        /// This rotation matrix is given by the following quaternion operations
        /// qrot(lam, [0,0,1]) * qrot(phi, [0,-1,0]) * [1,1,1,1]/2
        /// or
        /// qrot(pi/2 + lam, [0,0,1]) * qrot(-pi/2 + phi , [-1,0,0])
        /// where
        /// qrot(t,v) = [cos(t/2), sin(t/2)*v[1], sin(t/2)*v[2], sin(t/2)*v[3]]
        /// </summary>
        /// <param name="sphi"></param>
        /// <param name="cphi"></param>
        /// <param name="slam"></param>
        /// <param name="clam"></param>
        /// <param name="M"></param>
        internal static void Rotation(double sphi, double cphi, double slam, double clam, Span<double> M)
        {
            // Local X axis (east) in geocentric coords
            M[0] = -slam; M[3] = clam; M[6] = 0;
            // Local Y axis (north) in geocentric coords
            M[1] = -clam * sphi; M[4] = -slam * sphi; M[7] = cphi;
            // Local Z axis (up) in geocentric coords
            M[2] = clam * cphi; M[5] = slam * cphi; M[8] = sphi;
        }

        internal (double X, double Y, double Z) IntForward(double lat, double lon, double h, Span<double> M)
        {
            SinCosd(LatFix(lat), out var sphi, out var cphi);
            SinCosd(lon, out var slam, out var clam);
            var n = _a / Sqrt(1 - _e2 * Sq(sphi));
            var Z = (_e2m * n + h) * sphi;
            var X = (n + h) * cphi;
            var Y = X * slam;
            X *= clam;

            if (M != default) Rotation(sphi, cphi, slam, clam, M);

            return (X, Y, Z);
        }

        internal (double lat, double lon, double h) IntReverse(double X, double Y, double Z, Span<double> M)
        {
            double
              R = Hypot(X, Y),
              slam = R != 0 ? Y / R : 0,
              clam = R != 0 ? X / R : 1;
            var h = Hypot(R, Z);      // Distance to center of earth
            double sphi, cphi;
            if (h > _maxrad)
            {
                // We really far away (> 12 million light years); treat the earth as a
                // point and h, above, is an acceptable approximation to the height.
                // This avoids overflow, e.g., in the computation of disc below.  It's
                // possible that h has overflowed to inf; but that's OK.
                //
                // Treat the case X, Y finite, but R overflows to +inf by scaling by 2.
                R = Hypot(X / 2, Y / 2);
                slam = R != 0 ? (Y / 2) / R : 0;
                clam = R != 0 ? (X / 2) / R : 1;
                var H = Hypot(Z / 2, R);
                sphi = (Z / 2) / H;
                cphi = R / H;
            }
            else if (_e4a == 0)
            {
                // Treat the spherical case.  Dealing with underflow in the general case
                // with _e2 = 0 is difficult.  Origin maps to N pole same as with
                // ellipsoid.
                var H = Hypot(h == 0 ? 1 : Z, R);
                sphi = (h == 0 ? 1 : Z) / H;
                cphi = R / H;
                h -= _a;
            }
            else
            {
                // Treat prolate spheroids by swapping R and Z here and by switching
                // the arguments to phi = atan2(...) at the end.
                double
                  p = Sq(R / _a),
                  q = _e2m * Sq(Z / _a),
                  r = (p + q - _e4a) / 6;
                if (_f < 0) Swap(ref p, ref q);
                if (!(_e4a * q == 0 && r <= 0))
                {
                    double
                      // Avoid possible division by zero when r = 0 by multiplying
                      // equations for s and t by r^3 and r, resp.
                      S = _e4a * p * q / 4, // S = r^3 * s
                      r2 = Sq(r),
                      r3 = r * r2,
                      disc = S * (2 * r3 + S);
                    var u = r;
                    if (disc >= 0)
                    {
                        var T3 = S + r3;
                        // Pick the sign on the sqrt to maximize abs(T3).  This minimizes
                        // loss of precision due to cancellation.  The result is unchanged
                        // because of the way the T is used in definition of u.
                        T3 += T3 < 0 ? -Sqrt(disc) : Sqrt(disc); // T3 = (r * t)^3
                                                                 // N.B. cbrt always returns the real root.  cbrt(-8) = -2.
                        var T = Cbrt(T3); // T = r * t
                                          // T can be zero; but then r2 / T -> 0.
                        u += T + (T != 0 ? r2 / T : 0);
                    }
                    else
                    {
                        // T is complex, but the way u is defined the result is real.
                        var ang = Atan2(Sqrt(-disc), -(S + r3));
                        // There are three possible cube roots.  We choose the root which
                        // avoids cancellation.  Note that disc < 0 implies that r < 0.
                        u += 2 * r * Cos(ang / 3);
                    }
                    double
                      v = Sqrt(Sq(u) + _e4a * q), // guaranteed positive
                                                  // Avoid loss of accuracy when u < 0.  Underflow doesn't occur in
                                                  // e4 * q / (v - u) because u ~ e^4 when q is small and u < 0.
                      uv = u < 0 ? _e4a * q / (v - u) : u + v, // u+v, guaranteed positive
                                                               // Need to guard against w going negative due to roundoff in uv - q.
                      w = Max(0, _e2a * (uv - q) / (2 * v)),
                      // Rearrange expression for k to avoid loss of accuracy due to
                      // subtraction.  Division by 0 not possible because uv > 0, w >= 0.
                      k = uv / (Sqrt(uv + Sq(w)) + w),
                      k1 = _f >= 0 ? k : k - _e2,
                      k2 = _f >= 0 ? k + _e2 : k,
                      d = k1 * R / k2,
                      H = Hypot(Z / k1, R / k2);
                    sphi = (Z / k1) / H;
                    cphi = (R / k2) / H;
                    h = (1 - _e2m / k1) * Hypot(d, Z);
                }
                else
                {                  // e4 * q == 0 && r <= 0
                                   // This leads to k = 0 (oblate, equatorial plane) and k + e^2 = 0
                                   // (prolate, rotation axis) and the generation of 0/0 in the general
                                   // formulas for phi and h.  using the general formula and division by 0
                                   // in formula for h.  So handle this case by taking the limits:
                                   // f > 0: z -> 0, k      ->   e2 * sqrt(q)/sqrt(e4 - p)
                                   // f < 0: R -> 0, k + e2 -> - e2 * sqrt(q)/sqrt(e4 - p)
                    double
                      zz = Sqrt((_f >= 0 ? _e4a - p : p) / _e2m),
                      xx = Sqrt(_f < 0 ? _e4a - p : p),
                      H = Hypot(zz, xx);
                    sphi = zz / H;
                    cphi = xx / H;
                    if (Z < 0) sphi = -sphi; // for tiny negative Z (not for prolate)
                    h = -_a * (_f >= 0 ? _e2m : 1) * H / _e2a;
                }
            }
            var lat = Atan2d(sphi, cphi);
            var lon = Atan2d(slam, clam);

            if (M != null) Rotation(sphi, cphi, slam, clam, M);

            return (lat, lon, h);
        }

        /// <summary>
        /// Gets a value representing the equatorial radius (<i>a</i>) of the ellipsoid.
        /// </summary>
        public double EquatorialRadius => _a;

        /// <summary>
        /// Gets a value representing the flattening (<i>f</i>) of the ellipsoid.
        /// </summary>
        public double Flattening => _f;

        /// <summary>
        /// A global instantiation of <see cref="Geocentric"/> with the parameters for the WGS84 ellipsoid.
        /// </summary>
        public static Geocentric WGS84 { get; } = new Geocentric(Ellipsoid.WGS84);
    }
}
