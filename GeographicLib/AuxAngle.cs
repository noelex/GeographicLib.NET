using System.Runtime.CompilerServices;
using static GeographicLib.MathEx;
using static System.Math;

namespace GeographicLib
{
    /// <summary>
    /// An accurate representation of angles.
    /// </summary>
    /// <remarks>
    /// This class is an implementation of the methods described in
    /// 
    /// <list type="bullet">
    /// <item>
    /// C. F. F. Karney,
    /// <a href = "https://doi.org/10.1080/00396265.2023.2217604">On auxiliary latitudes</a>,
    /// Survey Review(2023);
    /// preprint <a href = "https://arxiv.org/abs/2212.05818">arXiv:2212.05818</a>.
    /// </item>
    /// </list>
    /// 
    /// An angle is represented by the <i>y</i> and <i>x</i> coordinates of a point in the
    /// 2d plane.The two coordinates are proportional to the sine and cosine of
    /// the angle.This allows angles close to the cardinal points to be
    /// represented accurately.It also saves on unnecessary recomputations of
    /// trigonometric functions of the angle.Only angles in [-180°,180°] can be represented.
    /// (A possible extension would be to keep count of the number of turns.)
    /// </remarks>
    public readonly struct AuxAngle
    {
        /// <summary>
        /// A "NaN" <see cref="AuxAngle"/>.
        /// </summary>
        public static readonly AuxAngle NaN = new AuxAngle(double.NaN, double.NaN);

        private readonly double _y, _x;

        /// <summary>
        /// Initializes a new instance of the <see cref="AuxAngle"/> struct with the specified coordinates. 
        /// </summary>
        /// <param name="y">The <i>y</i> coordinate.</param>
        /// <param name="x">The <i>x</i> coordinate.</param>
        public AuxAngle(double y = 0, double x = 1)
        {
            _x = x; _y = y;
        }

        /// <summary>
        /// Gets the <i>x</i> coordinate of the <see cref="AuxAngle"/>.
        /// </summary>
        public double X
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => _x;
        }

        /// <summary>
        /// Gets the <i>y</i> coordinate of the <see cref="AuxAngle"/>.
        /// </summary>
        public double Y
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => _y;
        }

        /// <summary>
        /// Gets the tangent of the angle.
        /// </summary>
        public double Tan
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => _y / _x;
        }

        /// <summary>
        /// Gets the <see cref="AuxAngle"/> converted to the conventional angle measured in degrees.
        /// </summary>
        public double Degrees
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => Atan2d(_y, _x);
        }

        /// <summary>
        /// Gets the <see cref="AuxAngle"/> converted to the conventional angle measured in radians.
        /// </summary>
        public double Radians
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => Atan2(_y, _x);
        }

        /// <summary>
        /// Gets the lambertian of the <see cref="AuxAngle"/>.
        /// </summary>
        public double Lam
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => Asinh(Tan);
        }

        /// <summary>
        /// Gets the lambertian of the <see cref="AuxAngle"/> in degrees.
        /// </summary>
        public double Lamd
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => Lam / Degree;
        }

        /// <summary>
        /// Gets a new normalized <see cref="AuxAngle"/> with the point lying on the unit
        /// circle and the <i>y</i> and <i>x</i> components are equal to the sine and
        /// cosine of the angle.
        /// </summary>
        /// <returns>A normalized copy of current <see cref="AuxAngle"/> instance.</returns>
        public AuxAngle Normalized()
        {
            if (double.IsNaN(Tan) ||
               (Abs(_y) > double.MaxValue / 2 &&
                Abs(_x) > double.MaxValue / 2))
            {
                // deal with
                // (0,0), (inf,inf), (nan,nan), (nan,x), (y,nan), (toobig,toobig)
                return NaN;
            }

            double r = Hypot(_y, _x),
              y = _y / r, x = _x / r;

            // deal with r = inf, then one of y,x becomes 1
            if (double.IsNaN(y)) y = CopySign(1, _y);
            if (double.IsNaN(x)) x = CopySign(1, _x);
            return new AuxAngle(y, x);
        }

        /// <summary>
        /// Copy the quadrant from specified <see cref="AuxAngle"/> instance.
        /// </summary>
        /// <param name="p">The <see cref="AuxAngle"/> instance from which the quadrant information is taken.</param>
        /// <returns>A copy of current <see cref="AuxAngle"/> instance in the same quadrant as <paramref name="p"/>.</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public AuxAngle CopyQuadrant(AuxAngle p)
            => new AuxAngle(CopySign(_y, p._y), CopySign(_x, p._x));

        /// <summary>
        /// Compute the absolute error in another angle.
        /// </summary>
        /// <param name="p">The other angle.</param>
        /// <returns>The absolute error between <paramref name="p"/> and <see langword="this"/> considered as angles in radians.</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double AbsError(AuxAngle p)
           => Abs((new AuxAngle(-p._y, p._x) + this).Radians);

        /// <summary>
        /// Compute the relative error in another angle.
        /// </summary>
        /// <param name="p">The other angle.</param>
        /// <returns>The relative error between <paramref name="p"/>.Tan  and <see langword="this"/>.Tan.</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double RelError(AuxAngle p)
           => Abs(((p._y / p._x) - Tan) / Tan);

        /// <summary>
        /// Create a copy of current <see cref="AuxAngle"/> instance with <see cref="X"/> or <see cref="Y"/> replaced with given value.
        /// </summary>
        /// <param name="y">The <i>y</i> coordinate.</param>
        /// <param name="x">The <i>x</i> coordinate.</param>
        /// <returns>The copied <see cref="AuxAngle"/> instance.</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public AuxAngle CopyWith(double? y = null, double? x = null)
        {
            return new AuxAngle(y ?? Y, x ?? X);
        }

        /// <summary>
        /// Add two <see cref="AuxAngle"/> instances.
        /// </summary>
        /// <remarks>
        /// Neither <paramref name="a"/> nor <paramref name="b"/> should have an infinite component.
        /// If necessary, create normalized angles with <see cref="Normalized"/> first.
        /// </remarks>
        /// <param name="a">The <see cref="AuxAngle"/> to be added.</param>
        /// <param name="b">The <see cref="AuxAngle"/> to be added.</param>
        /// <returns>The sum of <paramref name="a"/> and <paramref name="b"/>.</returns>
        public static AuxAngle operator +(AuxAngle a, AuxAngle b)
        {
            // Do nothing if p.tan() == 0 to preserve signs of y() and x()
            if (b.Tan != 0)
            {
                double
                    x = a._x * b._x - a._y * b._y,
                    y = a._y * b._x + a._x * b._y;
                return new AuxAngle(y, x);
            }
            return a;
        }

        /// <summary>
        /// Convert degrees to an <see cref="AuxAngle"/>.
        /// </summary>
        /// <param name="d">The angle measured in degrees.</param>
        /// <returns>The corresponding <see cref="AuxAngle"/>.</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static AuxAngle FromDegrees(double d)
        {
            SinCosd(d, out var y, out var x);
            return new AuxAngle(y, x);
        }

        /// <summary>
        /// Convert radians to an <see cref="AuxAngle"/>.
        /// </summary>
        /// <param name="r">The angle measured in radians.</param>
        /// <returns>The corresponding <see cref="AuxAngle"/>.</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static AuxAngle FromRadians(double r)
            => new AuxAngle(Sin(r), Cos(r));

        /// <summary>
        /// Convert lambertian to an <see cref="AuxAngle"/>.
        /// </summary>
        /// <param name="psi">The lambertian of the angle.</param>
        /// <returns>The corresponding <see cref="AuxAngle"/>.</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static AuxAngle FromLam(double psi)
            => new AuxAngle(Sinh(psi));

        /// <summary>
        /// Convert lambertian in degrees to an <see cref="AuxAngle"/>.
        /// </summary>
        /// <param name="psid">The lambertian of the angle in degrees.</param>
        /// <returns>The corresponding <see cref="AuxAngle"/>.</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static AuxAngle FromLamd(double psid)
            => new AuxAngle(Sinh(psid * Degree));

        /// <summary>
        /// Gets the string representation of current <see cref="AuxAngle"/> instance.
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return $"X = {X}, Y = {Y}, {Degrees:F2} degs, {Radians:F2} rads";
        }
    }
}