using System;

namespace GeographicLib
{
    /// <summary>
    /// Represents an ordered pair of double precision floating-point x- and y-coordinates that defines a point in a two-dimensional plane.
    /// </summary>
    public readonly struct Point : IEquatable<Point>
    {
        /// <summary>
        /// A <see cref="Point"/> with <i>X</i> = 0 and <i>Y</i> = 0.
        /// </summary>
        public static readonly Point Zero = new Point(0, 0);

        /// <summary>
        /// Initializes a new instance of the <see cref="Point"/> struct with the specified <paramref name="x"/> and <paramref name="y"/> components. 
        /// </summary>
        /// <param name="x">The <i>X</i> component of the <see cref="Point"/>.</param>
        /// <param name="y">The <i>Y</i> component of the <see cref="Point"/>.</param>
        public Point(double x, double y)
        {
            X = x;
            Y = y;
        }

        /// <summary>
        /// The <i>X</i> component of the <see cref="Point"/>.
        /// </summary>
        public double X { get; }

        /// <summary>
        /// The <i>Y</i> component of the <see cref="Point"/>.
        /// </summary>
        public double Y { get; }

        /// <inheritdoc/>
        public bool Equals(Point other)
        {
            return X == other.X && Y == other.Y;
        }

        /// <inheritdoc/>
        public override bool Equals(object obj)
        {
            if (obj is Point p)
            {
                return Equals(p);
            }
            return false;
        }

        /// <inheritdoc/>
        public static bool operator ==(Point p1, Point p2)
        {
            return p1.Equals(p2);
        }

        /// <inheritdoc/>
        public static bool operator !=(Point p1, Point p2)
        {
            return !(p1 == p2);
        }

        /// <inheritdoc/>
        public override int GetHashCode()
        {
            long x = ((Bit64)X).Int64,
                 y = ((Bit64)Y).Int64;

            unchecked
            {
                return ((int)(x & 0xffffffff)) ^ ((int)(x >> 32)) ^
                       ((int)(y & 0xffffffff)) ^ ((int)(y >> 32));
            }
        }

        /// <summary>
        /// Gets the string representation of current <see cref="Point"/> instance.
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return $"X={X}, Y={Y}";
        }
    }
}