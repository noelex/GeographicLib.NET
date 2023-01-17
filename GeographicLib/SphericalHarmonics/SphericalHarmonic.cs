using System;
using System.Collections.Generic;
using System.Text;

namespace GeographicLib.SphericalHarmonics
{
    /// <summary>
    /// Spherical harmonic series.
    /// </summary>
    /// <remarks>
    /// This class evaluates the spherical harmonic sum
    /// <code>
    /// V(x, y, z) = sum(n = 0..N)[ q^(n+1) * sum(m = 0..n)[
    ///   (C[n, m] * cos(m* lambda) + S[n, m] * sin(m* lambda)) *
    ///   P[n, m](cos(theta)) ] ]
    /// </code>
    /// <para>
    /// where
    /// <list type="bullet">
    /// <item><i>p</i>^2 = <i>x</i>^2 + <i>y</i>^2,</item>
    /// <item><i>r</i>^2 = <i>p</i>^2 + <i>z</i>^2,</item>
    /// <item><i>q</i> = <i>a</i> / <i>r</i>,</item>
    /// <item>θ = atan2(<i>p</i>, <i>z</i>) = the spherical <i>colatitude</i>,</item>
    /// <item>λ = atan2(<i>y</i>, <i>x</i>) = the longitude.</item>
    /// <item>P<i>nm</i>(<i>t</i>) is the associated Legendre polynomial of degree <i>n</i> and order <i>m</i>.</item>
    /// </list>
    /// </para>
    /// <para>
    /// Two normalizations are supported for P<i>nm</i>
    /// <list type="bullet">
    /// <item>fully normalized denoted by <see cref="Normalization.Full"/>.</item>
    /// <item>Schmidt semi-normalized denoted by <see cref="Normalization.Schmidt"/>.</item>
    /// </list>
    /// </para>
    /// <para>
    /// Clenshaw summation is used for the sums over both <i>n</i> and <i>m</i>.
    /// This allows the computation to be carried out without the need for any temporary arrays.
    /// See <see cref="SphericalEngine"/> for more information on the implementation.
    /// </para>
    /// <para>
    /// References:
    /// <list type="bullet">
    /// <item>
    /// C. W. Clenshaw, <a href="https://doi.org/10.1090/S0025-5718-1955-0071856-0">
    /// A note on the summation of Chebyshev series</a>, Math. Tables Aids Comput. 9(51), 118–120 (1955).
    /// </item>
    /// <item>
    /// R. E. Deakin, Derivatives of the earth's potentials, Geomatics Research Australasia 68, 31–60, (June 1998).
    /// </item>
    /// <item>
    /// <a href="https://archive.org/details/HeiskanenMoritz1967PhysicalGeodesy">Physical Geodesy (Freeman, San Francisco, 1967)</a>. (See Sec. 1-14, for a definition of Pbar.)
    /// </item>
    /// <item>
    /// S. A. Holmes and W. E. Featherstone, 
    /// <a href="https://doi.org/10.1007/s00190-002-0216-2">
    /// A unified approach to the Clenshaw summation and the recursive computation of very high degree and order normalised associated Legendre functions</a>,
    /// J. Geodesy 76(5), 279–299 (2002).
    /// </item>
    /// <item>
    /// C. C. Tscherning and K. Poder, Some geodetic applications of Clenshaw summation, Boll. Geod. Sci. Aff. 41(4), 349–375 (1982).
    /// </item>
    /// </list>
    /// </para>
    /// </remarks>
    public class SphericalHarmonic
    {
        private readonly Memory<SphericalEngine.Coeff> _c;
        private readonly double _a;
        private readonly Normalization _norm;

        /// <summary>
        /// Constructor with a full set of coefficients specified.
        /// </summary>
        /// <param name="C">the coefficients <i>Cnm</i>.</param>
        /// <param name="S">the coefficients <i>Snm</i>.</param>
        /// <param name="N">the maximum degree and order of the sum</param>
        /// <param name="a">the reference radius appearing in the definition of the sum.</param>
        /// <param name="norm">the normalization for the associated Legendre polynomials,
        /// either <see cref="Normalization.Full"/> (the default) or <see cref="Normalization.Schmidt"/>.</param>
        /// <remarks>
        /// The coefficients <i>Cnm</i> and <i>Snm</i> are stored in the one-dimensional vectors
        /// <i>C</i> and <i>S</i> which must contain (<i>N</i> + 1)(<i>N</i> + 2)/2 and <i>N</i> (<i>N</i> + 1)/2 elements,
        /// respectively, stored in "column-major" order. 
        /// Thus for <i>N</i> = 3, the order would be: <i>C</i>00, <i>C</i>10, <i>C</i>20, <i>C</i>30, <i>C</i>11, <i>C</i>21, <i>C</i>31, <i>C</i>22, <i>C</i>32, <i>C</i>33.
        /// In general the (<i>n</i>,<i>m</i>) element is at index <i>m</i> <i>N</i> − <i>m</i> (<i>m</i> − 1)/2 + <i>n</i>.
        /// The layout of <i>S</i> is the same except that the first column is omitted (since the <i>m</i> = 0 terms never contribute to the sum)
        /// and the 0th element is <i>S</i>11.
        /// <para>
        /// The class stores <i>pointers</i> to the first elements of <i>C</i> and <i>S</i>.
        /// These arrays should not be altered or destroyed during the lifetime of a <see cref="SphericalHarmonic"/> object.
        /// </para>
        /// </remarks>
        public SphericalHarmonic(ReadOnlyMemory<double> C, ReadOnlyMemory<double> S,
                                 int N, double a, Normalization norm = Normalization.Full)
            :this(new SphericalEngine.Coeff(C, S, N), a, norm) { }

        /// <summary>
        /// Constructor with a subset of coefficients specified.
        /// </summary>
        /// <param name="C">the coefficients <i>Cnm</i>.</param>
        /// <param name="S">the coefficients <i>Snm</i>.</param>
        /// <param name="N">the maximum degree and order of the sum</param>
        /// <param name="nmx">the maximum degree used in the sum. The sum over <i>n</i> is from <c>0</c> thru <i>nmx</i>.</param>
        /// <param name="mmx">the maximum order used in the sum. The sum over <i>m</i> is from <c>0</c> thru min(<i>n</i>, <i>mmx</i>).</param>
        /// <param name="a">the reference radius appearing in the definition of the sum.</param>
        /// <param name="norm">the normalization for the associated Legendre polynomials,
        /// either <see cref="Normalization.Full"/> (the default) or <see cref="Normalization.Schmidt"/>.</param>
        /// <remarks>
        /// <para>
        /// The class stores <i>pointers</i> to the first elements of <i>C</i> and <i>S</i>.
        /// These arrays should not be altered or destroyed during the lifetime of a <see cref="SphericalHarmonic"/> object.
        /// </para>
        /// </remarks>
        public SphericalHarmonic(ReadOnlyMemory<double> C, ReadOnlyMemory<double> S,
                                 int N, int nmx, int mmx, double a, Normalization norm = Normalization.Full)
            : this(new SphericalEngine.Coeff(C, S, N, nmx, mmx), a, norm) { }

        /// <summary>
        /// Constructor with a subset of coefficients specified.
        /// </summary>
        /// <param name="coeff">A <see cref="SphericalEngine.Coeff"/> object containing <i>Cnm</i> and <i>Snm</i> coefficients.</param>
        /// <param name="a">the reference radius appearing in the definition of the sum.</param>
        /// <param name="norm">the normalization for the associated Legendre polynomials,
        /// either <see cref="Normalization.Full"/> (the default) or <see cref="Normalization.Schmidt"/>.</param>
        public SphericalHarmonic(SphericalEngine.Coeff coeff, double a, Normalization norm = Normalization.Full)
        {
            _a = a;
            _norm = norm;
            _c = new[] { coeff };
        }

        /// <summary>
        /// Compute the spherical harmonic sum.
        /// </summary>
        /// <param name="x"><i>x</i> component of the cartesian coordinate.</param>
        /// <param name="y"><i>y</i> component of the cartesian coordinate.</param>
        /// <param name="z"><i>z</i> component of the cartesian coordinate.</param>
        /// <returns><i>V</i>, the spherical harmonic sum.</returns>
        /// <remarks>
        /// This routine requires constant memory and thus never throws an exception.
        /// </remarks>
        public double Evaluate(double x, double y, double z)
        {
            Span<double> f = stackalloc double[] { 1 };

            return SphericalEngine.Value(false, _norm, _c.Span, f, x, y, z, _a, out _, out _, out _);
        }

        /// <summary>
        /// Compute the spherical harmonic sum and its gradient.
        /// </summary>
        /// <param name="x"><i>x</i> component of the cartesian coordinate.</param>
        /// <param name="y"><i>y</i> component of the cartesian coordinate.</param>
        /// <param name="z"><i>z</i> component of the cartesian coordinate.</param>
        /// <param name="gradx"><i>x</i> component of the gradient.</param>
        /// <param name="grady"><i>y</i> component of the gradient.</param>
        /// <param name="gradz"><i>z</i> component of the gradient.</param>
        /// <returns><i>V</i>, the spherical harmonic sum.</returns>
        /// <remarks>
        /// This is the same as <see cref="Evaluate(double, double, double)"/>,
        /// except that the components of the gradients of the sum in the <i>x</i>, <i>y</i>, and <i>z</i> directions are computed.
        /// This routine requires constant memory and thus never throws an exception.
        /// </remarks>
        public double Evaluate(double x, double y, double z,
                               out double gradx, out double grady, out double gradz)
        {
            Span<double> f = stackalloc double[] { 1 };

            return SphericalEngine.Value(true, _norm, _c.Span, f, x, y, z, _a, out gradx, out grady, out gradz);
        }

        /// <summary>
        /// Create a <see cref="CircularEngine"/> to allow the efficient evaluation of several points on a circle of latitude.
        /// </summary>
        /// <param name="p">the radius of the circle.</param>
        /// <param name="z">the height of the circle above the equatorial plane.</param>
        /// <param name="gradp">if <see langword="true"/> the returned object will be able to compute the gradient of the sum.</param>
        /// <returns>A <see cref="CircularEngine"/> instance.</returns>
        /// <remarks>
        /// <see cref="Evaluate(double, double, double)"/> exchanges the order of the sums in the definition,
        /// i.e., ∑<i>n</i> = 0..<i>N</i> ∑<i>m</i> = 0..<i>n</i> becomes ∑<i>m</i> = 0..<i>N</i> ∑<i>n</i> = <i>m</i>..<i>N</i>.
        /// <see cref="Circle(double, double, bool)"/> performs the inner sum over degree n (which entails about <i>N</i>^2 operations).
        /// Calling <see cref="CircularEngine.Evaluate(double)"/> on the returned object performs the outer sum over the order <i>m</i> (about <i>N</i> operations).
        /// </remarks>
        public CircularEngine Circle(double p, double z, bool gradp)
        {
            Span<double> f = stackalloc double[] { 1 };

            return SphericalEngine.Circle(true, _norm, _c.Span, f, p, z, _a);
        }

        /// <summary>
        /// Gets the zeroth <see cref="SphericalEngine.Coeff"/> object.
        /// </summary>
        public ref SphericalEngine.Coeff Coefficients => ref _c.Span[0];
    }
}
