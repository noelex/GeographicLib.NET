using System;
using System.Collections.Generic;
using System.Text;

namespace GeographicLib.SphericalHarmonics
{
    /// <summary>
    /// Spherical harmonic series with two corrections to the coefficients.
    /// </summary>
    /// <remarks>
    /// This classes is similar to <see cref="SphericalHarmonic"/>, 
    /// except that the coefficients <i>Cnm</i> are replaced by <i>Cnm</i> + <i>tau</i> ' <i>C'nm</i> + <i>tau</i>'' <i>C''nm</i> 
    /// (and similarly for <i>Snm</i>).
    /// </remarks>
    public class SphericalHarmonic2
    {
        private readonly Memory<SphericalEngine.Coeff> _c;
        private readonly double _a;
        private readonly Normalization _norm;

        /// <summary>
        /// Constructor with a subset of coefficients specified.
        /// </summary>
        /// <param name="C">the coefficients <i>Cnm</i>.</param>
        /// <param name="S">the coefficients <i>Snm</i>.</param>
        /// <param name="N">the maximum degree and order of the sum.</param>
        /// <param name="C1">the coefficients <i>C'nm</i>.</param>
        /// <param name="S1">the coefficients <i>S'nm</i>.</param>
        /// <param name="N1">the maximum degree and order of the first correction coefficients <i>C'nm</i> and <i>S'nm</i>.</param>
        /// <param name="C2">the coefficients <i>C''nm</i>.</param>
        /// <param name="S2">the coefficients <i>S''nm</i>.</param>
        /// <param name="N2">the maximum degree and order of the first correction coefficients <i>C''nm</i> and <i>S''nm</i>.</param>
        /// <param name="a">the reference radius appearing in the definition of the sum.</param>
        /// <param name="norm">the normalization for the associated Legendre polynomials,
        /// either <see cref="Normalization.Full"/> (the default) or <see cref="Normalization.Schmidt"/>.</param>
        /// <remarks>
        /// See <see cref="SphericalHarmonic"/> for the way the coefficients should be stored.
        /// <para>
        /// The class stores pointers to the first elements of <i>C</i>, <i>S</i>, <i>C</i> ', <i>S</i> ', <i>C</i> '' and <i>S</i> ''.
        /// These arrays should not be altered or destroyed during the lifetime of a <see cref="SphericalHarmonic"/> object.
        /// </para>
        /// </remarks>
        public SphericalHarmonic2(ReadOnlyMemory<double> C,
               ReadOnlyMemory<double> S,
               int N,
               ReadOnlyMemory<double> C1,
               ReadOnlyMemory<double> S1,
               int N1,
               ReadOnlyMemory<double> C2,
               ReadOnlyMemory<double> S2,
               int N2,
               double a, Normalization norm = Normalization.Full)
            :this(new SphericalEngine.Coeff(C, S, N),
                  new SphericalEngine.Coeff(C1, S1, N1),
                  new SphericalEngine.Coeff(C2, S2, N2), a, norm)
        { }

        /// <summary>
        /// Constructor with a subset of coefficients specified.
        /// </summary>
        /// <param name="C">the coefficients <i>Cnm</i>.</param>
        /// <param name="S">the coefficients <i>Snm</i>.</param>
        /// <param name="N">the degree used to determine the layout of <i>C</i> and <i>S</i>.</param>
        /// <param name="nmx">the maximum degree used in the sum. The sum over <i>n</i> is from <c>0</c> thru <i>nmx</i>.</param>
        /// <param name="mmx">the maximum order used in the sum. The sum over <i>m</i> is from <c>0</c> thru min(<i>n</i>, <i>mmx</i>).</param>
        /// <param name="C1">the coefficients <i>C'nm</i>.</param>
        /// <param name="S1">the coefficients <i>S'nm</i>.</param>
        /// <param name="N1">the degree used to determine the layout of <i>C</i> ' and <i>S</i> '.</param>
        /// <param name="nmx1">the maximum degree used for <i>C</i> ' and <i>S</i> '.</param>
        /// <param name="mmx1">the maximum order used for <i>C</i> ' and <i>S</i> '.</param>
        /// <param name="C2">the coefficients <i>C''nm</i>.</param>
        /// <param name="S2">the coefficients <i>S''nm</i>.</param>
        /// <param name="N2">the degree used to determine the layout of <i>C</i> '' and <i>S</i> ''.</param>
        /// <param name="nmx2">the maximum degree used for <i>C</i> '' and <i>S</i> ''.</param>
        /// <param name="mmx2">the maximum order used for <i>C</i> '' and <i>S</i> ''.</param>
        /// <param name="a">the reference radius appearing in the definition of the sum.</param>
        /// <param name="norm">the normalization for the associated Legendre polynomials,
        /// either <see cref="Normalization.Full"/> (the default) or <see cref="Normalization.Schmidt"/>.</param>
        /// <remarks>
        /// See <see cref="SphericalHarmonic"/> for the way the coefficients should be stored.
        /// <para>
        /// The class stores pointers to the first elements of <i>C</i>, <i>S</i>, <i>C</i> ', <i>S</i> ', <i>C</i> '' and <i>S</i> ''.
        /// These arrays should not be altered or destroyed during the lifetime of a <see cref="SphericalHarmonic"/> object.
        /// </para>
        /// </remarks>
        public SphericalHarmonic2(ReadOnlyMemory<double> C,
               ReadOnlyMemory<double> S,
               int N, int nmx, int mmx,
               ReadOnlyMemory<double> C1,
               ReadOnlyMemory<double> S1,
               int N1, int nmx1, int mmx1,
               ReadOnlyMemory<double> C2,
               ReadOnlyMemory<double> S2,
               int N2, int nmx2, int mmx2,
               double a, Normalization norm = Normalization.Full)
            : this(new SphericalEngine.Coeff(C, S, N, nmx, mmx),
                   new SphericalEngine.Coeff(C1, S1, N1, nmx1, mmx1),
                   new SphericalEngine.Coeff(C2, S2, N2, nmx2, mmx2), a, norm)
        { }

        /// <summary>
        /// Constructor with a subset of coefficients specified.
        /// </summary>
        /// <param name="coeff0">A <see cref="SphericalEngine.Coeff"/> object containing <i>Cnm</i> and <i>Snm</i> coefficients.</param>
        /// <param name="coeff1">A <see cref="SphericalEngine.Coeff"/> object containing <i>C'nm</i> and <i>S'nm</i> coefficients.</param>
        /// <param name="coeff2">A <see cref="SphericalEngine.Coeff"/> object containing <i>C''nm</i> and <i>S''nm</i> coefficients.</param>
        /// <param name="a">the reference radius appearing in the definition of the sum.</param>
        /// <param name="norm">the normalization for the associated Legendre polynomials,
        /// either <see cref="Normalization.Full"/> (the default) or <see cref="Normalization.Schmidt"/>.</param>
        /// <remarks>See <see cref="SphericalHarmonic"/> for the way the coefficients should be stored.</remarks>
        public SphericalHarmonic2(SphericalEngine.Coeff coeff0, 
                SphericalEngine.Coeff coeff1, SphericalEngine.Coeff coeff2, double a, Normalization norm = Normalization.Full)
        {
            _a = a;
            _norm = norm;

            if (!(coeff1.N <= coeff0.N && coeff2.N <= coeff0.N))
                throw new GeographicException("N1 and N2 cannot be larger that N");
            if (!(coeff1.Nmx <= coeff0.Nmx && coeff2.Nmx <= coeff0.Nmx))
                throw new GeographicException("nmx1 and nmx2 cannot be larger that nmx");
            if (!(coeff1.Mmx <= coeff0.Mmx && coeff2.Mmx <= coeff0.Mmx))
                throw new GeographicException("mmx1 and mmx2 cannot be larger that mmx");

            _c = new[] { coeff0, coeff1, coeff2 };
        }

        /// <summary>
        /// Compute a spherical harmonic sum with a correction term.
        /// </summary>
        /// <param name="tau1">multiplier for correction coefficients <i>C</i> ' and <i>S</i> '.</param>
        /// <param name="tau2">multiplier for correction coefficients <i>C</i> '' and <i>S</i> ''.</param>
        /// <param name="x"><i>x</i> component of the cartesian coordinate.</param>
        /// <param name="y"><i>y</i> component of the cartesian coordinate.</param>
        /// <param name="z"><i>z</i> component of the cartesian coordinate.</param>
        /// <returns><i>V</i>, the spherical harmonic sum.</returns>
        /// <remarks>
        /// This routine requires constant memory and thus never throws an exception.
        /// </remarks>
        public double Evaluate(double tau1, double tau2, double x, double y, double z)
        {
            Span<double> f = stackalloc double[] { 1, tau1, tau2 };

            return SphericalEngine.Value(false, _norm, _c.Span, f, x, y, z, _a, out _, out _, out _);
        }

        /// <summary>
        /// Compute a spherical harmonic sum with a correction term and its gradient.
        /// </summary>
        /// <param name="tau1">multiplier for correction coefficients <i>C</i> ' and <i>S</i> '.</param>
        /// <param name="tau2">multiplier for correction coefficients <i>C</i> '' and <i>S</i> ''.</param>
        /// <param name="x"><i>x</i> component of the cartesian coordinate.</param>
        /// <param name="y"><i>y</i> component of the cartesian coordinate.</param>
        /// <param name="z"><i>z</i> component of the cartesian coordinate.</param>
        /// <param name="gradx"><i>x</i> component of the gradient.</param>
        /// <param name="grady"><i>y</i> component of the gradient.</param>
        /// <param name="gradz"><i>z</i> component of the gradient.</param>
        /// <returns><i>V</i>, the spherical harmonic sum.</returns>
        /// <remarks>
        /// This is the same as <see cref="Evaluate(double, double, double, double, double)"/>,
        /// except that the components of the gradients of the sum in the <i>x</i>, <i>y</i>, and <i>z</i> directions are computed.
        /// This routine requires constant memory and thus never throws an exception.
        /// </remarks>
        public double Evaluate(double tau1, double tau2, double x, double y, double z,
                               out double gradx, out double grady, out double gradz)
        {
            Span<double> f = stackalloc double[] { 1, tau1, tau2 };

            return SphericalEngine.Value(true, _norm, _c.Span, f, x, y, z, _a, out gradx, out grady, out gradz);
        }

        /// <summary>
        /// Create a <see cref="CircularEngine"/> to allow the efficient evaluation of several points on a circle of latitudeat a fixed value of <i>tau</i>.
        /// </summary>
        /// <param name="tau1">multiplier for correction coefficients <i>C</i> ' and <i>S</i> '.</param>
        /// <param name="tau2">multiplier for correction coefficients <i>C</i> '' and <i>S</i> ''.</param>
        /// <param name="p">the radius of the circle.</param>
        /// <param name="z">the height of the circle above the equatorial plane.</param>
        /// <param name="gradp">if <see langword="true"/> the returned object will be able to compute the gradient of the sum.</param>
        /// <returns>A <see cref="CircularEngine"/> instance.</returns>
        /// <remarks>
        /// <see cref="Evaluate(double, double, double, double, double)"/> exchanges the order of the sums in the definition,
        /// i.e., ∑<i>n</i> = 0..<i>N</i> ∑<i>m</i> = 0..<i>n</i> becomes ∑<i>m</i> = 0..<i>N</i> ∑<i>n</i> = <i>m</i>..<i>N</i>.
        /// <see cref="Circle(double, double, double, double, bool)"/> performs the inner sum over degree n (which entails about <i>N</i>^2 operations).
        /// Calling <see cref="CircularEngine.Evaluate(double)"/> on the returned object performs the outer sum over the order <i>m</i> (about <i>N</i> operations).
        /// </remarks>
        public CircularEngine Circle(double tau1, double tau2, double p, double z, bool gradp)
        {
            Span<double> f = stackalloc double[] { 1, tau1, tau2 };

            return SphericalEngine.Circle(true, _norm, _c.Span, f, p, z, _a);
        }

        /// <summary>
        /// Gets the zeroth <see cref="SphericalEngine.Coeff"/> object.
        /// </summary>
        public ref SphericalEngine.Coeff Coefficients => ref _c.Span[0];

        /// <summary>
        /// Gets the first <see cref="SphericalEngine.Coeff"/> object.
        /// </summary>
        public ref SphericalEngine.Coeff Coefficients1 => ref _c.Span[1];

        /// <summary>
        /// Gets the second <see cref="SphericalEngine.Coeff"/> object.
        /// </summary>
        public ref SphericalEngine.Coeff Coefficients2 => ref _c.Span[2];
    }
}
