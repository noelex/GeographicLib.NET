namespace GeographicLib
{
    /// <summary>
    /// Provides elliptic integrals and functions.
    /// </summary>
    /// <remarks>
    /// This class provides the elliptic functions and integrals needed for <see cref="Ellipsoid"/>,
    /// <see cref="GeodesicExact"/>, and <see cref="Projections.TransverseMercatorExact"/>. Two categories of function
    /// are provided:
    /// <list type="bullet">
    /// <item>static functions to compute symmetric elliptic integrals(https://dlmf.nist.gov/19.16.i)</item>
    /// <item>member functions to compute Legrendre's elliptic integrals(https://dlmf.nist.gov/19.2.ii) and the
    /// Jacobi elliptic functions (https://dlmf.nist.gov/22.2).</item>
    /// </list>
    /// <para>
    /// In the latter case, an object is constructed giving the modulus <i>k</i> (and
    ///  optionally the parameter α^2).  The modulus is always
    /// passed as its square<i>k</i>^2 which allows <i>k</i> to be pure
    /// imaginary (<i>k</i>^2 &lt; 0).  (Confusingly, Abramowitz and
    /// Stegun call <i>m</i> = <i>k</i>^2 the "parameter" and <i>n</i>  =
    /// α^2 the "characteristic".)
    /// </para>
    /// <para>
    /// In geodesic applications, it is convenient to separate the incomplete
    /// integrals into secular and periodic components, e.g.,
    /// E(φ, k) = (2 E(k) / π) [ φ + Δ E(φ, k)]
    /// where Δ<i>E</i>(φ, <i>k</i>) is an odd periodic function with period
    /// π.
    /// The computation of the elliptic integrals uses the algorithms given in
    /// - B.C.Carlson,
    ///   <a href = "https://doi.org/10.1007/BF02198293" > Computation of real or
    /// complex elliptic integrals</a>, Numerical Algorithms 10, 13--26 (1995);
    /// <a href="https://arxiv.org/abs/math/9409227">preprint</a>.
    /// with the additional optimizations given in https://dlmf.nist.gov/19.36.i.
    /// The computation of the Jacobi elliptic functions uses the algorithm given
    /// in
    /// - R.Bulirsch,
    ///   <a href = "https://doi.org/10.1007/BF01397975" > Numerical Calculation of
    /// Elliptic Integrals and Elliptic Functions</a>, Numericshe Mathematik 7,
    ///   78--90 (1965).
    /// The notation follows <a href = "https://dlmf.nist.gov/19" />  and <a href = " https://dlmf.nist.gov/22" />.
    /// </para>
    /// </remarks>
    public partial class EllipticFunction
    {
        // Max depth required for sncndn; probably 5 is enough.
        private Priv _priv;

        /// <summary>
        /// Constructor specifying the modulus and parameter.
        /// </summary>
        /// <param name="k2">the square of the modulus <i>k</i>^2. <i>k</i>^2 must lie in (-∞, 1].</param>
        /// <param name="alpha2">the parameter α^2. α^2 must lie in (-∞, 1].</param>
        /// <remarks>
        /// If only elliptic integrals of the first and second kinds are needed,
        /// then set α^2 = 0 (the default value); in this case, we
        /// have Π(φ, 0, <i>k</i>) = <i>F</i>(φ, <i>k</i>), <i>G</i>(φ, 0, <i>k</i>) = <i>E</i>(φ, <i>k</i>),
        /// and <i>H</i>(φ, 0, <i>k</i>) = <i>F</i>(φ, <i>k</i>) - <i>D</i>(φ, <i>k</i>).
        /// </remarks>
        public EllipticFunction(double k2 = 0, double alpha2 = 0) => Reset(k2, alpha2);

        /// <summary>
        /// Constructor specifying the modulus and parameter and their complements.
        /// </summary>
        /// <param name="k2">the square of the modulus <i>k</i>^2. <i>k</i>^2 must lie in (-∞, 1].</param>
        /// <param name="alpha2">the parameter α^2. α^2 must lie in (-∞, 1].</param>
        /// <param name="kp2">the complementary modulus squared <i>k'</i>^2 = 1 - <i>k</i>^2.  This must lie in [0, ∞).</param>
        /// <param name="alphap2">the complementary parameter α'^2 = 1 - α^2.  This must lie in [0, ∞).</param>
        /// <remarks>
        /// The arguments must satisfy <paramref name="k2"/> + <paramref name="kp2"/> = 1 and <paramref name="alpha2"/> + <paramref name="alphap2"/>
        /// = 1.  (No checking is done that these conditions are met.)  This
        /// constructor is provided to enable accuracy to be maintained, e.g., when
        /// <i>k</i> is very close to unity.
        /// </remarks>
        public EllipticFunction(double k2, double alpha2, double kp2, double alphap2) => Reset(k2, alpha2, kp2, alphap2);

        /// <summary>
        /// Gets a value representing the square of the modulus <i>k</i>^2.
        /// </summary>
        public double K2 => _priv.K2;

        /// <summary>
        /// Gets a value representing the complementary modulus squared <i>k'</i>^2 = 1 - <i>k</i>^2
        /// </summary>
        public double Kp2 => _priv.Kp2;

        /// <summary>
        /// Gets a value representing the parameter α^2.
        /// </summary>
        public double Alpha2 => _priv.Alpha2;

        /// <summary>
        /// Gets a value representing the complementary parameter α'^2 = 1 - α^2
        /// </summary>
        public double Alphap2 => _priv.Alphap2;

        #region Complete elliptic integrals

        /// <summary>
        /// The complete integral of the first kind.
        /// </summary>
        /// <returns><i>K</i>(<i>k</i>).</returns>
        public double K() => _priv.K();

        /// <summary>
        /// The complete integral of the second kind.
        /// </summary>
        /// <returns><i>E</i>(<i>k</i>).</returns>
        public double E() => _priv.E();

        /// <summary>
        /// Jahnke's complete integral.
        /// </summary>
        /// <returns><i>D</i>(<i>k</i>).</returns>
        public double D() => _priv.D();

        /// <summary>
        /// The difference between the complete integrals of the first and second kinds.
        /// </summary>
        /// <returns><i>K</i>(<i>k</i>) - <i>E</i>(<i>k</i>).</returns>
        public double KE() => _priv.KE();

        /// <summary>
        /// The complete integral of the third kind.
        /// </summary>
        /// <returns>Π(α^2, <i>k</i>).</returns>
        public double Pi() => _priv.Pi();

        /// <summary>
        /// Legendre's complete geodesic longitude integral.
        /// </summary>
        /// <returns><i>G</i>(α^2, <i>k</i>).</returns>
        public double G() => _priv.G();

        /// <summary>
        /// Cayley's complete geodesic longitude difference integral.
        /// </summary>
        /// <returns><i>H</i>(α^2, <i>k</i>).</returns>
        public double H() => _priv.H();

        #endregion

        #region Incomplete elliptic integrals

        /// <summary>
        /// The incomplete integral of the first kind.
        /// </summary>
        /// <param name="phi">φ</param>
        /// <returns><i>F</i>(φ, <i>k</i>).</returns>
        /// <remarks>
        /// <i>F</i>(φ, <i>k</i>) is defined in <a href="https://dlmf.nist.gov/19.2.E4"/>.
        /// </remarks>
        public double F(double phi) => _priv.F(phi);

        /// <summary>
        /// The incomplete integral of the second kind.
        /// </summary>
        /// <param name="phi">φ</param>
        /// <returns><i>E</i>(φ, <i>k</i>).</returns>
        /// <remarks>
        /// <i>E</i>(φ, <i>k</i>) is defined in <a href="https://dlmf.nist.gov/19.2.E5"/>.
        /// </remarks>
        public double E(double phi) => _priv.E(phi);

        /// <summary>
        /// The incomplete integral of the second kind with the argument given in degrees.
        /// </summary>
        /// <param name="ang">in <i>degrees</i>.</param>
        /// <returns><i>E</i>(π <i>ang</i>/180, <i>k</i>).</returns>
        public double Ed(double ang) => _priv.Ed(ang);

        /// <summary>
        /// The inverse of the incomplete integral of the second kind.
        /// </summary>
        /// <param name="x"></param>
        /// <returns>φ = <i>E</i>^-1(<i>x</i>, <i>k</i>); i.e., the solution of such that <i>E</i>(φ, <i>k</i>) = <i>x</i>.</returns>
        public double Einv(double x) => _priv.Einv(x);

        /// <summary>
        /// The incomplete integral of the third kind.
        /// </summary>
        /// <param name="phi">φ</param>
        /// <returns>Π(φ, α^2, <i>k</i>).</returns>
        /// <remarks>
        /// Π(φ, α^2, <i>k</i>) is defined in <a href="https://dlmf.nist.gov/19.2.E7"/>.
        /// </remarks>
        public double Pi(double phi) => _priv.Pi(phi);

        /// <summary>
        /// Jahnke's incomplete elliptic integral.
        /// </summary>
        /// <param name="phi">φ</param>
        /// <returns><i>D</i>(φ, <i>k</i>).</returns>
        /// <remarks>
        /// <i>D</i>(φ, <i>k</i>) is defined in <a href="https://dlmf.nist.gov/19.2.E4"/>.
        /// </remarks>
        public double D(double phi) => _priv.D(phi);

        /// <summary>
        /// Legendre's geodesic longitude integral.
        /// </summary>
        /// <param name="phi">φ</param>
        /// <returns><i>G</i>(φ, α^2, <i>k</i>).</returns>
        /// <remarks>
        /// Legendre expresses the longitude of a point on the geodesic in terms of
        /// this combination of elliptic integrals in
        /// <a href="https://books.google.com/books?id=riIOAAAAQAAJ&amp;pg=PA181">
        /// Exercices de Calcul Intégral, Vol. 1 (1811), p. 181</a>.
        /// See geodellip for the expression for the longitude in terms of this
        /// function.
        /// </remarks>
        public double G(double phi) => _priv.G(phi);

        /// <summary>
        /// Cayley's geodesic longitude difference integral.
        /// </summary>
        /// <param name="phi">φ</param>
        /// <returns><i>H</i>(φ, α^2, <i>k</i>).</returns>
        /// <remarks>
        /// Cayley expresses the longitude difference of a point on the geodesic in
        /// terms of this combination of elliptic integrals in <a href="https://books.google.com/books?id=Zk0wAAAAIAAJ&amp;pg=PA333">
        /// Phil.Mag. <b>40</b> (1870), p. 333.</a> 
        /// See geodellip for the expression for the longitude in terms of this
        /// function.
        /// </remarks>
        public double H(double phi) => _priv.H(phi);

        /// <summary>
        /// The incomplete integral of the first kind in terms of Jacobi elliptic functions.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns><i>F</i>(φ, <i>k</i>) as though φ ∈ (-π, π].</returns>
        public double F(double sn, double cn, double dn) => _priv.F(sn, cn, dn);

        /// <summary>
        /// The incomplete integral of the second kind in terms of Jacobi elliptic functions.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns><i>E</i>(φ, <i>k</i>) as though φ ∈ (-π, π].</returns>
        public double E(double sn, double cn, double dn) => _priv.E(sn, cn, dn);

        /// <summary>
        /// The incomplete integral of the third kind in terms of Jacobi elliptic functions.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns>Π(φ, α^2, <i>k</i>) as though φ ∈ (-π, π].</returns>
        public double Pi(double sn, double cn, double dn) => _priv.Pi(sn, cn, dn);

        /// <summary>
        /// Jahnke's incomplete elliptic integral in terms of Jacobi elliptic functions.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns><i>D</i>(φ, <i>k</i>) as though φ ∈ (-π, π].</returns>
        public double D(double sn, double cn, double dn) => _priv.D(sn, cn, dn);

        /// <summary>
        /// Legendre's geodesic longitude integral in terms of Jacobi elliptic functions.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns><i>G</i>(φ, α^2, <i>k</i>) as though φ ∈ (-π, π].</returns>
        public double G(double sn, double cn, double dn) => _priv.G(sn, cn, dn);

        /// <summary>
        /// Cayley's geodesic longitude difference integral in terms of Jacobi elliptic functions.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns><i>H</i>(φ, α^2, <i>k</i>) as though φ ∈ (-π, π].</returns>
        public double H(double sn, double cn, double dn) => _priv.H(sn, cn, dn);

        /// <summary>
        /// The periodic incomplete integral of the first kind.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns>the periodic function π <i>F</i>(φ, <i>k</i>) / (2 <i>K</i>(<i>k</i>)) - φ.</returns>
        public double DeltaF(double sn, double cn, double dn) => _priv.DeltaF(sn, cn, dn);

        /// <summary>
        /// The periodic incomplete integral of the second kind.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns>the periodic function π <i>E</i>(φ, <i>k</i>) / (2 <i>E</i>(<i>k</i>)) - φ.</returns>
        public double DeltaE(double sn, double cn, double dn) => _priv.DeltaE(sn, cn, dn);

        /// <summary>
        /// The periodic inverse of the incomplete integral of the second kind.
        /// </summary>
        /// <param name="stau">sinτ</param>
        /// <param name="ctau">cosτ</param>
        /// <returns>the periodic function <i>E</i>^-1(τ (2 <i>E</i>(<i>k</i>)/π), <i>k</i>) - τ.</returns>
        public double DeltaEinv(double stau, double ctau) => _priv.DeltaEinv(stau, ctau);

        /// <summary>
        /// The periodic incomplete integral of the third kind.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns>the periodic function π Π(φ, α^2, <i>k</i>) / (2 Π(α^2, <i>k</i>)) - φ.</returns>
        public double DeltaPi(double sn, double cn, double dn) => _priv.DeltaPi(sn, cn, dn);

        /// <summary>
        /// The periodic Jahnke's incomplete elliptic integral.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns>the periodic function π <i>D</i>(φ, <i>k</i>) / (2 <i>D</i>(<i>k</i>)) - φ.</returns>
        public double DeltaD(double sn, double cn, double dn) => _priv.DeltaD(sn, cn, dn);

        /// <summary>
        /// Legendre's periodic geodesic longitude integral.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns>the periodic function π <i>G</i>(φ, <i>k</i>) / (2 <i>G</i>(<i>k</i>)) - φ.</returns>
        public double DeltaG(double sn, double cn, double dn) => _priv.DeltaG(sn, cn, dn);

        /// <summary>
        /// Cayley's periodic geodesic longitude difference integral.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns>the periodic function π <i>H</i>(φ, <i>k</i>) / (2 <i>H</i>(<i>k</i>)) -  φ.</returns>
        public double DeltaH(double sn, double cn, double dn) => _priv.DeltaH(sn, cn, dn);

        /// <summary>
        /// The Jacobi elliptic functions.
        /// </summary>
        /// <param name="x">the argument.</param>
        /// <param name="sn">sn(<i>x</i>, <i>k</i>)</param>
        /// <param name="cn">cn(<i>x</i>, <i>k</i>)</param>
        /// <param name="dn">dn(<i>x</i>, <i>k</i>)</param>
        /// <remarks>
        /// Implementation of methods given in
        /// <para>
        /// R. Bulirsch, Numerical Calculation of Elliptic Integrals and Elliptic Functions, Numericshe Mathematik 7, 78-90 (1965)
        /// </para>
        /// </remarks>
        public void Sncndn(double x, out double sn, out double cn, out double dn) => _priv.Sncndn(x, out sn, out cn, out dn);

        /// <summary>
        /// The Δ amplitude function.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <returns>Δ = sqrt(1 - <i>k</i>^2 sin^2φ)</returns>
        public double Delta(double sn, double cn) => _priv.Delta(sn, cn);

        #endregion

        #region Symmetric elliptic integrals

        /// <summary>
        /// Symmetric integral of the first kind <i>Rf</i>.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        /// <returns><i>Rf</i>(<i>x</i>, <i>y</i>, <i>z</i>).</returns>
        /// <remarks>
        /// <i>Rf</i> is defined in <a href="https://dlmf.nist.gov/19.16.E1"/>.
        /// <para>
        /// At most one of arguments, <paramref name="x"/>, <paramref name="y"/>, <paramref name="z"/>, can be zero and those
        /// arguments that are nonzero must be positive. 
        /// </para>
        /// <para>
        /// If one of the arguments is
        /// zero, it is more efficient to call the two-argument version of this
        /// function with the non-zero arguments.
        /// </para>
        /// </remarks>
        public static double RF(double x, double y, double z) => Priv.RF(x, y, z);

        /// <summary>
        /// Complete symmetric integral of the first kind, <i>Rf</i> with one argument zero.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns><i>Rf</i>(<i>x</i>, <i>y</i>, 0).</returns>
        /// <remarks>
        /// The arguments <paramref name="x"/> and <paramref name="y"/> must be positive.
        /// </remarks>
        public static double RF(double x, double y) => Priv.RF(x, y);

        /// <summary>
        /// Degenerate symmetric integral of the first kind <i>Rc</i>.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns><i>Rc</i>(x, y) = <i>Rf</i>(<i>x</i>, <i>y</i>, <i>y</i>).</returns>
        /// <remarks>
        /// <i>Rc</i> is defined in <a href="https://dlmf.nist.gov/19.16.E17"/>.
        /// <para>
        /// Requires <paramref name="x"/> >= 0 and <paramref name="y"/> &gt; 0.
        /// </para>
        /// </remarks>
        public static double RC(double x, double y) => Priv.RC(x, y);

        /// <summary>
        /// Symmetric integral of the second kind <i>Rg</i>.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        /// <returns><i>Rg</i>(<i>x</i>, <i>y</i>, <i>z</i>).</returns>
        /// <remarks>
        /// <i>Rg</i> is defined in Carlson, eq 1.5. See also <a href="https://dlmf.nist.gov/19.23.E6_5"/>.
        /// <para>
        /// At most one of arguments, <paramref name="x"/>, <paramref name="y"/>, <paramref name="z"/>, can be zero and those
        /// arguments that are nonzero must be positive.
        /// </para>
        /// If one of the arguments is zero, it is more efficient to call the
        /// two-argument version of this function with the non-zero arguments.
        /// </remarks>
        public static double RG(double x, double y, double z) => Priv.RG(x, y, z);

        /// <summary>
        /// Complete symmetric integral of the second kind, <i>Rg</i> with one argument zero.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns><i>Rg</i>(<i>x</i>, <i>y</i>, 0).</returns>
        /// <remarks>
        /// The arguments <paramref name="x"/> and <paramref name="y"/> must be positive.
        /// </remarks>
        public static double RG(double x, double y) => Priv.RG(x, y);

        /// <summary>
        /// Symmetric integral of the third kind <i>Rj</i>.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        /// <param name="p"></param>
        /// <returns><i>Rj</i>(<i>x</i>, <i>y</i>, <i>z</i>, <i>p</i>).</returns>
        /// <remarks>
        /// <i>Rd</i> is defined in <a href="https://dlmf.nist.gov/19.16.E2"/>.
        /// <para>
        /// Requires <paramref name="p"/> &gt; 0, and <paramref name="x"/>, <paramref name="y"/>, <paramref name="z"/>
        /// are nonnegative with at most one of them being 0.
        /// </para>
        /// </remarks>
        public static double RJ(double x, double y, double z, double p) => Priv.RJ(x, y, z, p);

        /// <summary>
        /// Degenerate symmetric integral of the third kind <i>Rd</i>.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        /// <returns><i>Rd</i>(<i>x</i>, <i>y</i>, <i>z</i>) = <i>Rj(<i>x</i>, <i>y</i>, <i>z</i>)</i>.</returns>
        /// <remarks>
        /// <i>Rd</i> is defined in <a href="https://dlmf.nist.gov/19.16.E5"/>.
        /// <para>
        /// Requires <paramref name="x"/>, <paramref name="y"/>, <paramref name="z"/> to be positive
        /// except that at most one of <paramref name="x"/> and <paramref name="y"/> can be 0.
        /// </para>
        /// </remarks>
        public static double RD(double x, double y, double z) => Priv.RD(x, y, z);

        #endregion

        /// <summary>
        /// Reset the modulus and parameter.
        /// </summary>
        /// <param name="k2">the new value of square of the modulus <i>k</i>^2 which must lie in (-∞, 1].</param>
        /// <param name="alpha2">the new value of parameter α^2. α^2 must lie in (-∞, 1].</param>
        public void Reset(double k2, double alpha2) => _priv.Reset(k2, alpha2, 1 - k2, 1 - alpha2);

        /// <summary>
        /// Reset the modulus and parameter supplying also their complements.
        /// </summary>
        /// <param name="k2">the square of the modulus <i>k</i>^2. <i>k</i>^2 must lie in (-∞, 1].</param>
        /// <param name="alpha2">the parameter α^2. α^2 must lie in (-∞, 1].</param>
        /// <param name="kp2">the complementary modulus squared <i>k'</i>^2 = 1 - <i>k</i>^2.  This must lie in [0, ∞).</param>
        /// <param name="alphap2">the complementary parameter α'^2 = 1 - α^2.  This must lie in [0, ∞).</param>
        /// <remarks>
        /// The arguments must satisfy <paramref name="k2"/> + <paramref name="kp2"/> = 1 and <paramref name="alpha2"/> + <paramref name="alphap2"/>
        /// = 1.  (No checking is done that these conditions are met.)  This
        /// constructor is provided to enable accuracy to be maintained, e.g., when
        /// is very small.
        /// </remarks>
        public void Reset(double k2, double alpha2, double kp2, double alphap2) => _priv.Reset(k2, alpha2, kp2, alphap2);
    }
}
