using System;
using System.Collections.Generic;
using System.Text;

using static System.Math;
using static GeographicLib.Macros;
using static GeographicLib.MathEx;

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
    public class EllipticFunction
    {
        // Max depth required for sncndn; probably 5 is enough.
        private const int num_ = 13;

        private double _k2, _kp2, _alpha2, _alphap2, _eps;
        private double _Kc, _Ec, _Dc, _Pic, _Gc, _Hc;

        private static readonly double
            tolRF = Pow(3 * DBL_EPSILON * 0.01, 1 / 8d),
            tolRG0 = 2.7 * Sqrt(DBL_EPSILON * 0.01),
            tolRD = Pow(0.2 * (DBL_EPSILON * 0.01), 1 / 8d),
            tolJAC = Sqrt(DBL_EPSILON * 0.01);

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
        public double K2 => _k2;

        /// <summary>
        /// Gets a value representing the complementary modulus squared <i>k'</i>^2 = 1 - <i>k</i>^2
        /// </summary>
        public double Kp2 => _kp2;

        /// <summary>
        /// Gets a value representing the parameter α^2.
        /// </summary>
        public double Alpha2 => _alpha2;

        /// <summary>
        /// Gets a value representing the complementary parameter α'^2 = 1 - α^2
        /// </summary>
        public double Alphap2 => _alphap2;

        #region Complete elliptic integrals

        /// <summary>
        /// The complete integral of the first kind.
        /// </summary>
        /// <returns><i>K</i>(<i>k</i>).</returns>
        public double K() => _Kc;

        /// <summary>
        /// The complete integral of the second kind.
        /// </summary>
        /// <returns><i>E</i>(<i>k</i>).</returns>
        public double E() => _Ec;

        /// <summary>
        /// Jahnke's complete integral.
        /// </summary>
        /// <returns><i>D</i>(<i>k</i>).</returns>
        public double D() => _Dc;

        /// <summary>
        /// The difference between the complete integrals of the first and second kinds.
        /// </summary>
        /// <returns><i>K</i>(<i>k</i>) - <i>E</i>(<i>k</i>).</returns>
        public double KE() => _k2 * _Dc;

        /// <summary>
        /// The complete integral of the third kind.
        /// </summary>
        /// <returns>Π(α^2, <i>k</i>).</returns>
        public double Pi() => _Pic;

        /// <summary>
        /// Legendre's complete geodesic longitude integral.
        /// </summary>
        /// <returns><i>G</i>(α^2, <i>k</i>).</returns>
        public double G() => _Gc;

        /// <summary>
        /// Cayley's complete geodesic longitude difference integral.
        /// </summary>
        /// <returns><i>H</i>(α^2, <i>k</i>).</returns>
        public double H() => _Hc;

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
        public double F(double phi)
        {
            double sn = Sin(phi), cn = Cos(phi), dn = Delta(sn, cn);

            return Abs(phi) < PI ? F(sn, cn, dn) :
              (DeltaF(sn, cn, dn) + phi) * K() / (PI / 2);
        }

        /// <summary>
        /// The incomplete integral of the second kind.
        /// </summary>
        /// <param name="phi">φ</param>
        /// <returns><i>E</i>(φ, <i>k</i>).</returns>
        /// <remarks>
        /// <i>E</i>(φ, <i>k</i>) is defined in <a href="https://dlmf.nist.gov/19.2.E5"/>.
        /// </remarks>
        public double E(double phi)
        {
            double sn = Sin(phi), cn = Cos(phi), dn = Delta(sn, cn);

            return Abs(phi) < PI ? E(sn, cn, dn) :
              (DeltaE(sn, cn, dn) + phi) * E() / (PI / 2);
        }

        /// <summary>
        /// The incomplete integral of the second kind with the argument given in degrees.
        /// </summary>
        /// <param name="ang">in <i>degrees</i>.</param>
        /// <returns><i>E</i>(π <i>ang</i>/180, <i>k</i>).</returns>
        public double Ed(double ang)
        {
            double n = Ceiling(ang / 360 - 0.5);
            ang -= 360 * n;
            SinCosd(ang, out var sn, out var cn);
            return E(sn, cn, Delta(sn, cn)) + 4 * E() * n;
        }

        /// <summary>
        /// The inverse of the incomplete integral of the second kind.
        /// </summary>
        /// <param name="x"></param>
        /// <returns>φ = <i>E</i>^-1(<i>x</i>, <i>k</i>); i.e., the solution of such that <i>E</i>(φ, <i>k</i>) = <i>x</i>.</returns>
        public double Einv(double x)
        {
            double n = Floor(x / (2 * _Ec) + 0.5);
            x -= 2 * _Ec * n;           // x now in [-ec, ec)
                                        // Linear approximation
            var phi = PI * x / (2 * _Ec); // phi in [-pi/2, pi/2)
                                          // First order correction
            phi -= _eps * Sin(2 * phi) / 2;
            // For kp2 close to zero use asin(x/_Ec) or
            // J. P. Boyd, Applied Math. and Computation 218, 7005-7013 (2012)
            // https://doi.org/10.1016/j.amc.2011.12.021
            for (int i = 0; i < num_ || GEOGRAPHICLIB_PANIC; ++i)
            {
                double
                  sn = Sin(phi),
                  cn = Cos(phi),
                  dn = Delta(sn, cn),
                  err = (E(sn, cn, dn) - x) / dn;

                phi -= err;
                if (!(Abs(err) > tolJAC))
                    break;
            }
            return n * PI + phi;
        }

        /// <summary>
        /// The incomplete integral of the third kind.
        /// </summary>
        /// <param name="phi">φ</param>
        /// <returns>Π(φ, α^2, <i>k</i>).</returns>
        /// <remarks>
        /// Π(φ, α^2, <i>k</i>) is defined in <a href="https://dlmf.nist.gov/19.2.E7"/>.
        /// </remarks>
        public double Pi(double phi)
        {
            double sn = Sin(phi), cn = Cos(phi), dn = Delta(sn, cn);
            return Abs(phi) < PI ? Pi(sn, cn, dn) :
              (DeltaPi(sn, cn, dn) + phi) * Pi() / (PI / 2);
        }

        /// <summary>
        /// Jahnke's incomplete elliptic integral.
        /// </summary>
        /// <param name="phi">φ</param>
        /// <returns><i>D</i>(φ, <i>k</i>).</returns>
        /// <remarks>
        /// <i>D</i>(φ, <i>k</i>) is defined in <a href="https://dlmf.nist.gov/19.2.E4"/>.
        /// </remarks>
        public double D(double phi)
        {
            double sn = Sin(phi), cn = Cos(phi), dn = Delta(sn, cn);
            return Abs(phi) < PI ? D(sn, cn, dn) :
              (DeltaD(sn, cn, dn) + phi) * D() / (PI / 2);
        }

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
        public double G(double phi)
        {
            double sn = Sin(phi), cn = Cos(phi), dn = Delta(sn, cn);
            return Abs(phi) < PI ? G(sn, cn, dn) :
              (DeltaG(sn, cn, dn) + phi) * G() / (PI / 2);
        }

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
        public double H(double phi)
        {
            double sn = Sin(phi), cn = Cos(phi), dn = Delta(sn, cn);
            return Abs(phi) < PI ? H(sn, cn, dn) :
              (DeltaH(sn, cn, dn) + phi) * H() / (PI / 2);
        }

        /// <summary>
        /// The incomplete integral of the first kind in terms of Jacobi elliptic functions.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns><i>F</i>(φ, <i>k</i>) as though φ ∈ (-π, π].</returns>
        public double F(double sn, double cn, double dn)
        {
            // Carlson, eq. 4.5 and
            // https://dlmf.nist.gov/19.25.E5
            double cn2 = cn * cn, dn2 = dn * dn,
              fi = cn2 != 0 ? Abs(sn) * RF(cn2, dn2, 1) : K();

            // Enforce usual trig-like symmetries
            if (SignBit(cn))
                fi = 2 * K() - fi;

            return CopySign(fi, sn);
        }

        /// <summary>
        /// The incomplete integral of the second kind in terms of Jacobi elliptic functions.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns><i>E</i>(φ, <i>k</i>) as though φ ∈ (-π, π].</returns>
        public double E(double sn, double cn, double dn)
        {
            double
              cn2 = cn * cn, dn2 = dn * dn, sn2 = sn * sn,
              ei = cn2 != 0 ?
              Abs(sn) * (_k2 <= 0 ?
                  // Carlson, eq. 4.6 and
                  // https://dlmf.nist.gov/19.25.E9
                  RF(cn2, dn2, 1) - _k2 * sn2 * RD(cn2, dn2, 1) / 3 :
                  (_kp2 >= 0 ?
                    // https://dlmf.nist.gov/19.25.E10
                    _kp2 * RF(cn2, dn2, 1) +
                    _k2 * _kp2 * sn2 * RD(cn2, 1, dn2) / 3 +
                    _k2 * Abs(cn) / dn :
                    // https://dlmf.nist.gov/19.25.E11
                    -_kp2 * sn2 * RD(dn2, 1, cn2) / 3 +
                    dn / Abs(cn))) :
                E();

            // Enforce usual trig-like symmetries
            if (SignBit(cn))
                ei = 2 * E() - ei;

            return CopySign(ei, sn);
        }

        /// <summary>
        /// The incomplete integral of the third kind in terms of Jacobi elliptic functions.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns>Π(φ, α^2, <i>k</i>) as though φ ∈ (-π, π].</returns>
        public double Pi(double sn, double cn, double dn)
        {
            // Carlson, eq. 4.7 and
            // https://dlmf.nist.gov/19.25.E14
            double
              cn2 = cn * cn, dn2 = dn * dn, sn2 = sn * sn,
              pii = cn2 != 0 ? Abs(sn) * (RF(cn2, dn2, 1) +
                                          _alpha2 * sn2 *
                                          RJ(cn2, dn2, 1, cn2 + _alphap2 * sn2) / 3) :
              Pi();

            // Enforce usual trig-like symmetries
            if (SignBit(cn))
                pii = 2 * Pi() - pii;

            return CopySign(pii, sn);
        }

        /// <summary>
        /// Jahnke's incomplete elliptic integral in terms of Jacobi elliptic functions.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns><i>D</i>(φ, <i>k</i>) as though φ ∈ (-π, π].</returns>
        public double D(double sn, double cn, double dn)
        {
            // Carlson, eq. 4.8 and
            // https://dlmf.nist.gov/19.25.E13
            double
              cn2 = cn * cn, dn2 = dn * dn, sn2 = sn * sn,
              di = cn2 != 0 ? Abs(sn) * sn2 * RD(cn2, dn2, 1) / 3 : D();

            // Enforce usual trig-like symmetries
            if (SignBit(cn))
                di = 2 * D() - di;

            return CopySign(di, sn);
        }

        /// <summary>
        /// Legendre's geodesic longitude integral in terms of Jacobi elliptic functions.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns><i>G</i>(φ, α^2, <i>k</i>) as though φ ∈ (-π, π].</returns>
        public double G(double sn, double cn, double dn)
        {
            double
              cn2 = cn * cn, dn2 = dn * dn, sn2 = sn * sn,
              gi = cn2 != 0 ? Abs(sn) * (RF(cn2, dn2, 1) +
                                         (_alpha2 - _k2) * sn2 *
                                         RJ(cn2, dn2, 1, cn2 + _alphap2 * sn2) / 3) :
              G();

            // Enforce usual trig-like symmetries
            if (SignBit(cn))
                gi = 2 * G() - gi;

            return CopySign(gi, sn);
        }

        /// <summary>
        /// Cayley's geodesic longitude difference integral in terms of Jacobi elliptic functions.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns><i>H</i>(φ, α^2, <i>k</i>) as though φ ∈ (-π, π].</returns>
        public double H(double sn, double cn, double dn)
        {
            double
              cn2 = cn * cn, dn2 = dn * dn, sn2 = sn * sn,
              // WARNING: large cancellation if k2 = 1, alpha2 = 0, and phi near pi/2
              hi = cn2 != 0 ? Abs(sn) * (RF(cn2, dn2, 1) -
                                         _alphap2 * sn2 *
                                         RJ(cn2, dn2, 1, cn2 + _alphap2 * sn2) / 3) :
              H();

            // Enforce usual trig-like symmetries
            if (SignBit(cn))
                hi = 2 * H() - hi;

            return CopySign(hi, sn);
        }

        /// <summary>
        /// The periodic incomplete integral of the first kind.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns>the periodic function π <i>F</i>(φ, <i>k</i>) / (2 <i>K</i>(<i>k</i>)) - φ.</returns>
        public double DeltaF(double sn, double cn, double dn)
        {
            // Function is periodic with period pi
            if (SignBit(cn)) { cn = -cn; sn = -sn; }
            return F(sn, cn, dn) * (PI / 2) / K() - Atan2(sn, cn);
        }

        /// <summary>
        /// The periodic incomplete integral of the second kind.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns>the periodic function π <i>E</i>(φ, <i>k</i>) / (2 <i>E</i>(<i>k</i>)) - φ.</returns>
        public double DeltaE(double sn, double cn, double dn)
        {
            // Function is periodic with period pi
            if (SignBit(cn)) { cn = -cn; sn = -sn; }
            return E(sn, cn, dn) * (PI / 2) / E() - Atan2(sn, cn);
        }

        /// <summary>
        /// The periodic inverse of the incomplete integral of the second kind.
        /// </summary>
        /// <param name="stau">sinτ</param>
        /// <param name="ctau">cosτ</param>
        /// <returns>the periodic function <i>E</i>^-1(τ (2 <i>E</i>(<i>k</i>)/π), <i>k</i>) - τ.</returns>
        public double DeltaEinv(double stau, double ctau)
        {
            // Function is periodic with period pi
            if (SignBit(ctau)) { ctau = -ctau; stau = -stau; }
            var tau = Atan2(stau, ctau);
            return Einv(tau * E() / (PI / 2)) - tau;
        }

        /// <summary>
        /// The periodic incomplete integral of the third kind.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns>the periodic function π Π(φ, α^2, <i>k</i>) / (2 Π(α^2, <i>k</i>)) - φ.</returns>
        public double DeltaPi(double sn, double cn, double dn)
        {
            // Function is periodic with period pi
            if (SignBit(cn)) { cn = -cn; sn = -sn; }
            return Pi(sn, cn, dn) * (PI / 2) / Pi() - Atan2(sn, cn);
        }

        /// <summary>
        /// The periodic Jahnke's incomplete elliptic integral.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns>the periodic function π <i>D</i>(φ, <i>k</i>) / (2 <i>D</i>(<i>k</i>)) - φ.</returns>
        public double DeltaD(double sn, double cn, double dn)
        {
            // Function is periodic with period pi
            if (SignBit(cn)) { cn = -cn; sn = -sn; }
            return D(sn, cn, dn) * (PI / 2) / D() - Atan2(sn, cn);
        }

        /// <summary>
        /// Legendre's periodic geodesic longitude integral.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns>the periodic function π <i>G</i>(φ, <i>k</i>) / (2 <i>G</i>(<i>k</i>)) - φ.</returns>
        public double DeltaG(double sn, double cn, double dn)
        {
            // Function is periodic with period pi
            if (SignBit(cn)) { cn = -cn; sn = -sn; }
            return G(sn, cn, dn) * (PI / 2) / G() - Atan2(sn, cn);
        }

        /// <summary>
        /// Cayley's periodic geodesic longitude difference integral.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <param name="dn">sqrt(1 - <i>k</i>^2 sin^2φ)</param>
        /// <returns>the periodic function π <i>H</i>(φ, <i>k</i>) / (2 <i>H</i>(<i>k</i>)) -  φ.</returns>
        public double DeltaH(double sn, double cn, double dn)
        {
            // Function is periodic with period pi
            if (SignBit(cn)) { cn = -cn; sn = -sn; }
            return H(sn, cn, dn) * (PI / 2) / H() - Atan2(sn, cn);
        }

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
        public void Sncndn(double x, out double sn, out double cn, out double dn)
        {
            // Bulirsch's sncndn routine, p 89.
            if (_kp2 != 0)
            {
                double mc = _kp2, d = 0;
                if (SignBit(_kp2))
                {
                    d = 1 - mc;
                    mc /= -d;
                    d = Sqrt(d);
                    x *= d;
                }

                var c = 0d;           // To suppress warning about uninitialized variable
                Span<double>
                    m = stackalloc double[num_],
                    n = stackalloc double[num_];

                var l = 0;
                for (var a = 1d; l < num_ || GEOGRAPHICLIB_PANIC; ++l)
                {
                    // This converges quadratically.  Max 5 trips
                    m[l] = a;
                    n[l] = mc = Sqrt(mc);
                    c = (a + mc) / 2;
                    if (!(Abs(a - mc) > tolJAC * a))
                    {
                        ++l;
                        break;
                    }
                    mc *= a;
                    a = c;
                }
                x *= c;
                sn = Sin(x);
                cn = Cos(x);
                dn = 1;
                if (sn != 0)
                {
                    var a = cn / sn;
                    c *= a;
                    while (l-- > 0)
                    {
                        var b = m[l];
                        a *= c;
                        c *= dn;
                        dn = (n[l] + a) / (b + a);
                        a = c / b;
                    }
                    a = 1 / Sqrt(c * c + 1);
                    sn = SignBit(sn) ? -a : a;
                    cn = c * sn;
                    if (SignBit(_kp2))
                    {
                        Swap(ref cn, ref dn);
                        sn /= d;
                    }
                }
            }
            else
            {
                sn = Tanh(x);
                dn = cn = 1 / Cosh(x);
            }
        }

        /// <summary>
        /// The Δ amplitude function.
        /// </summary>
        /// <param name="sn">sinφ</param>
        /// <param name="cn">cosφ</param>
        /// <returns>Δ = sqrt(1 - <i>k</i>^2 sin^2φ)</returns>
        public double Delta(double sn, double cn) => Sqrt(_k2 < 0 ? 1 - _k2 * sn * sn : _kp2 + _k2 * cn * cn);

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
        public static double RF(double x, double y, double z)
        {
            // Carlson, eqs 2.2 - 2.7
            double
              A0 = (x + y + z) / 3,
              An = A0,
              Q = Max(Max(Abs(A0 - x), Abs(A0 - y)), Abs(A0 - z)) / tolRF,
              x0 = x,
              y0 = y,
              z0 = z,
              mul = 1;

            while (Q >= mul * Abs(An))
            {
                // Max 6 trips
                var lam = Sqrt(x0) * Sqrt(y0) + Sqrt(y0) * Sqrt(z0) + Sqrt(z0) * Sqrt(x0);
                An = (An + lam) / 4;
                x0 = (x0 + lam) / 4;
                y0 = (y0 + lam) / 4;
                z0 = (z0 + lam) / 4;
                mul *= 4;
            }

            double
              X = (A0 - x) / (mul * An),
              Y = (A0 - y) / (mul * An),
              Z = -(X + Y),
              E2 = X * Y - Z * Z,
              E3 = X * Y * Z;

            // https://dlmf.nist.gov/19.36.E1
            // Polynomial is
            // (1 - E2/10 + E3/14 + E2^2/24 - 3*E2*E3/44
            //    - 5*E2^3/208 + 3*E3^2/104 + E2^2*E3/16)
            // convert to Horner form...
            return (E3 * (6930 * E3 + E2 * (15015 * E2 - 16380) + 17160) +
                    E2 * ((10010 - 5775 * E2) * E2 - 24024) + 240240) /
              (240240 * Sqrt(An));
        }

        /// <summary>
        /// Complete symmetric integral of the first kind, <i>Rf</i> with one argument zero.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns><i>Rf</i>(<i>x</i>, <i>y</i>, 0).</returns>
        /// <remarks>
        /// The arguments <paramref name="x"/> and <paramref name="y"/> must be positive.
        /// </remarks>
        public static double RF(double x, double y)
        {
            // Carlson, eqs 2.36 - 2.38
            double xn = Sqrt(x), yn = Sqrt(y);

            if (xn < yn) Swap(ref xn, ref yn);

            while (Abs(xn - yn) > tolRG0 * xn)
            {
                // Max 4 trips
                var t = (xn + yn) / 2;
                yn = Sqrt(xn * yn);
                xn = t;
            }

            return PI / (xn + yn);
        }

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
        public static double RC(double x, double y) =>
            // Defined only for y != 0 and x >= 0.
            (!(x >= y) ?   // x < y  and catch nans
                           // https://dlmf.nist.gov/19.2.E18
                 Atan(Sqrt((y - x) / x)) / Sqrt(y - x) :
                 (x == y ? 1 / Sqrt(y) :
                   Asinh(y > 0 ?
                          // https://dlmf.nist.gov/19.2.E19
                          // atanh(sqrt((x - y) / x))
                          Sqrt((x - y) / y) :
                          // https://dlmf.nist.gov/19.2.E20
                          // atanh(sqrt(x / (x - y)))
                          Sqrt(-x / y)) / Sqrt(x - y)));

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
        public static double RG(double x, double y, double z)
        {
            if (z == 0)
                Swap(ref y, ref z);

            // Carlson, eq 1.7
            return (z * RF(x, y, z) - (x - z) * (y - z) * RD(x, y, z) / 3
                    + Sqrt(x * y / z)) / 2;
        }

        /// <summary>
        /// Complete symmetric integral of the second kind, <i>Rg</i> with one argument zero.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns><i>Rg</i>(<i>x</i>, <i>y</i>, 0).</returns>
        /// <remarks>
        /// The arguments <paramref name="x"/> and <paramref name="y"/> must be positive.
        /// </remarks>
        public static double RG(double x, double y)
        {
            // Carlson, eqs 2.36 - 2.39
            double
              x0 = Sqrt(Max(x, y)),
              y0 = Sqrt(Min(x, y)),
              xn = x0,
              yn = y0,
              s = 0,
              mul = 0.25;

            while (Abs(xn - yn) > tolRG0 * xn)
            {
                // Max 4 trips
                var t = (xn + yn) / 2;
                yn = Sqrt(xn * yn);
                xn = t;
                mul *= 2;
                t = xn - yn;
                s += mul * t * t;
            }

            return (Sq((x0 + y0) / 2) - s) * PI / (2 * (xn + yn));
        }

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
        public static double RJ(double x, double y, double z, double p)
        {
            // Carlson, eqs 2.17 - 2.25
            double
              A0 = (x + y + z + 2 * p) / 5,
              An = A0,
              delta = (p - x) * (p - y) * (p - z),
              Q = Max(Max(Abs(A0 - x), Abs(A0 - y)), Max(Abs(A0 - z), Abs(A0 - p))) / tolRD,
              x0 = x,
              y0 = y,
              z0 = z,
              p0 = p,
              mul = 1,
              mul3 = 1,
              s = 0;

            while (Q >= mul * Abs(An))
            {
                // Max 7 trips
                double
                  lam = Sqrt(x0) * Sqrt(y0) + Sqrt(y0) * Sqrt(z0) + Sqrt(z0) * Sqrt(x0),
                  d0 = (Sqrt(p0) + Sqrt(x0)) * (Sqrt(p0) + Sqrt(y0)) * (Sqrt(p0) + Sqrt(z0)),
                  e0 = delta / (mul3 * Sq(d0));

                s += RC(1, 1 + e0) / (mul * d0);
                An = (An + lam) / 4;
                x0 = (x0 + lam) / 4;
                y0 = (y0 + lam) / 4;
                z0 = (z0 + lam) / 4;
                p0 = (p0 + lam) / 4;
                mul *= 4;
                mul3 *= 64;
            }

            double
              X = (A0 - x) / (mul * An),
              Y = (A0 - y) / (mul * An),
              Z = (A0 - z) / (mul * An),
              P = -(X + Y + Z) / 2,
              E2 = X * Y + X * Z + Y * Z - 3 * P * P,
              E3 = X * Y * Z + 2 * P * (E2 + 2 * P * P),
              E4 = (2 * X * Y * Z + P * (E2 + 3 * P * P)) * P,
              E5 = X * Y * Z * P * P;

            // https://dlmf.nist.gov/19.36.E2
            // Polynomial is
            // (1 - 3*E2/14 + E3/6 + 9*E2^2/88 - 3*E4/22 - 9*E2*E3/52 + 3*E5/26
            //    - E2^3/16 + 3*E3^2/40 + 3*E2*E4/20 + 45*E2^2*E3/272
            //    - 9*(E3*E4+E2*E5)/68)
            return ((471240 - 540540 * E2) * E5 +
                    (612612 * E2 - 540540 * E3 - 556920) * E4 +
                    E3 * (306306 * E3 + E2 * (675675 * E2 - 706860) + 680680) +
                    E2 * ((417690 - 255255 * E2) * E2 - 875160) + 4084080) /
              (4084080 * mul * An * Sqrt(An)) + 6 * s;
        }

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
        public static double RD(double x, double y, double z)
        {
            // Carlson, eqs 2.28 - 2.34
            double
              A0 = (x + y + 3 * z) / 5,
              An = A0,
              Q = Max(Max(Abs(A0 - x), Abs(A0 - y)), Abs(A0 - z)) / tolRD,
              x0 = x,
              y0 = y,
              z0 = z,
              mul = 1,
              s = 0;

            while (Q >= mul * Abs(An))
            {
                // Max 7 trips
                var lam = Sqrt(x0) * Sqrt(y0) + Sqrt(y0) * Sqrt(z0) + Sqrt(z0) * Sqrt(x0);
                s += 1 / (mul * Sqrt(z0) * (z0 + lam));
                An = (An + lam) / 4;
                x0 = (x0 + lam) / 4;
                y0 = (y0 + lam) / 4;
                z0 = (z0 + lam) / 4;
                mul *= 4;
            }

            double
              X = (A0 - x) / (mul * An),
              Y = (A0 - y) / (mul * An),
              Z = -(X + Y) / 3,
              E2 = X * Y - 6 * Z * Z,
              E3 = (3 * X * Y - 8 * Z * Z) * Z,
              E4 = 3 * (X * Y - Z * Z) * Z * Z,
              E5 = X * Y * Z * Z * Z;

            // https://dlmf.nist.gov/19.36.E2
            // Polynomial is
            // (1 - 3*E2/14 + E3/6 + 9*E2^2/88 - 3*E4/22 - 9*E2*E3/52 + 3*E5/26
            //    - E2^3/16 + 3*E3^2/40 + 3*E2*E4/20 + 45*E2^2*E3/272
            //    - 9*(E3*E4+E2*E5)/68)
            return ((471240 - 540540 * E2) * E5 +
                    (612612 * E2 - 540540 * E3 - 556920) * E4 +
                    E3 * (306306 * E3 + E2 * (675675 * E2 - 706860) + 680680) +
                    E2 * ((417690 - 255255 * E2) * E2 - 875160) + 4084080) /
              (4084080 * mul * An * Sqrt(An)) + 3 * s;
        }

        #endregion

        /// <summary>
        /// Reset the modulus and parameter.
        /// </summary>
        /// <param name="k2">the new value of square of the modulus <i>k</i>^2 which must lie in (-∞, 1].</param>
        /// <param name="alpha2">the new value of parameter α^2. α^2 must lie in (-∞, 1].</param>
        public void Reset(double k2, double alpha2) => Reset(k2, alpha2, 1 - k2, 1 - alpha2);

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
        public void Reset(double k2, double alpha2, double kp2, double alphap2)
        {
            // Accept nans here (needed for GeodesicExact)
            if (k2 > 1)
                throw new GeographicException("Parameter k2 is not in (-inf, 1]");
            if (alpha2 > 1)
                throw new GeographicException("Parameter alpha2 is not in (-inf, 1]");

            // TODO: NaNs are allowed below as SignBit alone is failing GeodSolve94_CheckFixFor_lat2_eq_NaN_BeingTreatedAs_lat2_eq_Zero_Exact.
            if (kp2 < 0) 
                throw new GeographicException("Parameter kp2 is not in [0, inf)");
            if (alphap2 < 0)
                throw new GeographicException("Parameter alphap2 is not in [0, inf)");

            _k2 = k2;
            _kp2 = kp2;
            _alpha2 = alpha2;
            _alphap2 = alphap2;
            _eps = _k2 / Sq(Sqrt(_kp2) + 1);

            // Values of complete elliptic integrals for k = 0,1 and alpha = 0,1
            //         K     E     D
            // k = 0:  pi/2  pi/2  pi/4
            // k = 1:  inf   1     inf
            //                    Pi    G     H
            // k = 0, alpha = 0:  pi/2  pi/2  pi/4
            // k = 1, alpha = 0:  inf   1     1
            // k = 0, alpha = 1:  inf   inf   pi/2
            // k = 1, alpha = 1:  inf   inf   inf
            //
            // Pi(0, k) = K(k)
            // G(0, k) = E(k)
            // H(0, k) = K(k) - D(k)
            // Pi(0, k) = K(k)
            // G(0, k) = E(k)
            // H(0, k) = K(k) - D(k)
            // Pi(alpha2, 0) = pi/(2*sqrt(1-alpha2))
            // G(alpha2, 0) = pi/(2*sqrt(1-alpha2))
            // H(alpha2, 0) = pi/(2*(1 + sqrt(1-alpha2)))
            // Pi(alpha2, 1) = inf
            // H(1, k) = K(k)
            // G(alpha2, 1) = H(alpha2, 1) = RC(1, alphap2)

            if (_k2 != 0)
            {
                // Complete elliptic integral K(k), Carlson eq. 4.1
                // https://dlmf.nist.gov/19.25.E1
                _Kc = _kp2 != 0 ? RF(_kp2, 1) : double.PositiveInfinity;
                // Complete elliptic integral E(k), Carlson eq. 4.2
                // https://dlmf.nist.gov/19.25.E1
                _Ec = _kp2 != 0 ? 2 * RG(_kp2, 1) : 1;
                // D(k) = (K(k) - E(k))/k^2, Carlson eq.4.3
                // https://dlmf.nist.gov/19.25.E1
                _Dc = _kp2 != 0 ? RD(0, _kp2, 1) / 3 : double.PositiveInfinity;
            }
            else
            {
                _Kc = _Ec = PI / 2; _Dc = _Kc / 2;
            }

            if (_alpha2 != 0)
            {
                // https://dlmf.nist.gov/19.25.E2
                double
                    rj = (_kp2 != 0 && _alphap2 != 0) ? RJ(0, _kp2, 1, _alphap2) : double.PositiveInfinity,
                    // Only use rc if _kp2 = 0.
                    rc = _kp2 != 0 ? 0 : (_alphap2 != 0 ? RC(1, _alphap2) : double.PositiveInfinity);
                // Pi(alpha^2, k)
                _Pic = _kp2 != 0 ? _Kc + _alpha2 * rj / 3 : double.PositiveInfinity;
                // G(alpha^2, k)
                _Gc = _kp2 != 0 ? _Kc + (_alpha2 - _k2) * rj / 3 : rc;
                // H(alpha^2, k)
                _Hc = _kp2 != 0 ? _Kc - (_alphap2 != 0 ? _alphap2 * rj : 0) / 3 : rc;
            }
            else
            {
                _Pic = _Kc; _Gc = _Ec;
                // Hc = Kc - Dc but this involves large cancellations if k2 is close to
                // 1.  So write (for alpha2 = 0)
                //   Hc = int(cos(phi)^2/sqrt(1-k2*sin(phi)^2),phi,0,pi/2)
                //      = 1/sqrt(1-k2) * int(sin(phi)^2/sqrt(1-k2/kp2*sin(phi)^2,...)
                //      = 1/kp * D(i*k/kp)
                // and use D(k) = RD(0, kp2, 1) / 3
                // so Hc = 1/kp * RD(0, 1/kp2, 1) / 3
                //       = kp2 * RD(0, 1, kp2) / 3
                // using https://dlmf.nist.gov/19.20.E18
                // Equivalently
                //   RF(x, 1) - RD(0, x, 1)/3 = x * RD(0, 1, x)/3 for x > 0
                // For k2 = 1 and alpha2 = 0, we have
                //   Hc = int(cos(phi),...) = 1
                _Hc = _kp2 != 0 ? _kp2 * RD(0, 1, _kp2) / 3 : 1;
            }
        }


    }
}
