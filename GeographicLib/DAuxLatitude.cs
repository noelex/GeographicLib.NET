using System;
using static GeographicLib.MathEx;
using static System.Math;

namespace GeographicLib
{

    /// <summary>
    /// Divided differences of auxiliary latitudes.
    /// </summary>
    /// <remarks>
    /// This class computed the divided differences of auxiliary latitudes and
    /// some other divided differences needed to support rhumb line calculations.
    /// </remarks>
    public class DAuxLatitude : AuxLatitude
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="DAuxLatitude"/> class
        /// with the equatorial radius and flattening. 
        /// </summary>
        /// <param name="a">Equatorial radius.</param>
        /// <param name="f">
        /// Flattening of ellipsoid. Setting <paramref name="f"/> = 0 gives a sphere.
        /// Negative <paramref name="f"/> gives a prolate ellipsoid.
        /// </param>
        /// <exception cref="GeographicException"/>
        public DAuxLatitude(double a, double f) : base(a, f)
        {
        }

        /// <summary>
        /// The divided difference of one auxiliary latitude with respect to another.
        /// </summary>
        /// <param name="auxin">The type of auxiliary latitude <i>zeta</i>.</param>
        /// <param name="auxout">The type of auxiliary latitude <i>eta</i>.</param>
        /// <param name="zeta1">The first of the input auxiliary latitudes.</param>
        /// <param name="zeta2">The second of the input auxiliary latitudes.</param>
        /// <returns>The divided difference (<i>eta2</i> - <i>eta1</i>) / (<i>zeta2</i> - <i>zeta1</i>).</returns>
        /// <remarks>
        /// <para>This routine uses the series method.</para>
        ///
        /// <para>
        /// In the expression for the divided difference above, the angle quantities
        /// should be understood as the conventional measure of angle(either in
        /// radians or in degrees).
        /// </para>
        ///
        /// <para>
        /// The Fourier coefficients for a specific <paramref name="auxin"/> and <paramref name="auxout"/> are
        /// computed and saved on the first call; the saved coefficients are used on
        /// subsequent calls.The series method is accurate for abs(<i>f</i>) ≤ 1/150.
        /// </para>
        /// </remarks>
        public double DConvert(AuxLatitudeType auxin, AuxLatitudeType auxout,
                            AuxAngle zeta1, AuxAngle zeta2)
        {
            int k = ind(auxout, auxin);
            if (k < 0) return double.NaN;
            if (auxin == auxout) return 1;
            if (double.IsNaN(_c[Lmax * (k + 1) - 1]))
                fillcoeff(auxin, auxout, k);
            AuxAngle zeta1n = zeta1.Normalized(), zeta2n = zeta2.Normalized();
            return 1 + DClenshaw(true, zeta2n.Radians - zeta1n.Radians,
                                 zeta1n.Y, zeta1n.X, zeta2n.Y, zeta2n.X,
                                  _c.AsSpan().Slice(Lmax * k), Lmax);
        }

        /// <summary>
        /// The divided difference of the parametric latitude with respect to the geographic latitude.
        /// </summary>
        /// <param name="phi1">The first geographic latitude as an <see cref="AuxAngle"/>.</param>
        /// <param name="phi2">The second geographic latitude as an <see cref="AuxAngle"/>.</param>
        /// <returns>
        /// The divided difference (<i>beta2</i> - <i>beta1</i>) / (<i>phi2</i> - <i>phi1</i>), where <i>beta</i> is the parametric latitude.
        /// </returns>
        /// <remarks>
        /// This routine uses the exact formulas and is valid for arbitrary latitude.
        /// </remarks>
        public double DParametric(AuxAngle phi1, AuxAngle phi2)
        {
            double tx = phi1.Tan, ty = phi2.Tan, r;
            // DbetaDphi = Datan(fm1*tx, fm1*ty) * fm1 / Datan(tx, ty)
            // Datan(x, y) = 1/(1 + x^2),                       for x = y
            //             = (atan(y) - atan(x)) / (y-x),       for x*y < 0
            //             = atan( (y-x) / (1 + x*y) ) / (y-x), for x*y > 0
            if (!(tx * ty >= 0))        // This includes, e.g., tx = 0, ty = inf
                r = (Atan(_fm1 * ty) - Atan(_fm1 * tx)) /
                  (Atan(ty) - Atan(tx));
            else if (tx == ty)
            {        // This includes the case tx = ty = inf
                tx *= tx;
                if (tx <= 1)
                    r = _fm1 * (1 + tx) / (1 + _e2m1 * tx);
                else
                {
                    tx = 1 / tx;
                    r = _fm1 * (1 + tx) / (_e2m1 + tx);
                }
            }
            else
            {
                if (tx * ty <= 1)
                    r = Atan2(_fm1 * (ty - tx), 1 + _e2m1 * tx * ty)
                      / Atan2(ty - tx, 1 + tx * ty);
                else
                {
                    tx = 1 / tx; ty = 1 / ty;
                    r = Atan2(_fm1 * (ty - tx), _e2m1 + tx * ty)
                      / Atan2(ty - tx, 1 + tx * ty);
                }
            }
            return r;
        }

        /// <summary>
        /// The divided difference of the rectifying latitude with respect to the geographic latitude.
        /// </summary>
        /// <param name="phi1">The first geographic latitude as an <see cref="AuxAngle"/>.</param>
        /// <param name="phi2">The second geographic latitude as an <see cref="AuxAngle"/>.</param>
        /// <returns>
        /// The divided difference (<i>mu2</i> - <i>mu1</i>) / (<i>phi2</i> - <i>phi1</i>), where <i>mu</i> is the rectifying latitude.
        /// </returns>
        /// <remarks>
        /// This routine uses the exact formulas and is valid for arbitrary latitude.
        /// </remarks>
        public double DRectifying(AuxAngle phi1, AuxAngle phi2)
        {
            // Stipulate that phi1 and phi2 are in [-90d, 90d]
            double x = phi1.Radians, y = phi2.Radians;
            if (x == y)
            {
                double d;
                AuxAngle mu1 = Rectifying(phi1, out d);
                double tphi1 = phi1.Tan, tmu1 = mu1.Tan;
                return
                  IsFinite(tphi1) ? d * Sq(sc(tphi1) / sc(tmu1)) : 1 / d;
            }
            else if (x * y < 0)
                return (Rectifying(phi2).Radians -
                        Rectifying(phi1).Radians) / (y - x);
            else
            {
                AuxAngle bet1 = Parametric(phi1), bet2 = Parametric(phi2);
                double dEdbet = DE(bet1, bet2), dbetdphi = DParametric(phi1, phi2);
                return _b * dEdbet / RectifyingRadius(true) * dbetdphi;
            }
        }

        /// <summary>
        /// The divided difference of the isometric latitude with respect to the geographic latitude.
        /// </summary>
        /// <param name="phi1">The first geographic latitude as an <see cref="AuxAngle"/>.</param>
        /// <param name="phi2">The second geographic latitude as an <see cref="AuxAngle"/>.</param>
        /// <returns>
        /// The divided difference (<i>psi2</i> - <i>psi1</i>) / (<i>phi2</i> - <i>phi1</i>), 
        /// where <i>psi</i> = asinh(tan(<i>chi</i>)) is the isometric latitude
        /// and <i>chi</i> is the conformal latitude.
        /// </returns>
        /// <remarks>
        /// This routine uses the exact formulas and is valid for arbitrary latitude.
        /// </remarks>
        public double DIsometric(AuxAngle phi1, AuxAngle phi2)
        {
            // psi = asinh(tan(phi)) - e^2 * atanhee(tan(phi))
            double tphi1 = phi1.Tan, tphi2 = phi2.Tan;
            return double.IsNaN(tphi1) || double.IsNaN(tphi2) ? double.NaN :
              (double.IsInfinity(tphi1) || double.IsInfinity(tphi2) ? double.PositiveInfinity :
               (Dasinh(tphi1, tphi2) - _e2 * Datanhee(tphi1, tphi2)) /
               Datan(tphi1, tphi2));
        }

        /// <summary>
        /// The divided difference of <see cref="AuxLatitude.Clenshaw(bool, double, double, ReadOnlySpan{double}, int)"/>.
        /// </summary>
        /// <param name="sinp">If <see langword="true"/> sum the sine series, else sum the cosine series.</param>
        /// <param name="Delta">Either 1 <b>or</b> (<i>zeta2</i> - <i>zeta1</i>) in radians.</param>
        /// <param name="szeta1">sin(<i>zeta1</i>).</param>
        /// <param name="czeta1">cos(<i>zeta1</i>).</param>
        /// <param name="szeta2">sin(<i>zeta2</i>).</param>
        /// <param name="czeta2">cos(<i>zeta2</i>).</param>
        /// <param name="c">The array of coefficients.</param>
        /// <param name="K">The number of coefficients.</param>
        /// <returns>
        /// The divided difference.
        /// <para>
        /// The result is
        /// <c>
        ///    (AuxLatitude.Clenshaw(sinp, szeta2, czeta2, c, K) -
        ///      AuxLatitude.Clenshaw(sinp, szeta1, czeta1, c, K) ) / Delta
        /// </c>
        /// </para>
        /// </returns>
        /// <remarks>
        /// <paramref name="Delta"/> **must** be either 1 or (<i>zeta2</i> - <i>zeta1</i>);
        /// other values will return nonsense.
        /// </remarks>
        public static double DClenshaw(bool sinp, double Delta,
                                    double szeta1, double czeta1,
                                    double szeta2, double czeta2,
                                    ReadOnlySpan<double> c, int K)
        {
            // Evaluate
            // (Clenshaw(sinp, szeta2, czeta2, c, K) -
            //  Clenshaw(sinp, szeta1, czeta1, c, K)) / Delta
            // or
            // sum(c[k] * (sin( (2*k+2) * zeta2) - sin( (2*k+2) * zeta2)), i, 0, K-1)
            //   / Delta
            // (if !sinp, then change sin->cos here.)
            //
            // Delta is EITHER 1, giving the plain difference OR (zeta2 - zeta1) in
            // radians, giving the divided difference.  Other values will give
            // nonsense.
            //
            int k = K;
            // suffices a b denote [1,1], [2,1] elements of matrix/vector
            double D2 = Delta * Delta,
              czetap = czeta2 * czeta1 - szeta2 * szeta1,
              szetap = szeta2 * czeta1 + czeta2 * szeta1,
              czetam = czeta2 * czeta1 + szeta2 * szeta1,
              // sin(zetam) / Delta
              szetamd = (Delta == 1 ? szeta2 * czeta1 - czeta2 * szeta1 :
                         (Delta != 0 ? Sin(Delta) / Delta : 1)),
              Xa = 2 * czetap * czetam,
              Xb = -2 * szetap * szetamd,
              u0a = 0, u0b = 0, u1a = 0, u1b = 0; // accumulators for sum
            for (--k; k >= 0; --k)
            {
                // temporary real = X . U0 - U1 + c[k] * I
                double ta = Xa * u0a + D2 * Xb * u0b - u1a + c[k],
                  tb = Xb * u0a + Xa * u0b - u1b;
                // U1 = U0; U0 = real
                u1a = u0a; u0a = ta;
                u1b = u0b; u0b = tb;
            }
            // P = U0 . F[0] - U1 . F[-1]
            // if sinp:
            //   F[0] = [ sin(2*zeta2) + sin(2*zeta1),
            //           (sin(2*zeta2) - sin(2*zeta1)) / Delta]
            //        = 2 * [ szetap * czetam, czetap * szetamd ]
            //   F[-1] = [0, 0]
            // else:
            //   F[0] = [ cos(2*zeta2) + cos(2*zeta1),
            //           (cos(2*zeta2) - cos(2*zeta1)) / Delta]
            //        = 2 * [ czetap * czetam, -szetap * szetamd ]
            //   F[-1] = [2, 0]
            double F0a = (sinp ? szetap : czetap) * czetam,
              F0b = (sinp ? czetap : -szetap) * szetamd,
              Fm1a = sinp ? 0 : 1;  // Fm1b = 0;
                                    // Don't both to compute sum...
                                    // divided difference (or difference if Delta == 1)
            return 2 * (F0a * u0b + F0b * u0a - Fm1a * u1b);
        }

        /// <summary>
        /// The divided difference of the isometric latitude with respect to the conformal latitude.
        /// </summary>
        /// <param name="x">tan(<i>chi1</i>).</param>
        /// <param name="y">tan(<i>chi2</i>).</param>
        /// <returns>
        /// The divided difference (<i>psi2</i> - <i>psi1</i>) / (<i>chi2</i> - <i>chi1</i>),
        /// where <i>psi</i> = asinh(tan(<i>chi</i>)).
        /// </returns>
        /// <remarks>
        /// The parameters for this routine are the <b>tangents</b> of conformal latitude.
        /// <para>
        /// This routine computes Dasinh(x, y) / Datan(x, y).
        /// </para>
        /// </remarks>
        public static double Dlam(double x, double y)
        {
            return x == y ? sc(x) :
              (double.IsNaN(x) || double.IsNaN(y) ? double.NaN :
               (double.IsInfinity(x) || double.IsInfinity(y) ? double.PositiveInfinity :
                Dasinh(x, y) / Datan(x, y)));
        }

        /// <summary>
        /// The divided difference of the spherical rhumb area term with respect to the isometric latitude.
        /// </summary>
        /// <param name="x">tan(<i>chi1</i>).</param>
        /// <param name="y">tan(<i>chi2</i>).</param>
        /// <returns>
        /// The divided difference (p0(<i>chi2</i>) - p0(<i>chi1</i>)) / (<i>psi2</i> - <i>psi1</i>), 
        /// where p0(<i>chi</i>) = log(sec(<i>chi</i>)) and <i>psi</i> = asinh(tan(<i>chi</i>)).
        /// </returns>
        /// <remarks>
        /// This parameters for this routine are the <b>tangents</b> of conformal latitude.
        /// </remarks>
        public static double Dp0Dpsi(double x, double y)
        {
            return x == y ? sn(x) :
              (double.IsNaN(x + y) ? x + y : // N.B. nan for inf-inf
               (double.IsInfinity(x) ? CopySign(1.0, x) :
                (double.IsInfinity(y) ? CopySign(1.0, y) :
                 Dasinh(h(x), h(y)) * Dh(x, y) / Dasinh(x, y))));
        }

        internal static double Dsn(double x, double y)
        {
            double sc1 = sc(x);
            if (x == y) return 1 / (sc1 * (1 + x * x));
            double sc2 = sc(y), sn1 = sn(x), sn2 = sn(y);
            return x * y > 0 ?
              (sn1 / sc2 + sn2 / sc1) / ((sn1 + sn2) * sc1 * sc2) :
              (sn2 - sn1) / (y - x);
        }

        internal static double Datan(double x, double y)
        {
            double d = y - x, xy = x * y;
            return x == y ? 1 / (1 + xy) :
              (double.IsInfinity(xy) && xy > 0 ? 0 :
               (2 * xy > -1 ? Atan(d / (1 + xy)) : Atan(y) - Atan(x)) / d);
        }

        internal static double Dasinh(double x, double y)
        {
            double d = y - x, xy = x * y, hx = sc(x), hy = sc(y);
            // KF formula for x*y < 0 is asinh(y*hx - x*hy) / (y - x)
            // but this has problem if x*y overflows to -inf
            return x == y ? 1 / hx :
              (double.IsInfinity(d) ? 0 :
               (xy > 0 ? Asinh(d * (x * y < 1 ? (x + y) / (x * hy + y * hx) :
                                    (1 / x + 1 / y) / (hy / y + hx / x))) :
                Asinh(y) - Asinh(x)) / d);
        }

        // h(tan(x)) = tan(x) * sin(x) / 2
        internal static double h(double x) { return x * sn(x) / 2; }

        internal static double Dh(double x, double y)
        {
            if (double.IsNaN(x + y))
                return x + y;           // N.B. nan for inf-inf
            if (double.IsInfinity(x))
                return CopySign(1 / (double)(2), x);
            if (double.IsInfinity(y))
                return CopySign(1 / (double)(2), y);
            double sx = sn(x), sy = sn(y), d = sx * x + sy * y;
            if (d / 2 == 0)
                return (x + y) / 2;     // Handle underflow
            if (x * y <= 0)
                return (h(y) - h(x)) / (y - x); // Does not include x = y = 0
            double scx = sc(x), scy = sc(y);
            return ((x + y) / (2 * d)) *
              (Sq(sx * sy) + Sq(sy / scx) + Sq(sx / scy));
        }

        internal double Datanhee(double x, double y)
        {
            // atan(e*sn(tphi))/e:
            //  Datan(e*sn(x),e*sn(y))*Dsn(x,y)/Datan(x,y)
            // asinh(e1*sn(fm1*tphi)):
            //  Dasinh(e1*sn(fm1*x)), e1*sn(fm1*y)) *
            // e1 * Dsn(fm1*x, fm1*y) *fm1 / (e * Datan(x,y))
            // = Dasinh(e1*sn(fm1*x)), e1*sn(fm1*y)) *
            //  Dsn(fm1*x, fm1*y) / Datan(x,y)
            return _f < 0 ?
              Datan(_e * sn(x), _e * sn(y)) * Dsn(x, y) :
              Dasinh(_e1 * sn(_fm1 * x),
                     _e1 * sn(_fm1 * y)) *
              Dsn(_fm1 * x, _fm1 * y);
        }

        private static double Dsin(double x, double y)
        {
            double d = (x - y) / 2;
            return Cos((x + y) / 2) * (d != 0 ? Sin(d) / d : 1);
        }

        // (E(x) - E(y)) / (x - y)
        private double DE(AuxAngle X, AuxAngle Y)
        {
            AuxAngle Xn = X.Normalized(), Yn = Y.Normalized();
            // We assume that X and Y are in [-90d, 90d] and have the same sign
            // If not we would include
            //    if (Xn.y() * Yn.y() < 0)
            //      return d != 0 ? (E(X) - E(Y)) / d : 1;

            // The general formula fails for x = y = 0d and x = y = 90d.  Probably this
            // is fixable (the formula works for other x = y.  But let's also stipulate
            // that x != y .

            // Make both positive, so we can do the swap a <-> b trick
            Xn = Xn.CopyWith(y: Abs(Xn.Y));
            Yn = Yn.CopyWith(y: Abs(Yn.Y));

            double x = Xn.Radians, y = Yn.Radians, d = y - x,
              sx = Xn.Y, sy = Yn.Y, cx = Xn.X, cy = Yn.X,
              k2;
            // Switch prolate to oblate; we then can use the formulas for k2 < 0
            if (false && _f < 0)
            {
                d = -d; Swap(ref sx, ref cx); Swap(ref sy, ref cy);
                k2 = _e2;
            }
            else
            {
                k2 = -_e12;
            }
            // See DLMF: Eqs (19.11.2) and (19.11.4) letting
            // theta -> x, phi -> -y, psi -> z
            //
            // (E(y) - E(x)) / d = E(z)/d - k2 * sin(x) * sin(y) * sin(z)/d
            //                   = (E(z)/sin(z) - k2 * sin(x) * sin(y)) * sin(z)/d
            // tan(z/2) = (sin(x)*Delta(y) - sin(y)*Delta(x)) / (cos(x) + cos(y))
            //          = d * Dsin(x,y) * (sin(x) + sin(y))/(cos(x) + cos(y)) /
            //             (sin(x)*Delta(y) + sin(y)*Delta(x))
            //          = t = d * Dt
            // Delta(x) = sqrt(1 - k2 * sin(x)^2)
            // sin(z) = 2*t/(1+t^2); cos(z) = (1-t^2)/(1+t^2)
            double Dt = Dsin(x, y) * (sx + sy) /
              ((cx + cy) * (sx * Sqrt(1 - k2 * sy * sy) + sy * Sqrt(1 - k2 * sx * sx))),
              t = d * Dt, Dsz = 2 * Dt / (1 + t * t),
              sz = d * Dsz, cz = (1 - t) * (1 + t) / (1 + t * t),
              sz2 = sz * sz, cz2 = cz * cz, dz2 = 1 - k2 * sz2,
              // E(z)/sin(z)
              Ezbsz = (EllipticFunction.RF(cz2, dz2, 1)
                       - k2 * sz2 * EllipticFunction.RD(cz2, dz2, 1) / 3);
            return (Ezbsz - k2 * sx * sy) * Dsz;
        }
    }
}