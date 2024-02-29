using System;
using static GeographicLib.MathEx;
using static System.Math;

namespace GeographicLib
{
    partial class GeodesicLineExact
    {
        internal unsafe struct Priv
        {
            private int _nC4;

            private static readonly double tiny_ = GeodesicExact.tiny_;

            public double _lat1, _lon1, _azi1;
            public double _a, _f, _b, _c2, _f1, _e2, _salp0, _calp0, _k2,
              _salp1, _calp1, _ssig1, _csig1, _dn1, _stau1, _ctau1,
              _somg1, _comg1, _cchi1,
              _A4, _B41, _E0, _D0, _H0, _E1, _D1, _H1;

            public EllipticFunction.Priv _E;
            public GeodesicFlags _caps;

            public double _a13, _s13;

            // maximum possible size of __C4a is 4096, see GeodesicExact._nC4.
            // all the elements of _C4a are used
            private fixed double __C4a[4096];

            public void Init(GeodesicExact g,
                    double lat1, double lon1,
                    double azi1, double salp1, double calp1,
                    GeodesicFlags caps)
            {
                _E.Reset();
                LineInit(g, lat1, lon1, azi1, salp1, calp1, caps);
            }

            /// <summary>
            /// Constructor for a exact geodesic line staring at latitude <i>lat1</i>, longitude <i>lon1</i>, 
            /// and azimuth <i>azi1</i> (all in degrees).
            /// </summary>
            /// <param name="g">
            /// A <see cref="GeodesicExact"/> object used to compute the necessary information about the <see cref="GeodesicLineExact"/>.
            /// </param>
            /// <param name="lat1">latitude of point 1 (degrees).</param>
            /// <param name="lon1">longitude of point 1 (degrees).</param>
            /// <param name="azi1">azimuth at point 1 (degrees).</param>
            /// <param name="caps">bitor'ed combination of <see cref="GeodesicFlags"/> values specifying the capabilities
            /// the <see cref="GeodesicLineExact"/> object should possess, i.e., which quantities can be returned in calls to
            /// <see cref="GeodesicLineBase.Position(double, out double, out double, out double, out double, out double, out double, out double)"/>.</param>
            /// <remarks>
            /// <i>lat1</i> should be in the range [−90°, 90°].
            /// <para>
            /// The <see cref="GeodesicFlags"/> values possible for <i>caps</i> are
            /// <list type="bullet">
            /// <item><i>caps</i> |= <see cref="GeodesicFlags.Latitude"/> for the latitude <i>lat2</i>; this is added automatically;</item>
            /// <item><i>caps</i> |= <see cref="GeodesicFlags.Longitude"/> for the latitude <i>lon2</i>;</item>
            /// <item><i>caps</i> |= <see cref="GeodesicFlags.Azimuth"/> for the latitude <i>azi2</i>;this is added automatically;</item>
            /// <item><i>caps</i> |= <see cref="GeodesicFlags.Distance"/> for the distance <i>s12</i>;</item>
            /// <item><i>caps</i> |= <see cref="GeodesicFlags.ReducedLength"/> for the reduced length <i>m12</i>;</item>
            /// <item><i>caps</i> |= <see cref="GeodesicFlags.GeodesicScale"/> for the geodesic scales <i>M12</i> and <i>M21</i>;</item>
            /// <item><i>caps</i> |= <see cref="GeodesicFlags.Area"/> for the area <i>S12</i>;</item>
            /// <item><i>caps</i> |= <see cref="GeodesicFlags.DistanceIn"/> 
            /// permits the length of the geodesic to be given in terms of <i>s12</i>;
            /// without this capability the length can only be specified in terms of arc length;</item>
            /// <item><i>caps</i> |= <see cref="GeodesicFlags.All"/> for all of the above.</item>
            /// </list>
            /// The default value of <i>caps</i> is <see cref="GeodesicFlags.All"/>.
            /// </para>
            /// <para>
            /// If the point is at a pole, the azimuth is defined by keeping lon1 fixed,
            /// writing <i>lat1</i> = ±(90° − ε), and taking the limit ε → 0+.
            /// </para>
            /// </remarks>
            public void Init(GeodesicExact g,
                                     double lat1, double lon1, double azi1,
                                     GeodesicFlags caps)
            {
                _E.Reset();
                azi1 = AngNormalize(azi1);
                // Guard against underflow in salp0.  Also -0 is converted to +0.
                SinCosd(AngRound(azi1), out var salp1, out var calp1);
                LineInit(g, lat1, lon1, azi1, salp1, calp1, caps);
            }

            /// <inheritdoc/>
            public double GenPosition(bool arcmode, double s12_a12, GeodesicFlags outmask, out double lat2,
                out double lon2, out double azi2, out double s12, out double m12, out double M12, out double M21, out double S12)
            {
                lat2 = lon2 = azi2 = s12 = m12 = M12 = M21 = S12 = double.NaN;

                outmask &= _caps.Flags();
                if (!(arcmode || _caps.HasAny(GeodesicFlags.DistanceIn.Flags())))
                    // Uninitialized or impossible distance calculation requested
                    return double.NaN;

                // Avoid warning about uninitialized B12.
                double sig12, ssig12, csig12, E2 = 0, AB1 = 0;
                if (arcmode)
                {
                    // Interpret s12_a12 as spherical arc length
                    sig12 = s12_a12 * Degree;
                    SinCosd(s12_a12, out ssig12, out csig12);
                }
                else
                {
                    // Interpret s12_a12 as distance
                    double
                      tau12 = s12_a12 / (_b * _E0),
                      s = Sin(tau12),
                      c = Cos(tau12);
                    // tau2 = tau1 + tau12
                    E2 = -_E.DeltaEinv(_stau1 * c + _ctau1 * s, _ctau1 * c - _stau1 * s);
                    sig12 = tau12 - (E2 - _E1);
                    ssig12 = Sin(sig12);
                    csig12 = Cos(sig12);
                }

                double ssig2, csig2, sbet2, cbet2, salp2, calp2;
                // sig2 = sig1 + sig12
                ssig2 = _ssig1 * csig12 + _csig1 * ssig12;
                csig2 = _csig1 * csig12 - _ssig1 * ssig12;
                var dn2 = _E.Delta(ssig2, csig2);
                if (outmask.HasAny(GeodesicFlags.Distance | GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale))
                {
                    if (arcmode)
                    {
                        E2 = _E.DeltaE(ssig2, csig2, dn2);
                    }
                    AB1 = _E0 * (E2 - _E1);
                }
                // sin(bet2) = cos(alp0) * sin(sig2)
                sbet2 = _calp0 * ssig2;
                // Alt: cbet2 = hypot(csig2, salp0 * ssig2);
                cbet2 = Hypot(_salp0, _calp0 * csig2);
                if (cbet2 == 0)
                    // I.e., salp0 = 0, csig2 = 0.  Break the degeneracy in this case
                    cbet2 = csig2 = tiny_;
                // tan(alp0) = cos(sig2)*tan(alp2)
                salp2 = _salp0; calp2 = _calp0 * csig2; // No need to normalize

                if (outmask.HasAny(GeodesicFlags.Distance))
                    s12 = arcmode ? _b * (_E0 * sig12 + AB1) : s12_a12;

                if (outmask.HasAny(GeodesicFlags.Longitude))
                {
                    double somg2 = _salp0 * ssig2, comg2 = csig2,  // No need to normalize
                      E = CopySign(1, _salp0);       // east-going?
                                                     // Without normalization we have schi2 = somg2.
                    var cchi2 = _f1 * dn2 * comg2;
                    var chi12 = outmask.HasAny(GeodesicFlags.LongUnroll)
                      ? E * (sig12
                             - (Atan2(ssig2, csig2) - Atan2(_ssig1, _csig1))
                             + (Atan2(E * somg2, cchi2) - Atan2(E * _somg1, _cchi1)))
                      : Atan2(somg2 * _cchi1 - cchi2 * _somg1,
                              cchi2 * _cchi1 + somg2 * _somg1);
                    var lam12 = chi12 -
                      _e2 / _f1 * _salp0 * _H0 *
                      (sig12 + (_E.DeltaH(ssig2, csig2, dn2) - _H1));
                    var lon12 = lam12 / Degree;
                    lon2 = outmask.HasAny(GeodesicFlags.LongUnroll) ? _lon1 + lon12 :
                      AngNormalize(AngNormalize(_lon1) +
                                         AngNormalize(lon12));
                }

                if (outmask.HasAny(GeodesicFlags.Latitude))
                    lat2 = Atan2d(sbet2, _f1 * cbet2);

                if (outmask.HasAny(GeodesicFlags.Azimuth))
                    azi2 = Atan2d(salp2, calp2);

                if (outmask.HasAny(GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale))
                {
                    var J12 = _k2 * _D0 * (sig12 + (_E.DeltaD(ssig2, csig2, dn2) - _D1));
                    if (outmask.HasAny(GeodesicFlags.ReducedLength))
                        // Add parens around (_csig1 * ssig2) and (_ssig1 * csig2) to ensure
                        // accurate cancellation in the case of coincident points.
                        m12 = _b * ((dn2 * (_csig1 * ssig2) - _dn1 * (_ssig1 * csig2))
                                  - _csig1 * csig2 * J12);
                    if (outmask.HasAny(GeodesicFlags.GeodesicScale))
                    {
                        var t = _k2 * (ssig2 - _ssig1) * (ssig2 + _ssig1) / (_dn1 + dn2);
                        M12 = csig12 + (t * ssig2 - csig2 * J12) * _ssig1 / _dn1;
                        M21 = csig12 - (t * _ssig1 - _csig1 * J12) * ssig2 / dn2;
                    }
                }

                if (outmask.HasAny(GeodesicFlags.Area))
                {
                    double B42;
                    if (_A4 == 0)
                        B42 = 0;
                    else
                    {
                        fixed (void* pC4a = __C4a)
                        {
                            var _C4a = new Span<double>(pC4a, _nC4);
                            B42 = DST.Integral(ssig2, csig2, _C4a);
                        }
                    }

                    double salp12, calp12;
                    if (_calp0 == 0 || _salp0 == 0)
                    {
                        // alp12 = alp2 - alp1, used in atan2 so no need to normalize
                        salp12 = salp2 * _calp1 - calp2 * _salp1;
                        calp12 = calp2 * _calp1 + salp2 * _salp1;
                        // We used to include here some patch up code that purported to deal
                        // with nearly meridional geodesics properly.  However, this turned out
                        // to be wrong once _salp1 = -0 was allowed (via
                        // GeodesicExact::InverseLine).  In fact, the calculation of {s,c}alp12
                        // was already correct (following the IEEE rules for handling signed
                        // zeros).  So the patch up code was unnecessary (as well as
                        // dangerous).
                    }
                    else
                    {
                        // tan(alp) = tan(alp0) * sec(sig)
                        // tan(alp2-alp1) = (tan(alp2) -tan(alp1)) / (tan(alp2)*tan(alp1)+1)
                        // = calp0 * salp0 * (csig1-csig2) / (salp0^2 + calp0^2 * csig1*csig2)
                        // If csig12 > 0, write
                        //   csig1 - csig2 = ssig12 * (csig1 * ssig12 / (1 + csig12) + ssig1)
                        // else
                        //   csig1 - csig2 = csig1 * (1 - csig12) + ssig12 * ssig1
                        // No need to normalize
                        salp12 = _calp0 * _salp0 *
                          (csig12 <= 0 ? _csig1 * (1 - csig12) + ssig12 * _ssig1 :
                           ssig12 * (_csig1 * ssig12 / (1 + csig12) + _ssig1));
                        calp12 = Sq(_salp0) + Sq(_calp0) * _csig1 * csig2;
                    }
                    S12 = _c2 * Atan2(salp12, calp12) + _A4 * (B42 - _B41);
                }

                return arcmode ? s12_a12 : sig12 / Degree;
            }

            private void LineInit(GeodesicExact g,
                      double lat1, double lon1,
                      double azi1, double salp1, double calp1,
                      GeodesicFlags caps)
            {
                _lat1 = LatFix(lat1);
                _lon1 = lon1;
                _azi1 = azi1;
                _salp1 = salp1;
                _calp1 = calp1;
                _a = g._a;
                _f = g._f;
                _b = g._b;
                _c2 = g._c2;
                _f1 = g._f1;
                _e2 = g._e2;
                _nC4 = g._nC4;
                // Always allow latitude and azimuth and unrolling of longitude
                caps = caps | GeodesicFlags.Latitude | GeodesicFlags.Azimuth | GeodesicFlags.LongUnroll;

                _caps |=
                    caps &
                    (GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth |
                     GeodesicFlags.Distance | GeodesicFlags.Area | GeodesicFlags.LongUnroll);

                if (caps.HasFlag(GeodesicFlags.DistanceIn))
                {
                    _caps |= GeodesicFlags.DistanceIn.Flags() | (GeodesicFlags)GeodesicCapability.E;
                }
                if (caps.HasFlag(GeodesicFlags.ReducedLength))
                {
                    _caps |= GeodesicFlags.ReducedLength.Flags() | (GeodesicFlags)GeodesicCapability.D;
                }
                if (caps.HasFlag(GeodesicFlags.GeodesicScale))
                {
                    _caps |= GeodesicFlags.GeodesicScale.Flags() | (GeodesicFlags)GeodesicCapability.D;
                }

                SinCosd(AngRound(_lat1), out var sbet1, out var cbet1); sbet1 *= _f1;
                // Ensure cbet1 = +epsilon at poles
                Norm(ref sbet1, ref cbet1); cbet1 = Max(tiny_, cbet1);
                _dn1 = (_f >= 0 ? Sqrt(1 + g._ep2 * Sq(sbet1)) :
                        Sqrt(1 - _e2 * Sq(cbet1)) / _f1);

                // Evaluate alp0 from sin(alp1) * cos(bet1) = sin(alp0),
                _salp0 = _salp1 * cbet1; // alp0 in [0, pi/2 - |bet1|]
                                         // Alt: calp0 = hypot(sbet1, calp1 * cbet1).  The following
                                         // is slightly better (consider the case salp1 = 0).
                _calp0 = Hypot(_calp1, _salp1 * sbet1);
                // Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
                // sig = 0 is nearest northward crossing of equator.
                // With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
                // With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
                // With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
                // Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
                // With alp0 in (0, pi/2], quadrants for sig and omg coincide.
                // No atan2(0,0) ambiguity at poles since cbet1 = +epsilon.
                // With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi.
                _ssig1 = sbet1; _somg1 = _salp0 * sbet1;
                _csig1 = _comg1 = sbet1 != 0 || _calp1 != 0 ? cbet1 * _calp1 : 1;
                // Without normalization we have schi1 = somg1.
                _cchi1 = _f1 * _dn1 * _comg1;
                Norm(ref _ssig1, ref _csig1); // sig1 in (-pi, pi]
                                              // Math::norm(_somg1, _comg1); -- don't need to normalize!
                                              // Math::norm(_schi1, _cchi1); -- don't need to normalize!

                _k2 = Sq(_calp0) * g._ep2;

                _E.Reset(-_k2, -g._ep2, 1 + _k2, 1 + g._ep2);

                if (_caps.Capabilities().HasFlag(GeodesicCapability.E))
                {
                    _E0 = _E.E() / (PI / 2);
                    _E1 = _E.DeltaE(_ssig1, _csig1, _dn1);
                    double s = Sin(_E1), c = Cos(_E1);
                    // tau1 = sig1 + B11
                    _stau1 = _ssig1 * c + _csig1 * s;
                    _ctau1 = _csig1 * c - _ssig1 * s;
                    // Not necessary because Einv inverts E
                    //    _E1 = -_E.deltaEinv(_stau1, _ctau1);
                }

                if (_caps.Capabilities().HasFlag(GeodesicCapability.D))
                {
                    _D0 = _E.D() / (PI / 2);
                    _D1 = _E.DeltaD(_ssig1, _csig1, _dn1);
                }

                if (_caps.Capabilities().HasFlag(GeodesicCapability.H))
                {
                    _H0 = _E.H() / (PI / 2);
                    _H1 = _E.DeltaH(_ssig1, _csig1, _dn1);
                }

                if (_caps.Capabilities().HasFlag(GeodesicCapability.C4))
                {
                    // Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0)
                    _A4 = Sq(_a) * _calp0 * _salp0 * _e2;
                    if (_A4 == 0)
                        _B41 = 0;
                    else
                    {
                        var i4 = new GeodesicExact.I4Integrand(g._ep2, _k2);

                        fixed (void* pC4a = __C4a)
                        {
                            var _C4a = new Span<double>(pC4a, _nC4);
                            g._fft.Transform(ref i4, _C4a);
                            _B41 = DST.Integral(_ssig1, _csig1, _C4a);
                        }
                    }
                }

                _a13 = _s13 = double.NaN;
            }
        }
    }
}
