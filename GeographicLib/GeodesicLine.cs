using System;
using System.Collections.Generic;
using System.Text;

using static System.Math;
using static GeographicLib.MathEx;

namespace GeographicLib
{
    /// <summary>
    /// Represents a geodesic line.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <see cref="GeodesicLine"/> facilitates the determination of a series of points on a single geodesic.
    /// The starting point (<i>lat1</i>, <i>lon1</i>) and the azimuth azi1 are specified in the constructor;
    /// alternatively, the <see cref="Geodesic.Line"/> method can be used to create a <see cref="GeodesicLine"/>.
    /// <see cref="IGeodesicLine.Position(double, out double, out double, out double, out double, out double, out double, out double)"/> returns the location of point 2 a distance s12 along the geodesic.
    /// In addition, <see cref="IGeodesicLine.ArcPosition(double, out double, out double, out double, out double, out double, out double, out double, out double)"/> gives the position of point 2 an arc length <i>a12</i> along the geodesic.
    /// </para>
    /// <para>
    /// You can register the position of a reference point 3 a distance (arc length), <i>s13</i> (<i>a13</i>)
    /// along the geodesic with the <see cref="Distance"/> (<see cref="Arc"/>) functions.
    /// Points a fractional distance along the line can be found by providing, for example, 0.5 * <see cref="Distance"/>
    /// as an argument to GeodesicLine.Position. The <see cref="Geodesic.InverseLine"/> or <see cref="GeodesicBase.DirectLine"/> methods return
    /// <see cref="GeodesicLine"/> objects with point 3 set to the point 2 of the corresponding geodesic problem.
    /// <see cref="GeodesicLine"/> objects created with the public constructor or with <see cref="Geodesic.Line"/>
    /// have <i>s13</i> and <i>a13</i> set to <see cref="double.NaN"/>.
    /// </para>
    /// <para>
    /// The default copy constructor and assignment operators work with this class. 
    /// Similarly, a vector can be used to hold <see cref="GeodesicLine"/> objects.
    /// </para>
    /// <para>
    /// The calculations are accurate to better than 15 nm (15 nanometers).
    /// See Sec. 9 of <a href="https://arxiv.org/abs/1102.1215v1">arXiv:1102.1215v1</a> for details.
    /// The algorithms used by this class are based on series expansions using the flattening f as a small parameter.
    /// These are only accurate for |<i>f</i>| &lt; 0.02; however reasonably accurate results will be obtained for
    /// |<i>f</i>| &lt; 0.2. For very eccentric ellipsoids, use <see cref="GeodesicLineExact"/> instead.
    /// </para>
    /// <para>
    /// The algorithms are described in
    /// <list type="bullet">
    /// <item>
    /// C. F. F. Karney, <a href="https://doi.org/10.1007/s00190-012-0578-z">Algorithms for geodesics</a>,
    /// J. Geodesy 87, 43–55 (2013); DOI: <a href="https://doi.org/10.1007/s00190-012-0578-z">10.1007/s00190-012-0578-z</a>; 
    /// addenda: <a href="https://geographiclib.sourceforge.io/geod-addenda.html">geod-addenda.html</a>.
    /// </item>
    /// </list>
    /// </para>
    /// <para>
    /// For more information on geodesics see <a href="https://geographiclib.sourceforge.io/html/geodesic.html">Geodesics on an ellipsoid of revolution</a>.
    /// </para>
    /// </remarks>
    public class GeodesicLine : GeodesicLineBase
    {
        private const int nC1_ = Geodesic.nC1_;
        private const int nC1p_ = Geodesic.nC1p_;
        private const int nC2_ = Geodesic.nC2_;
        private const int nC3_ = Geodesic.nC3_;
        private const int nC4_ = Geodesic.nC4_;

        private readonly double tiny_;
        private readonly double _lat1, _lon1, _azi1;
        private readonly double _a, _f, _b, _c2, _f1, _salp0, _calp0, _k2,
          _salp1, _calp1, _ssig1, _csig1, _dn1, _stau1, _ctau1, _somg1, _comg1,
          _A1m1, _A2m1, _A3c, _B11, _B21, _B31, _A4, _B41;

        private double _a13, _s13;

        // index zero elements of _C1a, _C1pa, _C2a, _C3a are unused
        private readonly Memory<double>
            _C1a = new double[nC1_ + 1],
            _C1pa = new double[nC1p_ + 1],
            _C2a = new double[nC2_ + 1],
            _C3a = new double[nC3_],
            _C4a = new double[nC4_];    // all the elements of _C4a are used

        private readonly GeodesicFlags _caps;

        internal GeodesicLine(Geodesic g,
                 double lat1, double lon1,
                 double azi1, double salp1, double calp1,
                 GeodesicFlags caps, bool arcmode, double s13_a13)
        {
            LineInit(g, lat1, lon1, azi1, salp1, calp1, caps,
                ref tiny_, ref _lat1, ref _lon1, ref _azi1, ref _salp1, ref _calp1, ref _a, ref _f, ref _b, ref _c2, ref _f1,
                ref _dn1, ref _salp0, ref _calp0, ref _ssig1, ref _somg1, ref _csig1, ref _comg1, ref _k2, ref _A1m1, ref _B11,
                ref _stau1, ref _ctau1, ref _A2m1, ref _B21, ref _A3c, ref _B31, ref _A4, ref _B41, ref _caps);
            SetDistance(arcmode, s13_a13);
        }

        /// <summary>
        /// Constructor for a geodesic line staring at latitude <i>lat1</i>, longitude <i>lon1</i>, 
        /// and azimuth <i>azi1</i> (all in degrees).
        /// </summary>
        /// <param name="g">
        /// A <see cref="Geodesic"/> object used to compute the necessary information about the <see cref="GeodesicLine"/>.
        /// </param>
        /// <param name="lat1">latitude of point 1 (degrees).</param>
        /// <param name="lon1">longitude of point 1 (degrees).</param>
        /// <param name="azi1">azimuth at point 1 (degrees).</param>
        /// <param name="caps">bitor'ed combination of <see cref="GeodesicFlags"/> values specifying the capabilities
        /// the <see cref="GeodesicLine"/> object should possess, i.e., which quantities can be returned in calls to
        /// <see cref="IGeodesicLine.Position(double, out double, out double, out double, out double, out double, out double, out double)"/>.</param>
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
        public GeodesicLine(Geodesic g, double lat1, double lon1, double azi1,
                 GeodesicFlags caps = GeodesicFlags.All)
        {
            azi1 = AngNormalize(azi1);

            // Guard against underflow in salp0.  Also -0 is converted to +0.
            SinCosd(AngRound(azi1), out var salp1, out var calp1);
            LineInit(g, lat1, lon1, azi1, salp1, calp1, caps,
                ref tiny_, ref _lat1, ref _lon1, ref _azi1, ref _salp1, ref _calp1, ref _a, ref _f, ref _b, ref _c2, ref _f1,
                ref _dn1, ref _salp0, ref _calp0, ref _ssig1, ref _somg1, ref _csig1, ref _comg1, ref _k2, ref _A1m1, ref _B11,
                ref _stau1, ref _ctau1, ref _A2m1, ref _B21, ref _A3c, ref _B31, ref _A4, ref _B41, ref _caps);
        }

        private void LineInit(Geodesic g,
                  double lat1, double lon1,
                  double azi1, double salp1, double calp1,
                  GeodesicFlags caps,
                  ref double tiny_, ref double _lat1, ref double _lon1, ref double _azi1, ref double _salp1, ref double _calp1,
                  ref double _a, ref double _f, ref double _b, ref double _c2, ref double _f1, ref double _dn1, ref double _salp0,
                  ref double _calp0, ref double _ssig1, ref double _somg1, ref double _csig1, ref double _comg1, ref double _k2,
                  ref double _A1m1, ref double _B11, ref double _stau1, ref double _ctau1, ref double _A2m1, ref double _B21,
                  ref double _A3c, ref double _B31, ref double _A4, ref double _B41, ref GeodesicFlags _caps)
        {
            tiny_ = g.tiny_;
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

            // Always allow latitude and azimuth and unrolling of longitude
            _caps = caps | GeodesicFlags.Latitude | GeodesicFlags.Azimuth | GeodesicFlags.LongUnroll;

            SinCosd(AngRound(_lat1), out var sbet1, out var cbet1);
            sbet1 *= _f1;

            // Ensure cbet1 = +epsilon at poles
            Norm(ref sbet1, ref cbet1);
            cbet1 = Max(tiny_, cbet1);
            _dn1 = Sqrt(1 + g._ep2 * Sq(sbet1));

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
            Norm(ref _ssig1, ref _csig1); // sig1 in (-pi, pi]
                                          // Math::norm(_somg1, _comg1); -- don't need to normalize!

            _k2 = Sq(_calp0) * g._ep2;
            var eps = _k2 / (2 * (1 + Sqrt(1 + _k2)) + _k2);

            if (_caps.Capabilities().HasFlag(GeodesicCapability.C1))
            {
                _A1m1 = Geodesic.A1m1f(eps);
                Geodesic.C1f(eps, _C1a.Span);
                _B11 = Geodesic.SinCosSeries(true, _ssig1, _csig1, _C1a.Span, nC1_);
                double s = Sin(_B11), c = Cos(_B11);
                // tau1 = sig1 + B11
                _stau1 = _ssig1 * c + _csig1 * s;
                _ctau1 = _csig1 * c - _ssig1 * s;
                // Not necessary because C1pa reverts C1a
                //    _B11 = -SinCosSeries(true, _stau1, _ctau1, _C1pa, nC1p_);
            }

            if (_caps.Capabilities().HasFlag(GeodesicCapability.C1p))
                Geodesic.C1pf(eps, _C1pa.Span);

            if (_caps.Capabilities().HasFlag(GeodesicCapability.C2))
            {
                _A2m1 = Geodesic.A2m1f(eps);
                Geodesic.C2f(eps, _C2a.Span);
                _B21 = Geodesic.SinCosSeries(true, _ssig1, _csig1, _C2a.Span, nC2_);
            }

            if (_caps.Capabilities().HasFlag(GeodesicCapability.C3))
            {
                g.C3f(eps, _C3a.Span);
                _A3c = -_f * _salp0 * g.A3f(eps);
                _B31 = Geodesic.SinCosSeries(true, _ssig1, _csig1, _C3a.Span, nC3_ - 1);
            }

            if (_caps.Capabilities().HasFlag(GeodesicCapability.C4))
            {
                g.C4f(eps, _C4a.Span);
                // Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0)
                _A4 = Sq(_a) * _calp0 * _salp0 * g._e2;
                _B41 = Geodesic.SinCosSeries(false, _ssig1, _csig1, _C4a.Span, nC4_);
            }

            _a13 = _s13 = double.NaN;
        }

        /// <inheritdoc/>
        public override double GenPosition(bool arcmode, double s12_a12, GeodesicFlags outmask,
                           out double lat2, out double lon2, out double azi2,
                           out double s12, out double m12, out double M12, out double M21,
                           out double S12)
        {
            lat2 = lon2 = azi2 = s12 = m12 = M12 = M21 = S12 = double.NaN;
            outmask &= _caps & (GeodesicFlags)GeodesicCapability.OutMask;
            if (!(arcmode || _caps.HasAny(GeodesicFlags.DistanceIn)))
                // Uninitialized or impossible distance calculation requested
                return double.NaN;

            // Avoid warning about uninitialized B12.
            double sig12, ssig12, csig12, B12 = 0, AB1 = 0;
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
                  tau12 = s12_a12 / (_b * (1 + _A1m1)),
                  s = Sin(tau12),
                  c = Cos(tau12);
                // tau2 = tau1 + tau12
                B12 = -Geodesic.SinCosSeries(true,
                                               _stau1 * c + _ctau1 * s,
                                               _ctau1 * c - _stau1 * s,
                                               _C1pa.Span, nC1p_);
                sig12 = tau12 - (B12 - _B11);
                ssig12 = Sin(sig12); csig12 = Cos(sig12);
                if (Abs(_f) > 0.01)
                {
                    // Reverted distance series is inaccurate for |f| > 1/100, so correct
                    // sig12 with 1 Newton iteration.  The following table shows the
                    // approximate maximum error for a = WGS_a() and various f relative to
                    // GeodesicExact.
                    //     erri = the error in the inverse solution (nm)
                    //     errd = the error in the direct solution (series only) (nm)
                    //     errda = the error in the direct solution
                    //             (series + 1 Newton) (nm)
                    //
                    //       f     erri  errd errda
                    //     -1/5    12e6 1.2e9  69e6
                    //     -1/10  123e3  12e6 765e3
                    //     -1/20   1110 108e3  7155
                    //     -1/50  18.63 200.9 27.12
                    //     -1/100 18.63 23.78 23.37
                    //     -1/150 18.63 21.05 20.26
                    //      1/150 22.35 24.73 25.83
                    //      1/100 22.35 25.03 25.31
                    //      1/50  29.80 231.9 30.44
                    //      1/20   5376 146e3  10e3
                    //      1/10  829e3  22e6 1.5e6
                    //      1/5   157e6 3.8e9 280e6
                    double
                      ssig2_ = _ssig1 * csig12 + _csig1 * ssig12,
                      csig2_ = _csig1 * csig12 - _ssig1 * ssig12;
                    B12 = Geodesic.SinCosSeries(true, ssig2_, csig2_, _C1a.Span, nC1_);
                    double serr = (1 + _A1m1) * (sig12 + (B12 - _B11)) - s12_a12 / _b;
                    sig12 = sig12 - serr / Sqrt(1 + _k2 * Sq(ssig2_));
                    ssig12 = Sin(sig12); csig12 = Cos(sig12);
                    // Update B12 below
                }
            }

            double ssig2, csig2, sbet2, cbet2, salp2, calp2;
            // sig2 = sig1 + sig12
            ssig2 = _ssig1 * csig12 + _csig1 * ssig12;
            csig2 = _csig1 * csig12 - _ssig1 * ssig12;
            var dn2 = Sqrt(1 + _k2 * Sq(ssig2));
            if (outmask.HasAny(GeodesicFlags.Distance | GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale))
            {
                if (arcmode || Abs(_f) > 0.01)
                    B12 = Geodesic.SinCosSeries(true, ssig2, csig2, _C1a.Span, nC1_);
                AB1 = (1 + _A1m1) * (B12 - _B11);
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
                s12 = arcmode ? _b * ((1 + _A1m1) * sig12 + AB1) : s12_a12;

            if (outmask.HasAny(GeodesicFlags.Longitude))
            {
                // tan(omg2) = sin(alp0) * tan(sig2)
                double somg2 = _salp0 * ssig2, comg2 = csig2,  // No need to normalize
                  E = CopySign(1, _salp0);       // east-going?
                                                 // omg12 = omg2 - omg1
                var omg12 = outmask.HasAny(GeodesicFlags.LongUnroll)
                  ? E * (sig12
                         - (Atan2(ssig2, csig2) - Atan2(_ssig1, _csig1))
                         + (Atan2(E * somg2, comg2) - Atan2(E * _somg1, _comg1)))
                  : Atan2(somg2 * _comg1 - comg2 * _somg1,
                          comg2 * _comg1 + somg2 * _somg1);
                var lam12 = omg12 + _A3c *
                  (sig12 + (Geodesic.SinCosSeries(true, ssig2, csig2, _C3a.Span, nC3_ - 1)
                             - _B31));
                var lon12 = lam12 / Degree;
                lon2 = outmask.HasAny(GeodesicFlags.LongUnroll) ? _lon1 + lon12 :
                  AngNormalize(AngNormalize(_lon1) + AngNormalize(lon12));
            }

            if (outmask.HasAny(GeodesicFlags.Latitude))
                lat2 = Atan2d(sbet2, _f1 * cbet2);

            if (outmask.HasAny(GeodesicFlags.Azimuth))
                azi2 = Atan2d(salp2, calp2);

            if (outmask.HasAny(GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale))
            {
                double
                  B22 = Geodesic.SinCosSeries(true, ssig2, csig2, _C2a.Span, nC2_),
                  AB2 = (1 + _A2m1) * (B22 - _B21),
                  J12 = (_A1m1 - _A2m1) * sig12 + (AB1 - AB2);
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
                var
                  B42 = Geodesic.SinCosSeries(false, ssig2, csig2, _C4a.Span, nC4_);
                double salp12, calp12;
                if (_calp0 == 0 || _salp0 == 0)
                {
                    // alp12 = alp2 - alp1, used in atan2 so no need to normalize
                    salp12 = salp2 * _calp1 - calp2 * _salp1;
                    calp12 = calp2 * _calp1 + salp2 * _salp1;
                    // We used to include here some patch up code that purported to deal
                    // with nearly meridional geodesics properly.  However, this turned out
                    // to be wrong once _salp1 = -0 was allowed (via
                    // Geodesic::InverseLine).  In fact, the calculation of {s,c}alp12
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

        /// <inheritdoc/>
        public override void SetDistance(bool arcmode, double s13_a13)
        {
            if (arcmode) Arc = s13_a13; else Distance = s13_a13;
        }

        /// <inheritdoc/>
        public override double GetDistance(bool arcmode) => (arcmode ? _a13 : _s13);

        #region Public properties

        /// <inheritdoc/>
        public override double Latitude => _lat1;

        /// <inheritdoc/>
        public override double Longitude => _lon1;

        /// <inheritdoc/>
        public override double Azimuth => _azi1;

        /// <inheritdoc/>
        public override double SineAzimuth => _salp1;

        /// <inheritdoc/>
        public override double CosineAzimuth => _calp1;

        /// <inheritdoc/>
        public override double EquatorialAzimuth => Atan2d(_salp0, _calp0);

        /// <inheritdoc/>
        public override double SineEquatorialAzimuth => _salp0;

        /// <inheritdoc/>
        public override double CosineEquatorialAzimuth => _calp0;

        /// <inheritdoc/>
        public override double EquatorialArc => Atan2d(_ssig1, _csig1);

        /// <inheritdoc/>
        public override double EquatorialRadius => _a;

        /// <inheritdoc/>
        public override double Flattening => _f;

        /// <inheritdoc/>
        public override GeodesicFlags Capabilities => _caps;

        /// <inheritdoc/>
        public override double Distance
        {
            get => GetDistance(false);
            set
            {
                _s13 = value;

                // This will set _a13 to NaN if the GeodesicLine doesn't have the
                // DISTANCE_IN capability.
                _a13 = GenPosition(false, _s13, 0u, out _, out _, out _, out _, out _, out _, out _, out _);
            }
        }

        /// <inheritdoc/>
        public override double Arc
        {
            get => GetDistance(true);
            set
            {
                _a13 = value;
                // In case the GeodesicLine doesn't have the DISTANCE capability.
                _s13 = double.NaN;

                GenPosition(true, _a13, GeodesicFlags.Distance, out _, out _, out _, out _s13, out _, out _, out _, out _);
            }
        }

        #endregion
    }
}
