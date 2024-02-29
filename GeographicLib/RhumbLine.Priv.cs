using static GeographicLib.MathEx;
using static System.Math;

namespace GeographicLib
{
    partial class RhumbLine
    {
        internal struct Priv
        {
            public readonly Rhumb _rh;
            public readonly double _lat1, _lon1, _azi12, _salp, _calp, _mu1, _psi1;
            public readonly AuxAngle _phi1, _chi1;

            public Priv(Rhumb rh, double lat1, double lon1, double azi12)
            {
                (_rh, _lat1, _lon1, _azi12) = (rh, LatFix(lat1), lon1, AngNormalize(azi12));

                SinCosd(_azi12, out _salp, out _calp);
                _phi1 = AuxAngle.FromDegrees(lat1);
                _mu1 = _rh._aux.Convert(AuxLatitudeType.Phi, AuxLatitudeType.Mu,
                                        _phi1, _rh.IsExact).Degrees;
                _chi1 = _rh._aux.Convert(AuxLatitudeType.Phi, AuxLatitudeType.Chi,
                                         _phi1, _rh.IsExact);
                _psi1 = _chi1.Lam;
            }

            /// <summary>
            /// The general position routine. <see cref="Position(double, out double, out double)"/> is defined in term so this function.
            /// </summary>
            /// <param name="s12">distance between point 1 and point 2 (meters); it can be negative.</param>
            /// <param name="outmask">a bitor'ed combination of <see cref="GeodesicFlags"/> values specifying which of the following parameters should be set.</param>
            /// <param name="lat2">latitude of point 2 (degrees).</param>
            /// <param name="lon2">longitude of point 2 (degrees).</param>
            /// <param name="S12">area under the rhumb line (meters^2).</param>
            /// <remarks>
            /// The <see cref="GeodesicFlags"/> values possible for <paramref name="outmask"/> are
            /// <list type="bullet">
            /// <item>outmask |= <see cref="GeodesicFlags.Latitude"/> for the latitude <paramref name="lat2"/>;</item>
            /// <item>outmask |= <see cref="GeodesicFlags.Longitude"/> for the longitude <paramref name="lon2"/>;</item>
            /// <item>outmask |= <see cref="GeodesicFlags.Area"/> for the area <paramref name="S12"/>;</item>
            /// <item>outmask |= <see cref="GeodesicFlags.All"/> for all of the above;</item>
            /// <item>outmask |= <see cref="GeodesicFlags.LongUnroll"/> to unroll <paramref name="lon2"/> instead of wrapping it into the range [−180°, 180°].</item>
            /// </list>
            /// With the <see cref="GeodesicFlags.LongUnroll"/> bit set, the quantity <paramref name="lon2"/> − <i>lon1</i> indicates
            /// how many times and in what sense the rhumb line encircles the ellipsoid.
            /// <para>
            /// If s12 is large enough that the rhumb line crosses a pole,
            /// the longitude of point 2 is indeterminate (a <see cref="double.NaN"/> is returned for <paramref name="lon2"/> and <paramref name="S12"/>).
            /// </para>
            /// </remarks>
            public void GenPosition(double s12, GeodesicFlags outmask, out double lat2, out double lon2, out double S12)
            {
                lat2 = lon2 = S12 = double.NaN;

                double
                  r12 = s12 / (_rh._rm * Degree), // scaled distance in degrees
                  mu12 = r12 * _calp,
                  mu2 = _mu1 + mu12;
                double lat2x, lon2x;
                if (Abs(mu2) <= QD)
                {
                    AuxAngle mu2a = AuxAngle.FromDegrees(mu2),
                             phi2 = _rh._aux.Convert(AuxLatitudeType.Mu, AuxLatitudeType.Phi,
                                                  mu2a, _rh.IsExact),
                             chi2 = _rh._aux.Convert(AuxLatitudeType.Phi, AuxLatitudeType.Chi,
                                                  phi2, _rh.IsExact);
                    lat2x = phi2.Degrees;
                    double dmudpsi = _rh.IsExact ?
                      _rh._aux.DRectifying(_phi1, phi2) / _rh._aux.DIsometric(_phi1, phi2) :
                      _rh._aux.DConvert(AuxLatitudeType.Chi, AuxLatitudeType.Mu, _chi1, chi2)
                      / DAuxLatitude.Dlam(_chi1.Tan, chi2.Tan);
                    lon2x = r12 * _salp / dmudpsi;
                    if (outmask.HasFlag(GeodesicFlags.Area))
                        S12 = _rh._c2 * lon2x * _rh.MeanSinXi(_chi1, chi2);
                    lon2x = outmask.HasFlag(GeodesicFlags.LongUnroll) ? _lon1 + lon2x :
                      AngNormalize(AngNormalize(_lon1) + lon2x);
                }
                else
                {
                    // Reduce to the interval [-180, 180)
                    mu2 = AngNormalize(mu2);
                    // Deal with points on the anti-meridian
                    if (Abs(mu2) > QD) mu2 = AngNormalize(HD - mu2);
                    lat2x = _rh._aux.Convert(AuxLatitudeType.Mu, AuxLatitudeType.Phi,
                                             AuxAngle.FromDegrees(mu2), _rh.IsExact).Degrees;
                    lon2x = double.NaN;
                    if (outmask.HasFlag(GeodesicFlags.Area))
                        S12 = double.NaN;
                }
                if (outmask.HasFlag(GeodesicFlags.Latitude)) lat2 = lat2x;
                if (outmask.HasFlag(GeodesicFlags.Longitude)) lon2 = lon2x;
            }
        }
    }
}
