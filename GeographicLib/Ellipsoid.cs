﻿using System;
using System.Collections.Generic;
using System.Text;
using GeographicLib.Projections;

using static System.Math;
using static GeographicLib.MathEx;

namespace GeographicLib
{
    /// <summary>
    /// Defines properties of an ellipsoid.
    /// </summary>
    /// <remarks>
    /// This class returns various properties of the ellipsoid and converts
    /// between various types of latitudes.The latitude conversions are also
    /// possible using the various projections supported by GeographicLib; but
    /// <see cref="Ellipsoid"/> provides more direct access(sometimes using private functions
    /// of the projection classes).  <see cref="RectifyingLatitude(double)"/>,
    /// <see cref="InverseRectifyingLatitude(double)"/>, and <see cref="MeridianDistance(double)"/>
    /// provide functionality which can be provided by the <see cref="Geodesic"/> class.
    /// However <see cref="Geodesic"/> uses a series approximation(valid for abs <i>f</i> &lt; 1/150),
    /// whereas <see cref="Ellipsoid"/> computes these quantities using EllipticFunction which
    /// provides accurate results even when <i>f</i> is large.Use of this class
    /// should be limited to -3 &lt; <i>f</i> &lt; 3/4 (i.e., 1/4 &lt; b/a &lt; 4).
    /// </remarks>
    public class Ellipsoid : IEllipsoid
    {
        private const int _numit = 10;

        private readonly double stol_;
        private readonly double _a, _f, _f1, _f12, _e2, _es, _e12, _n, _b;

        private readonly TransverseMercator _tm;
        private readonly AlbersEqualArea _au;
        private readonly EllipticFunction _ell;

        /// <summary>
        /// Constructor for a ellipsoid with equatorial radius and flattening.
        /// </summary>
        /// <param name="a">Equatorial radius (meters).</param>
        /// <param name="f">Flattening of ellipsoid.  Setting <i>f</i> = 0 gives a sphere.</param>
        public Ellipsoid(double a, double f)
        {
            stol_ = 0.01 * Sqrt(DBL_EPSILON);
            _a = a;
            _f = f;
            _f1 = 1 - _f;
            _f12 = Sq(_f1);
            _e2 = _f * (2 - _f);
            _es = (_f < 0 ? -1 : 1) * Sqrt(Abs(_e2));
            _e12 = _e2 / (1 - _e2);
            _n = _f / (2 - _f);
            _b = _a * _f1;

            _tm = new TransverseMercator(_a, _f, 1d);
            _ell = new EllipticFunction(-_e12);
            _au = new AlbersEqualArea(_a, _f, 0d, 1d, 0d, 1d, 1d);
        }

        #region Properties

        /// <summary>
        /// Gets a value representing the equatorial radius of the ellipsoid (meters).  This is the value used in the constructor.
        /// </summary>
        public double EquatorialRadius => _a;

        /// <summary>
        /// Gets a value representing the polar semi-axis (meters).
        /// </summary>
        public double MinorRadius => _b;

        /// <summary>
        /// Gets a value representing the distance between the equator and a pole along a
        /// meridian(meters).  For a sphere <i>L</i> = (φ/2) <i>a</i>.The radius
        /// of a sphere with the same meridian length is <i>L</i> / (φ/2).
        /// </summary>
        public double QuarterMeridian => _b * _ell.E();

        /// <summary>
        /// Gets a value representing the total area of the ellipsoid (meters^2).  For
        /// a sphere <i>A</i> = 4π <i>a</i>^2.  The radius of a sphere
        /// with the same area is sqrt(<i>A</i> / (4π))
        /// </summary>
        public double Area
            => 4 * PI *
              ((Sq(_a) + Sq(_b) *
                (_e2 == 0 ? 1 :
                 (_e2 > 0 ? Atanh(Sqrt(_e2)) : Atan(Sqrt(-_e2))) /
                 Sqrt(Abs(_e2)))) / 2);

        /// <summary>
        /// Gets a value representing the total volume of the ellipsoid (meters^3).
        /// For a sphere <i>V</i> = (4π / 3) <i>a</i>^3.  The radius of
        /// a sphere with the same volume is cbrt(<i>V</i> / (4π/3)).
        /// </summary>
        public double Volume => (4 * PI) * Sq(_a) * _b / 3;

        /// <summary>
        /// <i>f</i> = (<i>a</i> - <i>b</i>) / <i>a</i>, the flattening of the
        /// ellipsoid.This is the value used in the constructor.This is zero,
        /// positive, or negative for a sphere, oblate ellipsoid, or prolate
        /// ellipsoid.
        /// </summary>
        public double Flattening => _f;

        /// <summary>
        /// <i>f</i> ' = (<i>a</i> - <i>b</i>) / <i>b</i>, the second flattening of
        /// the ellipsoid.This is zero, positive, or negative for a sphere,
        /// oblate ellipsoid, or prolate ellipsoid.
        /// </summary>
        public double SecondFlatterning => _f / (1 - _f);

        /// <summary>
        /// <i>n</i> = (<i>a</i> - <i>b</i>) / (<i>a</i> + <i>b</i>), the third flattening
        /// of the ellipsoid.This is zero, positive, or negative for a sphere,
        /// oblate ellipsoid, or prolate ellipsoid.
        /// </summary>
        public double ThirdFlatterning => _n;

        /// <summary>
        /// <i>e</i>^2 = (<i>a</i>^2 - <i>b</i>^2) / <i>a</i>^2, the eccentricity squared
        /// of the ellipsoid.This is zero, positive, or negative for a sphere,
        /// oblate ellipsoid, or prolate ellipsoid.
        /// </summary>
        public double EccentricitySq => _e2;

        /// <summary>
        /// <i>e'</i> ^2 = (<i>a</i>^2 - <i>b</i>^2) / <i>b</i>^2, the second eccentricity
        /// squared of the ellipsoid.This is zero, positive, or negative for a
        /// sphere, oblate ellipsoid, or prolate ellipsoid.
        /// </summary>
        public double SecondEccentricitySq => _e12;

        /// <summary>
        /// <i>e''</i> ^2 = (<i>a</i>^2 - <i>b</i>^2) / (<i>a</i>^2 + <i>b</i>^2),
        /// the third eccentricity squared of the ellipsoid.This is zero,
        /// positive, or negative for a sphere, oblate ellipsoid, or prolate
        /// ellipsoid.
        /// </summary>
        public double ThirdEccentricitySq => _e2 / (2 - _e2);

        #endregion

        #region Latitude conversion methods

        /// <summary>
        /// Parametric latitude conversion.
        /// </summary>
        /// <param name="phi">the geographic latitude (degrees).</param>
        /// <returns>β, the parametric latitude (degrees).</returns>
        /// <remarks>
        /// <para>The geographic latitude, φ, is the angle between the equatorial
        /// plane and a vector normal to the surface of the ellipsoid.</para>
        /// <para>
        /// The parametric latitude (also called the reduced latitude), β,
        /// allows the cartesian coordinated of a meridian to be expressed
        /// conveniently in parametric form as
        /// </para>
        /// <para>
        /// <i>R</i> = <i>a</i> cos β;
        /// </para>
        /// <para>
        /// <i>Z</i> = <i>b</i> sin β;
        /// </para>
        /// <para>
        /// where <i>a</i> and <i>b</i> are the equatorial radius and the polar semi-axis.
        /// For a sphere β = φ.
        /// </para>
        /// <para>
        /// φ must lie in the range [-90°, 90°]; the
        /// result is undefined if this condition does not hold.The returned value
        /// β lies in [-90°, 90°].
        /// </para>
        /// </remarks>
        public double ParametricLatitude(double phi) => Atand(_f1 * Tand(LatFix(phi)));

        /// <summary>
        /// Inverse of <see cref="ParametricLatitude(double)"/>.
        /// </summary>
        /// <remarks>
        /// β must lie in the range [-90°, 90°]; the
        /// result is undefined if this condition does not hold.The returned value
        /// φ lies in [-90°, 90°].
        /// </remarks>
        /// <param name="beta">the parametric latitude (degrees).</param>
        /// <returns>φ, the geographic latitude (degrees).</returns>
        public double InverseParametricLatitude(double beta) => Atand(Tand(LatFix(beta)) / _f1);

        /// <summary>
        /// Geocentric latitude conversion.
        /// </summary>
        /// <remarks>
        /// <para>
        /// The geocentric latitude, θ, is the angle between the equatorial
        /// plane and a line between the center of the ellipsoid and a point on the
        /// ellipsoid.For a sphere θ = φ.
        /// </para>
        /// <para>
        /// φ must lie in the range [-90°, 90°]; the
        /// result is undefined if this condition does not hold.The returned value
        /// θ lies in [-90°, 90°].
        /// </para>
        /// </remarks>
        /// <param name="phi">the geographic latitude (degrees).</param>
        /// <returns>θ, the geocentric latitude (degrees).</returns>
        public double GeocentricLatitude(double phi) => Atand(_f12 * Tand(LatFix(phi)));

        /// <summary>
        /// Inverse of <see cref="GeocentricLatitude(double)"/>.
        /// </summary>
        /// <param name="theta">the geocentric latitude (degrees).</param>
        /// <returns>φ, the geographic latitude (degrees).</returns>
        /// <remarks>
        /// θ must lie in the range [-90°, 90°]; the
        /// result is undefined if this condition does not hold.The returned value
        /// φ lies in [-90°, 90°].
        /// </remarks>
        public double InverseGeocentricLatitude(double theta) => Atand(Tand(LatFix(theta)) / _f12);

        /// <summary>
        /// Rectifying latitude conversion.
        /// </summary>
        /// <param name="phi">the geographic latitude (degrees).</param>
        /// <returns>µ, the rectifying latitude (degrees).</returns>
        /// <remarks>
        /// <para>
        /// The rectifying latitude, µ, has the property that the distance along
        /// a meridian of the ellipsoid between two points with rectifying latitudes
        /// µ1 and µ2 is equal to
        /// (µ2 - µ1) <i>L</i> / 90°,
        /// where <i>L</i> = <see cref="QuarterMeridian"/>.For a sphere µ = φ.
        /// </para>
        /// <para>
        /// φ must lie in the range [-90°, 90°]; the
        /// result is undefined if this condition does not hold.The returned value
        /// µ lies in [-90°, 90°].
        /// </para>
        /// </remarks>
        public double RectifyingLatitude(double phi) => Abs(phi) == 90 ? phi : 90 * MeridianDistance(phi) / QuarterMeridian;

        /// <summary>
        /// Inverse of <see cref="RectifyingLatitude(double)"/>.
        /// </summary>
        /// <param name="mu">the rectifying latitude (degrees).</param>
        /// <returns>φ, the geographic latitude (degrees).</returns>
        /// <remarks>
        /// µ must lie in the range [-90°, 90°]; the
        /// result is undefined if this condition does not hold.The returned value
        /// φ lies in [-90°, 90°].
        /// </remarks>
        public double InverseRectifyingLatitude(double mu)
            => Abs(mu) == 90 ? mu
                : InverseParametricLatitude(_ell.Einv(mu * _ell.E() / 90) / Degree);

        /// <summary>
        /// Authalic latitude conversion.
        /// </summary>
        /// <param name="phi">the geographic latitude (degrees).</param>
        /// <returns>ξ, the authalic latitude (degrees).</returns>
        /// <remarks>
        /// The authalic latitude, ξ, has the property that the area of the
        /// ellipsoid between two circles with authalic latitudes
        /// ξ1 and ξ2 is equal to (sin ξ2 - sin ξ1) <i>A</i> / 2, where <i>A</i>
        /// = <see cref="Area"/>.  For a sphere ξ = φ.
        /// <para>
        /// φ must lie in the range [-90°, 90°]; the
        /// result is undefined if this condition does not hold.The returned value
        /// ξ lies in [-90°, 90°].
        /// </para>
        /// </remarks>
        public double AuthalicLatitude(double phi) => Atand(_au.Txif(Tand(LatFix(phi))));

        /// <summary>
        /// Inverse of <see cref="AuthalicLatitude(double)"/>.
        /// </summary>
        /// <param name="xi">the authalic latitude (degrees).</param>
        /// <returns>φ, the geographic latitude (degrees).</returns>
        /// <remarks>
        /// ξ must lie in the range [-90°, 90°]; the
        /// result is undefined if this condition does not hold.The returned value
        /// φ lies in [-90°, 90°].
        /// </remarks>
        public double InverseAuthalicLatitude(double xi) => Atand(_au.Tphif(Tand(LatFix(xi))));

        /// <summary>
        /// Conformal latitude conversion.
        /// </summary>
        /// <param name="phi">the geographic latitude (degrees).</param>
        /// <returns>χ, the conformal latitude (degrees).</returns>
        /// <remarks>
        /// The conformal latitude, χ, gives the mapping of the ellipsoid to a
        /// sphere which which is conformal(angles are preserved) and in which the
        /// equator of the ellipsoid maps to the equator of the sphere.For a
        /// sphere χ = φ.
        /// <para>
        /// φ must lie in the range [-90°, 90°]; the
        /// result is undefined if this condition does not hold.The returned value
        /// χ lies in [-90°, 90°].
        /// </para>
        /// </remarks>
        public double ConformalLatitude(double phi) => Atand(Taupf(Tand(LatFix(phi)), _es));

        /// <summary>
        /// Inverse of <see cref="ConformalLatitude(double)"/>.
        /// </summary>
        /// <param name="chi">the conformal latitude (degrees).</param>
        /// <returns>φ, the geographic latitude (degrees).</returns>
        /// <remarks>
        /// χ must lie in the range [-90°, 90°]; the
        /// result is undefined if this condition does not hold.The returned value
        /// φ lies in [-90°, 90°].
        /// </remarks>
        public double InverseConformalLatitude(double chi) => Atand(Tauf(Tand(LatFix(chi)), _es));

        /// <summary>
        /// Isometric latitude conversion.
        /// </summary>
        /// <param name="phi">the geographic latitude (degrees).</param>
        /// <returns>Ψ, the isometric latitude (degrees).</returns>
        /// <remarks>
        /// <para>
        /// The isometric latitude gives the mapping of the ellipsoid to a plane
        /// which which is conformal(angles are preserved) and in which the equator
        /// of the ellipsoid maps to a straight line of constant scale; this mapping
        /// defines the Mercator projection.For a sphere Ψ =
        /// sinh^-1 tan φ.
        /// </para>
        /// <para>
        /// φ must lie in the range [-90°, 90°]; the result is
        /// undefined if this condition does not hold.The value returned for φ
        /// = ±90° is some(positive or negative) large but finite value,
        /// such that <see cref="InverseIsometricLatitude"/> returns the original value of φ.
        /// </para>
        /// </remarks>
        public double IsometricLatitude(double phi) => Asinh(Taupf(Tand(LatFix(phi)), _es)) / Degree;

        /// <summary>
        /// Inverse of <see cref="IsometricLatitude(double)"/>.
        /// </summary>
        /// <param name="psi">the isometric latitude (degrees).</param>
        /// <returns>φ, the geographic latitude (degrees).</returns>
        /// <remarks>
        /// The returned value φ lies in [-90°, 90°].  For a
        /// sphere φ = tan^-1 sinh Ψ.
        /// </remarks>
        public double InverseIsometricLatitude(double psi) => Atand(Tauf(Sinh(psi * Degree), _es));

        #endregion

        #region Other quantities

        /// <summary>
        /// Calculates the radius of a circle of specified latitude.
        /// </summary>
        /// <param name="phi">the geographic latitude (degrees).</param>
        /// <returns>
        /// <i>R</i> = <i>a</i> cos β, the radius of a circle of latitude
        /// φ (meters). <i>R</i> (π/180°) gives meters per degree
        /// longitude measured along a circle of latitude.
        /// </returns>
        /// <remarks>
        /// φ must lie in the range [-90°, 90°]; the result is undefined if this condition does not hold.
        /// </remarks>
        public double CircleRadius(double phi) =>
            Abs(phi) == 90 ? 0 :
              // a * cos(beta)
              _a / Hypot(1d, _f1 * Tand(LatFix(phi)));

        /// <summary>
        /// Calculates the distance of a circle of specified latitude from the equator measured parallel to the ellipsoid axis.
        /// </summary>
        /// <param name="phi">the geographic latitude (degrees).</param>
        /// <returns>
        /// <i>Z</i> = <i>b</i> sin β, the distance of a circle of latitude
        ///   φ from the equator measured parallel to the ellipsoid axis
        /// (meters).
        /// </returns>
        /// <remarks>
        /// φ must lie in the range [-90°, 90°]; the result is undefined if this condition does not hold.
        /// </remarks>
        public double CircleHeight(double phi) =>
             // b * sin(beta)
             _b * (_f1 * Tand(phi)) / Hypot(1d, _f1 * Tand(LatFix(phi)));

        /// <summary>
        /// Calculates the distance along a meridian between the equator and a point of specified latitude.
        /// </summary>
        /// <param name="phi">the geographic latitude (degrees).</param>
        /// <returns>
        /// <i>s</i>, the distance along a meridian
        /// between the equator and a point of latitude φ (meters). <i>s</i> is
        ///  given by <i>s</i> = µ <i>L</i> / 90°, where <i>L</i> = <see cref="QuarterMeridian"/>.
        /// </returns>
        /// <remarks>
        /// φ must lie in the range [-90°, 90°]; the result is undefined if this condition does not hold.
        /// </remarks>
        public double MeridianDistance(double phi) => _b * _ell.Ed(ParametricLatitude(phi));

        /// <summary>
        /// Calculates the meridional radius of curvature of the ellipsoid at specified latitude.
        /// </summary>
        /// <param name="phi">the geographic latitude (degrees).</param>
        /// <returns>
        /// ρ, the meridional radius of curvature of the ellipsoid at
        /// latitude φ (meters); this is the curvature of the meridian. ρ is given by ρ = (180°/π) d<i>s</i> / dφ,
        /// where <i>s</i> = <see cref="MeridianDistance(double)"/>; thus ρ (π/180°)
        /// gives meters per degree latitude measured along a meridian.
        /// </returns>
        /// <remarks>
        /// φ must lie in the range [-90°, 90°]; the result is undefined if this condition does not hold.
        /// </remarks>
        public double MeridionalCurvatureRadius(double phi)
        {
            var v = 1 - _e2 * Sq(Sind(LatFix(phi)));
            return _a * (1 - _e2) / (v * Sqrt(v));
        }

        /// <summary>
        /// Calculate the transverse radius of curvature of the ellipsoid at specified latitude.
        /// </summary>
        /// <param name="phi">the geographic latitude (degrees).</param>
        /// <returns>
        /// ν, the transverse radius of curvature of the ellipsoid at
        /// latitude φ (meters); this is the curvature of a curve on the
        /// ellipsoid which also lies in a plane perpendicular to the ellipsoid
        /// and to the meridian.  ν is related to <i>R</i> = <see cref="CircleRadius(double)"/> by <i>R</i> = ν cos φ.
        /// </returns>
        /// <remarks>
        /// φ must lie in the range [-90°, 90°]; the result is undefined if this condition does not hold.
        /// </remarks>
        public double TransverseCurvatureRadius(double phi)
        {
            var v = 1 - _e2 * Sq(Sind(LatFix(phi)));
            return _a / Sqrt(v);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="phi">the geographic latitude (degrees).</param>
        /// <param name="azi">the angle between the meridian and the normal section (degrees).</param>
        /// <returns>
        /// the radius of curvature of the ellipsoid in the normal
        /// section at latitude φ inclined at an angle <paramref name="azi"/> to the
        /// meridian(meters).
        /// </returns>
        /// <remarks>
        /// φ must lie in the range [-90°, 90°]; the result is undefined if this condition does not hold.
        /// </remarks>
        public double NormalCurvatureRadius(double phi, double azi)
        {
            var v = 1 - _e2 * Sq(Sind(LatFix(phi)));
            SinCosd(azi, out var salp, out var calp);
            return _a / (Sqrt(v) * (Sq(calp) * v / (1 - _e2) + Sq(salp)));
        }

        #endregion

        #region Eccentricity conversions

        /// <summary>
        /// Converts second flatterning <i>f</i> ' to flatterning <i>f</i>.
        /// </summary>
        /// <param name="fp"><i>f</i> ' = (<i>a</i> - <i>b</i>) / <i>b</i>, the second flattening.</param>
        /// <returns><i>f</i> = (<i>a</i> - <i>b</i>) / <i>a</i>, the flattening.</returns>
        /// <remarks>
        /// <i>f</i> ' should lie in (-1, ∞). The returned value <i>f</i> lies in (-∞, 1).
        /// </remarks>
        public static double SecondFlatteningToFlattening(double fp) => fp / (1 + fp);

        /// <summary>
        /// Converts flatterning <i>f</i> to second flatterning <i>f</i> '.
        /// </summary>
        /// <param name="f"><i>f</i> = (<i>a</i> - <i>b</i>) / <i>a</i>, the flattening.</param>
        /// <returns><i>f</i> ' = (<i>a</i> - <i>b</i>) / <i>b</i>, the second flattening.</returns>
        /// <remarks>
        /// <i>f</i> should lie in (-∞, 1). The returned value <i>f</i> ' lies in (-1, ∞).
        /// </remarks>
        public static double FlatteningToSecondFlattening(double f) => f / (1 - f);

        /// <summary>
        /// Converts third flatterning <i>n</i> to flatterning <i>f</i>.
        /// </summary>
        /// <param name="n"><i>n</i> = (<i>a</i> - <i>b</i>) / (<i>a</i> + <i>b</i>), the third flattening.</param>
        /// <returns><i>f</i> = (<i>a</i> - <i>b</i>) / <i>a</i>, the flattening.</returns>
        /// <remarks>
        /// <i>n</i> should lie in (-1, 1). The returned value <i>f</i> lies in (-∞, 1).
        /// </remarks>
        public static double ThirdFlatteningToFlattening(double n) => 2 * n / (1 + n);

        /// <summary>
        /// Converts flatterning <i>f</i> to third flatterning <i>n</i>.
        /// </summary>
        /// <param name="f"><i>f</i> = (<i>a</i> - <i>b</i>) / <i>a</i>, the flattening.</param>
        /// <returns><i>n</i> = (<i>a</i> - <i>b</i>) / (<i>a</i> + <i>b</i>), the third flattening.</returns>
        /// <remarks>
        /// <i>f</i> should lie in (-∞, 1). The returned value <i>n</i> lies in (-1, 1).
        /// </remarks>
        public static double FlatteningToThirdFlattening(double f) => f / (2 - f);

        /// <summary>
        /// Converts eccentricity squared <i>e</i>^2 to flattening <i>f</i>.
        /// </summary>
        /// <param name="e2"><i>e</i>^2 = (<i>a</i>^2 - <i>b</i>^2) / <i>a</i>^2, the eccentricity squared.</param>
        /// <returns><i>f</i> = (<i>a</i> - <i>b</i>) / <i>a</i>, the flattening.</returns>
        /// <remarks>
        /// <i>e</i>^2 should lie in (-∞, 1). The returned value <i>f</i> lies in (-∞, 1).
        /// </remarks>
        public static double EccentricitySqToFlattening(double e2) => e2 / (Sqrt(1 - e2) + 1);

        /// <summary>
        /// Converts flattening <i>f</i> to eccentricity squared <i>e</i>^2.
        /// </summary>
        /// <param name="f"><i>f</i> = (<i>a</i> - <i>b</i>) / <i>a</i>, the flattening.</param>
        /// <returns><i>e</i>^2 = (<i>a</i>^2 - <i>b</i>^2) / <i>a</i>^2, the eccentricity squared.</returns>
        /// <remarks>
        /// <i>f</i> should lie in (-∞, 1). The returned value <i>e</i>^2 lies in (-∞, 1).
        /// </remarks>
        public static double FlatteningToEccentricitySq(double f) => f * (2 - f);

        /// <summary>
        /// Converts second eccentricity squared <i>e'</i>^2 to flattening <i>f</i>.
        /// </summary>
        /// <param name="ep2"><i>e'</i>^2 = (<i>a</i>^2 - <i>b</i>^2) / <i>b</i>^2, the second eccentricity squared.</param>
        /// <returns><i>f</i> = (<i>a</i> - <i>b</i>) / <i>a</i>, the flattening.</returns>
        /// <remarks>
        /// <i>e'</i>^2  should lie in (-1, ∞). The returned value <i>f</i> lies in (-∞, 1).
        /// </remarks>
        public static double SecondEccentricitySqToFlattening(double ep2) => ep2 / (Sqrt(1 + ep2) + 1 + ep2);

        /// <summary>
        /// Converts flattening <i>f</i> to second eccentricity squared <i>e'</i>^2.
        /// </summary>
        /// <param name="f"><i>f</i> = (<i>a</i> - <i>b</i>) / <i>a</i>, the flattening.</param>
        /// <returns><i>e'</i>^2 = (<i>a</i>^2 - <i>b</i>^2) / <i>b</i>^2, the second eccentricity squared.</returns>
        /// <remarks>
        /// <i>f</i> should lie in (-∞, 1). The returned value <i>e'</i>^2 lies in (-1, ∞).
        /// </remarks>
        public static double FlatteningToSecondEccentricitySq(double f) => f * (2 - f) / Sq(1 - f);

        /// <summary>
        /// Converts third eccentricity squared <i>e''</i> ^2 to flatterning <i>f</i>.
        /// </summary>
        /// <param name="epp2"><i>e''</i> ^2 = (<i>a</i>^2 - <i>b</i>^2) / (<i>a</i>^2 + <i>b</i>^2), the third eccentricity squared.</param>
        /// <returns><i>f</i> = (<i>a</i> - <i>b</i>) / <i>a</i>, the flattening.</returns>
        /// <remarks>
        /// <i>e''</i> ^2 should lie in (-1,1). The returned value <i>f</i> lies in (-∞, 1).
        /// </remarks>
        public static double ThirdEccentricitySqToFlattening(double epp2) => 2 * epp2 / (Sqrt((1 - epp2) * (1 + epp2)) + 1 + epp2);

        /// <summary>
        /// Converts flatterning <i>f</i> to third eccentricity squared <i>e''</i> ^2.
        /// </summary>
        /// <param name="f"><i>f</i> = (<i>a</i> - <i>b</i>) / <i>a</i>, the flattening.</param>
        /// <returns><i>e''</i> ^2 = (<i>a</i>^2 - <i>b</i>^2) / (<i>a</i>^2 + <i>b</i>^2), the third eccentricity squared.</returns>
        /// <remarks>
        /// <i>f</i> should lie in (-∞, 1). The returned value <i>e''</i> ^2 lies in (-1,1).
        /// </remarks>
        public static double FlatteningToThirdEccentricitySq(double f) => f * (2 - f) / (1 + Sq(1 - f));

        #endregion

        /// <summary>
        /// A global instantiation of Ellipsoid with the parameters for the WGS84 ellipsoid.
        /// </summary>
        /// <returns></returns>
        public static Ellipsoid WGS84 { get; } = new Ellipsoid(Constants.WGS84_a, Constants.WGS84_f);
    }
}
