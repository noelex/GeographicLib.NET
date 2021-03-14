using System;
using System.Collections.Generic;
using System.Text;

namespace GeographicLib
{
    /// <summary>
    /// Base class of various geodesic implementations.
    /// </summary>
    public abstract class GeodesicBase : IGeodesic
    {
        /// <inheritdoc/>
        public abstract double EquatorialRadius { get; }

        /// <inheritdoc/>
        public abstract double Flattening { get; }

        /// <inheritdoc/>
        public abstract double EllipsoidArea { get; }

        /// <inheritdoc/>
        public abstract double GenDirect(double lat1, double lon1, double azi1, bool arcmode, double s12_a12,
            GeodesicFlags outmask, out double lat2, out double lon2, out double azi2, out double s12, out double m12,
            out double M12, out double M21, out double S12);

        /// <inheritdoc/>
        public abstract IGeodesicLine GenDirectLine(double lat1, double lon1, double azi1, bool arcmode,
            double s12_a12, GeodesicFlags caps = GeodesicFlags.All);

        /// <inheritdoc/>
        public abstract double GenInverse(double lat1, double lon1, double lat2, double lon2, GeodesicFlags outmask, out double s12, out double azi1, 
            out double azi2, out double m12, out double M12, out double M21, out double S12);

        /// <inheritdoc/>
        public abstract IGeodesicLine InverseLine(double lat1, double lon1, double lat2, double lon2, GeodesicFlags caps = GeodesicFlags.All);

        /// <inheritdoc/>
        public abstract IGeodesicLine Line(double lat1, double lon1, double azi1, GeodesicFlags caps = GeodesicFlags.All);

        /// <inheritdoc/>
        public IGeodesicLine DirectLine(double lat1, double lon1, double azi1, double s12,
            GeodesicFlags caps = GeodesicFlags.All) => GenDirectLine(lat1, lon1, azi1, false, s12, caps);

        /// <inheritdoc/>
        public IGeodesicLine ArcDirectLine(double lat1, double lon1, double azi1, double a12,
                           GeodesicFlags caps = GeodesicFlags.All) => GenDirectLine(lat1, lon1, azi1, true, a12, caps);

        #region Direct geodesic problem specified in terms of arc length

        /// <inheritdoc/>
        public void ArcDirect(double lat1, double lon1, double azi1, double a12,
                out double lat2, out double lon2, out double azi2, out double s12,
                out double m12, out double M12, out double M21, out double S12)
            => GenDirect(lat1, lon1, azi1, true, a12,
                     GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth | GeodesicFlags.Distance |
                     GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale | GeodesicFlags.Area,
                     out lat2, out lon2, out azi2, out s12, out m12, out M12, out M21, out S12);

        /// <inheritdoc/>
        public void ArcDirect(double lat1, double lon1, double azi1, double a12,
                   out double lat2, out double lon2, out double azi2)
            => GenDirect(lat1, lon1, azi1, true, a12,
                      GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth,
                      out lat2, out lon2, out azi2, out _, out _, out _, out _, out _);

        /// <inheritdoc/>
        public void ArcDirect(double lat1, double lon1, double azi1, double a12,
                   out double lat2, out double lon2, out double azi2, out double s12)
            => GenDirect(lat1, lon1, azi1, true, a12,
                  GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth | GeodesicFlags.Distance,
                  out lat2, out lon2, out azi2, out s12, out _, out _, out _, out _);

        /// <inheritdoc/>
        public void ArcDirect(double lat1, double lon1, double azi1, double a12,
                       out double lat2, out double lon2, out double azi2,
                       out double s12, out double m12)
            => GenDirect(lat1, lon1, azi1, true, a12,
                      GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth |
                      GeodesicFlags.Distance | GeodesicFlags.ReducedLength,
                      out lat2, out lon2, out azi2, out s12, out m12, out _, out _, out _);

        /// <inheritdoc/>
        public void ArcDirect(double lat1, double lon1, double azi1, double a12,
                       out double lat2, out double lon2, out double azi2, out double s12,
                       out double M12, out double M21)
            => GenDirect(lat1, lon1, azi1, true, a12,
                      GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth |
                      GeodesicFlags.Distance | GeodesicFlags.GeodesicScale,
                      out lat2, out lon2, out azi2, out s12, out _, out M12, out M21, out _);

        /// <inheritdoc/>
        public void ArcDirect(double lat1, double lon1, double azi1, double a12,
                       out double lat2, out double lon2, out double azi2, out double s12,
                       out double m12, out double M12, out double M21)
            => GenDirect(lat1, lon1, azi1, true, a12,
                      GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth |
                      GeodesicFlags.Distance | GeodesicFlags.GeodesicScale | GeodesicFlags.ReducedLength,
                      out lat2, out lon2, out azi2, out s12, out m12, out M12, out M21, out _);

        #endregion

        #region Direct geodesic problem specified in terms of distance

        /// <inheritdoc/>
        public double Direct(double lat1, double lon1, double azi1, double s12, out double lat2, out double lon2)
            => GenDirect(lat1, lon1, azi1, false, s12,
                GeodesicFlags.Latitude | GeodesicFlags.Longitude,
                out lat2, out lon2, out _, out _, out _, out _, out _, out _);

        /// <inheritdoc/>
        public double Direct(double lat1, double lon1, double azi1, double s12, out double lat2, out double lon2, out double azi2)
            => GenDirect(lat1, lon1, azi1, false, s12,
                GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth,
                out lat2, out lon2, out azi2, out _, out _, out _, out _, out _);

        /// <inheritdoc/>
        public double Direct(double lat1, double lon1, double azi1, double s12, out double lat2, out double lon2, out double azi2, out double m12)
            => GenDirect(lat1, lon1, azi1, false, s12,
                GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth | GeodesicFlags.ReducedLength,
                out lat2, out lon2, out azi2, out _, out m12, out _, out _, out _);

        /// <inheritdoc/>
        public double Direct(double lat1, double lon1, double azi1, double s12, out double lat2, out double lon2, out double azi2, out double M12, out double M21)
            => GenDirect(lat1, lon1, azi1, false, s12,
                GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth | GeodesicFlags.GeodesicScale,
                out lat2, out lon2, out azi2, out _, out _, out M12, out M21, out _);

        /// <inheritdoc/>
        public double Direct(double lat1, double lon1, double azi1, double s12,
            out double lat2, out double lon2, out double azi2, out double m12, out double M12, out double M21)
            => GenDirect(lat1, lon1, azi1, false, s12,
                       GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth |
                       GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale,
                       out lat2, out lon2, out azi2, out _, out m12, out M12, out M21, out _);

        /// <inheritdoc/>
        public double Direct(double lat1, double lon1, double azi1, double s12,
            out double lat2, out double lon2, out double azi2, out double m12, out double M12, out double M21, out double S12)
             => GenDirect(lat1, lon1, azi1, false, s12,
                       GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth |
                       GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale | GeodesicFlags.Area,
                       out lat2, out lon2, out azi2, out _, out m12, out M12, out M21, out S12);
        #endregion

        #region Inverse geodesic problem

        /// <inheritdoc/>
        public double Inverse(double lat1, double lon1, double lat2, double lon2,
                       out double s12, out double azi1, out double azi2, out double m12,
                       out double M12, out double M21, out double S12)
            => GenInverse(lat1, lon1, lat2, lon2,
                       GeodesicFlags.Distance | GeodesicFlags.Azimuth |
                       GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale | GeodesicFlags.Area,
                       out s12, out azi1, out azi2, out m12, out M12, out M21, out S12);

        /// <inheritdoc/>
        public double Inverse(double lat1, double lon1, double lat2, double lon2,
                           out double s12)
            => GenInverse(lat1, lon1, lat2, lon2,
                       GeodesicFlags.Distance,
                       out s12, out _, out _, out _, out _, out _, out _);

        /// <inheritdoc/>
        public double Inverse(double lat1, double lon1, double lat2, double lon2,
                       out double azi1, out double azi2)
            => GenInverse(lat1, lon1, lat2, lon2,
                       GeodesicFlags.Azimuth,
                       out _, out azi1, out azi2, out _, out _, out _, out _);

        /// <inheritdoc/>
        public double Inverse(double lat1, double lon1, double lat2, double lon2,
                           out double s12, out double azi1, out double azi2)
            => GenInverse(lat1, lon1, lat2, lon2,
                       GeodesicFlags.Distance | GeodesicFlags.Azimuth,
                       out s12, out azi1, out azi2, out _, out _, out _, out _);

        /// <inheritdoc/>
        public double Inverse(double lat1, double lon1, double lat2, double lon2,
                           out double s12, out double azi1, out double azi2, out double m12)
            => GenInverse(lat1, lon1, lat2, lon2,
                       GeodesicFlags.Distance | GeodesicFlags.Azimuth | GeodesicFlags.ReducedLength,
                       out s12, out azi1, out azi2, out m12, out _, out _, out _);

        /// <inheritdoc/>
        public double Inverse(double lat1, double lon1, double lat2, double lon2,
                           out double s12, out double azi1, out double azi2,
                           out double M12, out double M21)
            => GenInverse(lat1, lon1, lat2, lon2,
                       GeodesicFlags.Distance | GeodesicFlags.Azimuth | GeodesicFlags.GeodesicScale,
                       out s12, out azi1, out azi2, out _, out M12, out M21, out _);

        /// <inheritdoc/>
        public double Inverse(double lat1, double lon1, double lat2, double lon2,
                           out double s12, out double azi1, out double azi2, out double m12,
                           out double M12, out double M21)
            => GenInverse(lat1, lon1, lat2, lon2,
                       GeodesicFlags.Distance | GeodesicFlags.Azimuth |
                       GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale,
                       out s12, out azi1, out azi2, out m12, out M12, out M21, out _);

        #endregion
    }
}
