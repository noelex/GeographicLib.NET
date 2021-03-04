using System;
using System.Collections.Generic;
using System.Text;

namespace GeographicLib
{
    /// <summary>
    /// Base class of various geodesic line implementations.
    /// </summary>
    public abstract class GeodesicLineBase:IGeodesicLine
    {
        #region Properties

        /// <inheritdoc/>
        public abstract double Arc { get; set; }

        /// <inheritdoc/>
        public abstract double Azimuth { get; }

        /// <inheritdoc/>
        public abstract GeodesicFlags Capabilities { get; }

        /// <inheritdoc/>
        public abstract double CosineAzimuth { get; }

        /// <inheritdoc/>
        public abstract double CosineEquatorialAzimuth { get; }

        /// <inheritdoc/>
        public abstract double Distance { get; set; }

        /// <inheritdoc/>
        public abstract double EquatorialArc { get; }

        /// <inheritdoc/>
        public abstract double EquatorialAzimuth { get; }

        /// <inheritdoc/>
        public abstract double Latitude { get; }

        /// <inheritdoc/>
        public abstract double Longitude { get; }

        /// <inheritdoc/>
        public abstract double SineAzimuth { get; }

        /// <inheritdoc/>
        public abstract double SineEquatorialAzimuth { get; }

        /// <inheritdoc/>
        public abstract double EquatorialRadius { get; }

        /// <inheritdoc/>
        public abstract double Flattening { get; }

        #endregion

        /// <inheritdoc/>
        public abstract double GenPosition(bool arcmode, double s12_a12, GeodesicFlags outmask, 
            out double lat2, out double lon2, out double azi2, out double s12, out double m12, out double M12, out double M21, out double S12);

        /// <inheritdoc/>
        public abstract double GetDistance(bool arcmode);

        /// <inheritdoc/>
        public abstract void SetDistance(bool arcmode, double s13_a13);

        /// <inheritdoc/>
        public bool HasCapability(GeodesicFlags testcaps) => Capabilities.Flags().HasFlag(testcaps);

        #region Position in terms of arc length

        /// <inheritdoc/>
        public void ArcPosition(double a12, out double lat2, out double lon2, out double azi2,
                     out double s12, out double m12, out double M12, out double M21,
                     out double S12) =>
            GenPosition(true, a12,
                        GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth | GeodesicFlags.Distance |
                        GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale | GeodesicFlags.Area,
                        out lat2, out lon2, out azi2, out s12, out m12, out M12, out M21, out S12);

        /// <inheritdoc/>
        public void ArcPosition(double a12, out double lat2, out double lon2) =>
            GenPosition(true, a12,
                        GeodesicFlags.Latitude | GeodesicFlags.Longitude,
                        out lat2, out lon2, out _, out _, out _, out _, out _, out _);

        /// <inheritdoc/>
        public void ArcPosition(double a12,
                         out double lat2, out double lon2, out double azi2) =>
            GenPosition(true, a12,
                        GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth,
                        out lat2, out lon2, out azi2, out _, out _, out _, out _, out _);

        /// <inheritdoc/>
        public void ArcPosition(double a12, out double lat2, out double lon2, out double azi2,
                          out double s12) =>
            GenPosition(true, a12,
                        GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth | GeodesicFlags.Distance,
                        out lat2, out lon2, out azi2, out s12, out _, out _, out _, out _);

        /// <inheritdoc/>
        public void ArcPosition(double a12, out double lat2, out double lon2, out double azi2,
                          out double s12, out double m12) =>
            GenPosition(true, a12,
                        GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth |
                        GeodesicFlags.Distance | GeodesicFlags.ReducedLength,
                        out lat2, out lon2, out azi2, out s12, out m12, out _, out _, out _);

        /// <inheritdoc/>
        public void ArcPosition(double a12, out double lat2, out double lon2, out double azi2,
                         out double s12, out double M12, out double M21) =>
            GenPosition(true, a12,
                        GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth |
                        GeodesicFlags.Distance | GeodesicFlags.GeodesicScale,
                        out lat2, out lon2, out azi2, out s12, out _, out M12, out M21, out _);

        /// <inheritdoc/>
        public void ArcPosition(double a12, out double lat2, out double lon2, out double azi2,
                         out double s12, out double m12, out double M12, out double M21) =>
            GenPosition(true, a12,
                        GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth |
                        GeodesicFlags.Distance | GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale,
                        out lat2, out lon2, out azi2, out s12, out m12, out M12, out M21, out _);

        #endregion

        #region Position in terms of distance

        /// <inheritdoc/>
        public double Position(double s12,
                        out double lat2, out double lon2, out double azi2,
                        out double m12, out double M12, out double M21,
                        out double S12) =>
            GenPosition(false, s12,
                         GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth |
                         GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale | GeodesicFlags.Area,
                         out lat2, out lon2, out azi2, out _, out m12, out M12, out M21, out S12);

        /// <inheritdoc/>
        public double Position(double s12, out double lat2, out double lon2) =>
            GenPosition(false, s12,
                        GeodesicFlags.Latitude | GeodesicFlags.Longitude,
                        out lat2, out lon2, out _, out _, out _, out _, out _, out _);

        /// <inheritdoc/>
        public double Position(double s12, out double lat2, out double lon2,
                            out double azi2) =>
            GenPosition(false, s12,
                        GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth,
                        out lat2, out lon2, out azi2, out _, out _, out _, out _, out _);

        /// <inheritdoc/>      /// <returns><i>a12</i>, arc length from point 1 to point 2 (degrees).</returns>
        public double Position(double s12, out double lat2, out double lon2,
                            out double azi2, out double m12) =>
            GenPosition(false, s12,
                        GeodesicFlags.Latitude | GeodesicFlags.Longitude |
                        GeodesicFlags.Azimuth | GeodesicFlags.ReducedLength,
                        out lat2, out lon2, out azi2, out _, out m12, out _, out _, out _);

        /// <inheritdoc/>
        public double Position(double s12, out double lat2, out double lon2,
                            out double azi2, out double M12, out double M21) =>
            GenPosition(false, s12,
                        GeodesicFlags.Latitude | GeodesicFlags.Longitude |
                        GeodesicFlags.Azimuth | GeodesicFlags.GeodesicScale,
                        out lat2, out lon2, out azi2, out _, out _, out M12, out M21, out _);

        /// <inheritdoc/>
        public double Position(double s12,
                            out double lat2, out double lon2, out double azi2,
                            out double m12, out double M12, out double M21) =>
            GenPosition(false, s12,
                        GeodesicFlags.Latitude | GeodesicFlags.Longitude | GeodesicFlags.Azimuth |
                        GeodesicFlags.ReducedLength | GeodesicFlags.GeodesicScale,
                        out lat2, out lon2, out azi2, out _, out m12, out M12, out M21, out _);

        #endregion
    }
}
