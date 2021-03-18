using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib
{
    /// <summary>
    /// Encapsulates the results of the solution of geodesic problems.
    /// </summary>
    public abstract class GeodesicResult
    {
        /// <summary>
        /// Gets a value representing <i>s12</i>, the distance between point 1 and point 2, in meters.
        /// </summary>
        public virtual double Distance { get; protected internal set; }

        /// <summary>
        /// Gets a value representing <i>m12</i>, reduced length of the geodesic, in meters.
        /// </summary>
        public virtual double ReducedLength { get; protected internal set; }

        /// <summary>
        /// Gets a value representing <i>M12</i>, the dimensionless geodesic scale of point 2 relative to point 1.
        /// </summary>
        public virtual double GeodesicScale12 { get; protected internal set; }

        /// <summary>
        /// Gets a value representing <i>M21</i>, the dimensionless geodesic scale of point 1 relative to point 2.
        /// </summary>
        public virtual double GeodesicScale21 { get; protected internal set; }

        /// <summary>
        /// Gets a value representing <i>S12</i>, the area under the geodesic in meters^2.
        /// </summary>
        public virtual double Area { get; protected internal set; }

        /// <summary>
        /// Gets a value representing <i>a12</i>, the arc length of between point 1 and point 2, in degrees.
        /// </summary>
        public virtual double ArcLength { get; protected internal set; }

        /// <summary>
        /// Gets a value representing <i>azi2</i>, the forward azimuth at point 2, in degrees.
        /// </summary>
        public virtual double Azimuth2 { get; protected internal set; }
    }

    /// <summary>
    /// Encapsulates the results of the solution of direct geodesic problems.
    /// </summary>
    public class DirectGeodesicResult : GeodesicResult
    {
        /// <summary>
        /// Gets a value representing <i>lat2</i>, latitude of point 2, in degrees.
        /// </summary>
        public virtual double Latitude { get; protected internal set; }

        /// <summary>
        /// Gets a value representing <i>lon2</i>, longitude of point 2, in degrees.
        /// </summary>
        public virtual double Longitude { get; protected internal set; }
    }

    /// <summary>
    /// Encapsulates the results of the solution of inverse geodesic problems.
    /// </summary>
    public class InverseGeodesicResult : GeodesicResult
    {
        /// <summary>
        /// Gets a value representing <i>azi1</i>, the forward azimuth at point 1, in degrees.
        /// </summary>
        public virtual double Azimuth1 { get; protected internal set; }
    }
}
