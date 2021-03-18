using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib
{
    /// <summary>
    /// Encapsulates the results of the rhumb problem solutions.
    /// </summary>
    public abstract class RhumbResult
    {
        /// <summary>
        /// Gets a value representing <i>S12</i>, the area under the rhumb line, in meters^2.
        /// </summary>
        public virtual double Area { get; protected internal set; }
    }

    /// <summary>
    /// Encapsulates the results of the direct rhumb problem solutions.
    /// </summary>
    public class DirectRhumbResult : RhumbResult
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
    /// Encapsulates the results of the inverse rhumb problem solutions.
    /// </summary>
    public class InverseRhumbResult : RhumbResult
    {
        /// <summary>
        /// Gets a value representing <i>azi12</i>, the azimuth of the rhumb line, in degrees.
        /// </summary>
        public virtual double Azimuth { get; protected internal set; }

        /// <summary>
        /// Gets a value representing <i>s12</i>, the distance between point 1 and point 2, in meters.
        /// </summary>
        public virtual double Distance { get; protected internal set; }
    }
}
