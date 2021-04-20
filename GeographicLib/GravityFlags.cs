using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib
{
    /// <summary>
    /// Bit masks for the capabilities to be given to the <see cref="GravityCircle"/> object produced by <see cref="GravityModel.Circle(double, double, GravityFlags)"/>.
    /// </summary>
    [Flags]
    public enum GravityFlags : uint
    {
        /// <summary>
        /// No capabilities.
        /// </summary>
        None = 0U,

        /// <summary>
        /// Allow calls to <see cref="GravityCircle.Gravity(double)"/>,
        /// <see cref="GravityCircle.W(double)"/>, and
        /// <see cref="GravityCircle.V(double)"/>.
        /// </summary>
        Gravity = GravityCapability.G,

        /// <summary>
        /// Allow calls to <see cref="GravityCircle.Disturbance(double)"/> and
        /// <see cref="GravityCircle.Td(double)"/>.
        /// </summary>
        Disturbance = GravityCapability.Delta,

        /// <summary>
        /// Allow calls to <see cref="GravityCircle.T(double)"/>
        /// (i.e., computing the disturbing potential and not the gravity disturbance vector).
        /// </summary>
        DisturbingPotential = GravityCapability.T,

        /// <summary>
        /// Allow calls to <see cref="GravityCircle.SphericalAnomaly(double)"/>.
        /// </summary>
        SphericalAnomaly = GravityCapability.Delta | GravityCapability.Gamma,

        /// <summary>
        /// Allow calls to <see cref="GravityCircle.GeoidHeight(double)"/>.
        /// </summary>
        GeoidHeight = GravityCapability.T | GravityCapability.C | GravityCapability.Gamma0,

        /// <summary>
        /// All capabilities.
        /// </summary>
        All = GravityCapability.All,
    }
}
