using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib
{
    [Flags]
    public enum GravityFlags : uint
    {
        None = 0U,

        Gravity = GravityCapability.G,

        Disturbance = GravityCapability.Delta,

        DisturbingPotential = GravityCapability.T,

        SphericalAnomaly = GravityCapability.Delta | GravityCapability.Gamma,

        GeoidHeight = GravityCapability.T | GravityCapability.C | GravityCapability.Gamma0,

        All = GravityCapability.All,
    }
}
