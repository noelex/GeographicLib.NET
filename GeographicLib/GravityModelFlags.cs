using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib
{
    public enum GravityModelFlags : uint
    {
        None = 0U,

        Gravity = GravityModelCapability.G,

        Disturbance = GravityModelCapability.Delta,

        DisturbingPotential = GravityModelCapability.T,

        SphericalAnomaly = GravityModelCapability.Delta | GravityModelCapability.Gamma,

        GeoidHeight = GravityModelCapability.T | GravityModelCapability.C | GravityModelCapability.Gamma0,

        All = GravityModelCapability.All,
    }
}
