using System;
using System.Collections.Generic;
using System.Text;

namespace GeographicLib
{
    /// <summary>
    /// Exposes equatorial radius (<i>a</i>) and flatterning (<i>f</i>) of an ellipsoid.
    /// </summary>
    public interface IEllipsoid
    {
        /// <summary>
        /// Gets a value representing the equatorial radius (<i>a</i>) of the ellipsoid.
        /// </summary>
        double EquatorialRadius { get; }

        /// <summary>
        /// Gets a value representing the flattening (<i>f</i>) of the ellipsoid.
        /// </summary>
        double Flattening { get; }
    }
}
