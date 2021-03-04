using System;
using System.Collections.Generic;
using System.Text;

namespace GeographicLib.SphericalHarmonics
{
    /// <summary>
    /// Supported normalizations for associated Legendre polynomials.
    /// </summary>
    public enum Normalization
    {
        /// <summary>
        /// Fully normalized associated Legendre polynomials.
        /// </summary>
        /// <remarks>
        /// These are defined by <i>P</i><i>nm</i>^full(<i>z</i>) = (−1)^<i>m</i> sqrt(<i>k</i> (2<i>n</i> + 1)
        /// (<i>n</i> − <i>m</i>)! / (<i>n</i> + <i>m</i>)!) <b>P</b><i>n</i>^<i>m</i>(<i>z</i>), where <b>P</b><i>n</i>^<i>m</i>(<i>z</i>)
        /// is Ferrers function (also known as the Legendre function on the cut or the associated Legendre polynomial)
        /// <a href="https://dlmf.nist.gov/14.7.E10"></a> and <i>k</i> = 1 for <i>m</i> = 0 and <i>k</i> = 2 otherwise.
        /// <para>
        /// The mean squared value of <i>P</i><i>nm</i>^full(cosθ) cos(<i>m</i>λ) and <i>P</i><i>nm</i>^full(cosθ) sin(<i>m</i>λ) over the sphere is 1.
        /// </para>
        /// </remarks>
        Full,

        /// <summary>
        /// Schmidt semi-normalized associated Legendre polynomials.
        /// </summary>
        /// <remarks>
        /// These are defined by <i>P</i><i>nm</i>^schmidt(<i>z</i>) = (−1)<i>m</i> sqrt(<i>k</i> (<i>n</i> − <i>m</i>)! / (<i>n</i> + <i>m</i>)!) 
        /// <b>P</b><i>n</i>^<i>m</i>(<i>z</i>),where <b>P</b><i>n</i>^<i>m</i>(<i>z</i>) 
        /// is Ferrers function (also known as the Legendre function on the cut or the associated Legendre polynomial)
        /// <a href="https://dlmf.nist.gov/14.7.E10"></a> and <i>k</i> = 1 for <i>m</i> = 0 and <i>k</i> = 2 otherwise.
        /// <para>
        /// The mean squared value of <i>P</i><i>nm</i>^schmidt(cosθ) cos(<i>m</i>λ) and <i>P</i><i>nm</i>^schmidt(cosθ) sin(<i>m</i>λ) 
        /// over the sphere is 1/(2<i>n</i> + 1).
        /// </para>
        /// </remarks>
        Schmidt
    }
}
