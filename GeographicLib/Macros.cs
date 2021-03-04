using System;
using System.Collections.Generic;
using System.Text;

namespace GeographicLib
{
    internal static class Macros
    {
        /// <summary>
        /// The precision of floating point numbers used in %GeographicLib.  1 means
        /// float (single precision); 2 (the default) means double; 3 means long double;
        /// 4 is reserved for quadruple precision.  Nearly all the testing has been
        /// carried out with doubles and that's the recommended configuration.  In order
        /// for long double to be used, GEOGRAPHICLIB_HAVE_LONG_DOUBLE needs to be
        /// defined.  Note that with Microsoft Visual Studio, long double is the same as
        /// double.
        /// </summary>
        public const int GEOGRAPHICLIB_PRECISION = 2;

        public static bool GEOGRAPHICLIB_PANIC =>
            // Signal a convergence failure with multiprec types by throwing an exception
            // at loop exit.
            // throw new GeographicException("Convergence failure")

            // Ignore convergence failures with standard floating points types by allowing
            // loop to exit cleanly.
            false;

        public const int GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER = GEOGRAPHICLIB_PRECISION == 2 ? 6 : (GEOGRAPHICLIB_PRECISION == 1 ? 4 : 8);

        /// <summary>
        /// The order of the expansions used by Geodesic. GEOGRAPHICLIB_GEODESIC_ORDER can be set to any integer in [3, 8].
        /// </summary>
        public const int GEOGRAPHICLIB_GEODESIC_ORDER = 
            (GEOGRAPHICLIB_PRECISION == 2 ? 6 : (GEOGRAPHICLIB_PRECISION == 1 ? 3 : (GEOGRAPHICLIB_PRECISION == 3 ? 7 : 8)));

        public const int GEOGRAPHICLIB_GEODESICEXACT_ORDER = 30;
    }
}
