using System;
using System.Collections.Generic;
using System.Text;

namespace GeographicLib
{
    /// <summary>
    /// Represents exception thrown by GeographicLib.
    /// </summary>
    public class GeographicException:Exception
    {
        /// <summary>
        /// Initialize a new <see cref="GeographicException"/> with specified error message.
        /// </summary>
        /// <param name="msg">Error message.</param>
        public GeographicException(string msg) : base(msg) { }
    }
}
