using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib
{
    /// <summary>
    /// 
    /// </summary>
    public static class DMS
    {
        private const string hemispheres_ = "SNWE";
        private const string signs_ = "-+";
        private const string digits_ = "0123456789";
        private const string dmsindicators_ = "D'\":";

        private static readonly string[] components_ = new[] { "degrees", "minutes", "seconds" };
    }
}
