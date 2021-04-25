#if !NET5_0_OR_GREATER
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib
{
    static class OperatingSystem
    {
        public static bool IsWindows()
        {
            var windir = Environment.GetEnvironmentVariable("windir");
            return !string.IsNullOrEmpty(windir) && windir.Contains(@"\") && Directory.Exists(windir);
        }

        public static bool IsMacOS() => File.Exists(@"/System/Library/CoreServices/SystemVersion.plist");
    }
}
#endif
