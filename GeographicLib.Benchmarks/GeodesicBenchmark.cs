using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Jobs;
using System;
using System.Collections.Generic;
using System.Text;

namespace GeographicLib.Benchmarks
{
    [SimpleJob(RuntimeMoniker.Net60)]
    [SimpleJob(RuntimeMoniker.Net50)]
    [SimpleJob(RuntimeMoniker.NetCoreApp31, baseline: true)]
    [SimpleJob(RuntimeMoniker.NetCoreApp21)]
    public class GeodesicBenchmark
    {
        [Params(true, false)]
        public bool UseManagedCMath { get => MathEx.UseManagedCMath; set => MathEx.UseManagedCMath = value; }


        [Benchmark]
        public double Direct()
            => Geodesic.WGS84.GenDirect(1, 2, 3, false, 4, 
                GeodesicFlags.Standard, out _, out _, out _, out _, out _, out _, out _, out _);

        [Benchmark]
        public double Inverse()
            => Geodesic.WGS84.GenInverse(1, 2, 3, 4, 
                GeodesicFlags.Standard, out _, out _, out _, out _, out _, out _, out _);
    }
}
