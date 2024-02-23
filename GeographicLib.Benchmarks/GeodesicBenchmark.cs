using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Jobs;

namespace GeographicLib.Benchmarks
{
    [SimpleJob(RuntimeMoniker.Net80)]
    [SimpleJob(RuntimeMoniker.Net70)]
    [SimpleJob(RuntimeMoniker.Net60, baseline: true)]
    [SimpleJob(RuntimeMoniker.Net50)]
    [SimpleJob(RuntimeMoniker.NetCoreApp31)]
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
