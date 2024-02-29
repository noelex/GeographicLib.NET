using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Jobs;

namespace GeographicLib.Benchmarks
{
    [MemoryDiagnoser]
    [SimpleJob(RuntimeMoniker.Net80)]
    [SimpleJob(RuntimeMoniker.Net70)]
    [SimpleJob(RuntimeMoniker.Net60, baseline: true)]
    public class GeodesicBenchmark
    {
        private static readonly IGeodesicLike
            _wgs84 = Geodesic.WGS84,
            _wgs84exact = new Geodesic(Geodesic.WGS84, exact: true),
            _rhumb = Rhumb.WGS84;

        [Params("Geodesic", "GeodesicExact", "Rhumb")]
        public string Target { get; set; }

        private IGeodesicLike Geod =>
            Target == "Geodesic" ?
                _wgs84 :
                    Target == "GeodesicExact" ? _wgs84exact : _rhumb;

        [Benchmark]
        public double Direct()
            => Geod.GenDirect(1, 2, 3, false, 4,
                GeodesicFlags.Standard, out _, out _, out _, out _, out _, out _, out _, out _);

        [Benchmark]
        public double Inverse()
            => Geod.GenInverse(1, 2, 3, 4,
                GeodesicFlags.Standard, out _, out _, out _, out _, out _, out _, out _);
    }
}
