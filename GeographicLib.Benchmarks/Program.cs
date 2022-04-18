using BenchmarkDotNet.Running;

namespace GeographicLib.Benchmarks
{
    class Program
    {
        public static void Main(string[] args)
        {
            BenchmarkRunner.Run<GeodesicBenchmark>();
        }
    }
}