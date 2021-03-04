# GeographicLib.NET

![ubuntu](https://github.com/noelex/GeographicLib.NET/workflows/ubuntu/badge.svg)
![windows](https://github.com/noelex/GeographicLib.NET/workflows/windows/badge.svg)
![NuGet Badge](https://buildstats.info/nuget/GeographicLib.NET?includePreReleases=true)]

GeographicLib.NET is a native .NET implementation of [GeographicLib](https://sourceforge.net/p/geographiclib) written in pure C#.

## What's different from NETGeographicLib

Unlike the original NETGeographicLib, GeographicLib.NET is implemented in pure C# without binding the original C++ GeograpbicLib library by using C++/CLI or P/Invoke, thus achieves higher level of portability.

You should be able to use GeographicLib.NET with any target framework and platform that supports .NET Standard 2.0.

## Features
Currently GeographicLib.NET has ported a few features such as Ellipsoid, Geodesic and various projection classes. More features will be ported and the ultimate goal is to implement all features provided by GeographicLib. 

Here's a simple list of what features are implemented and planned.
 - [x] Implement missing math functions
 - [x] Ellipoid and EllipticFunction
 - [x] Geodesic and GeodesicLine
 - [x] Projections
 - [x] Spherical harmonic series
 - [ ] 'Exact' version of Geodesic, GeodesicLine and TransverseMercator projection
 - [ ] Coordinate conversions
 - [ ] PolygonArea
 - [ ] Rhumb
 - [ ] Geoid
 - [ ] GravityModel
 - [x] MagneticModel
 - [x] Geocode conversions
 - [ ] DMS parser
 - [ ] Utility programs
 - [ ] More tests

Geodesic and GeodesicLine have been tested with the [test set for geodesic](https://zenodo.org/record/32156#.YCFzsFBLQ_0).
Output of Geodesic and GeodesicLine are verified against the output of the original GeodSolve utility at binary level.

## Nuget package
This is still a work in progress, thus no Nuget package is available yet. 
Packages will be released once all major features are implemented and the public API surface gets stablized.

For now, you'll have to build the library from source. Functionalities are not fully tested, use at your own risk.

## Documentation
GeographicLib.NET includes a detailed XML documentation for all public APIs.
Since the API surface of GeographicLib.NET is highly compatible with the original GeographicLib,
you can also refer the original documentation [here](https://geographiclib.sourceforge.io/html/index.html) for usage and explanation.
