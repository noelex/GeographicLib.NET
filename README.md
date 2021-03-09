# GeographicLib.NET

[![master](https://github.com/noelex/GeographicLib.NET/actions/workflows/master.yml/badge.svg)](https://github.com/noelex/GeographicLib.NET/actions/workflows/master.yml)
[![develop](https://github.com/noelex/GeographicLib.NET/actions/workflows/develop.yml/badge.svg)](https://github.com/noelex/GeographicLib.NET/actions/workflows/develop.yml)

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
 - [x] 'Exact' version of Geodesic, GeodesicLine and TransverseMercator projection
 - [x] Coordinate conversions
 - [ ] PolygonArea
 - [ ] Rhumb
 - [ ] Geoid
 - [ ] GravityModel
 - [x] MagneticModel
 - [x] Geocode conversions
 - [x] DMS encoding/decoding
 - [ ] More tests

Geodesic and GeodesicLine have been tested with the [test set for geodesic](https://zenodo.org/record/32156#.YCFzsFBLQ_0).
Output of Geodesic and GeodesicLine are verified against the output of the original GeodSolve utility at binary level.

## Installing
### Stable
The library is still under development, thus no stable package release is available currently.

### Nightly ![Nightly Builds](https://buildstats.info/nuget/GeographicLib.NET?includePreReleases=true)
Nightly builds of GeographicLib.NET are hosted on NuGet pre release channel.
You can install them using the following command:

```
dotnet add package GeographicLib.NET --prerelease
```

## Mathematical Functions
GeographicLib uses several C mathematical functions that are not present in or not available in all versions of .NET. These functions are:
 - remquo
 - hypot
 - log1p
 - expm1
 - scalbn (available since .NET 5.0)
 - copysign (available since .NET 5.0)
 - atanh (available since .NET Standard 2.1)
 - asinh (available since .NET Standard 2.1)
 - cbrt (available since .NET Standard 2.1)

GeographicLib.NET provides managed implemetations of these functions (ported from musl libc).

`GeographicLib.MathEx` class will use implemetations provided by .NET runtime whenenver possible, and will fallback to use managed implemetations when not available in .NET runtime. 

You can also force `GeographicLib.MathEx` to fallback to native implemations provided by system C runtime libraries, rather than managed implementaions.
These functions provide better performance, but may produce completely different results in some edge cases.

GeographicLib.NET fallbacks to managed implemtations by default to ensure compatibility between different platforms.
You can tell GeographicLib.NET to fallback to native implemetations by setting `GeographicLib.MathEx.UseManagedCMath` property to `false`.

## Documentation
GeographicLib.NET includes a detailed XML documentation for all public APIs.
Since the API surface of GeographicLib.NET is highly compatible with the original GeographicLib,
you can also refer the original documentation [here](https://geographiclib.sourceforge.io/html/index.html) for usage and explanation.
