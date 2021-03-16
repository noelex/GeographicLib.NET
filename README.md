# GeographicLib.NET

[![master](https://github.com/noelex/GeographicLib.NET/actions/workflows/master.yml/badge.svg)](https://github.com/noelex/GeographicLib.NET/actions/workflows/master.yml)
[![develop](https://github.com/noelex/GeographicLib.NET/actions/workflows/develop.yml/badge.svg)](https://github.com/noelex/GeographicLib.NET/actions/workflows/develop.yml)

[GeographicLib](https://sourceforge.net/p/geographiclib) is a small set of C++ classes for performing conversions between geographic,
UTM, UPS, MGRS, geocentric, and local cartesian coordinates,for gravity (e.g., EGM2008), geoid height and geomagnetic field (e.g., WMM2020) calculations,
and for solving geodesic problems.

GeographicLib.NET is a native .NET implementation of [GeographicLib](https://sourceforge.net/p/geographiclib) written in pure C#.

## What's different from NETGeographicLib

Unlike the original NETGeographicLib, GeographicLib.NET is implemented in pure C# without binding the original C++ GeograpbicLib library by using C++/CLI or P/Invoke, thus achieves higher level of portability.

You should be able to use GeographicLib.NET with any target framework and platform that supports .NET Standard 2.0 or above.

## Features
Bellow is a list of planned and implemented features.
 - [x] Projections (`AlbersEqualArea`, `AzimuthalEquidistant`, `CassiniSoldner`, `Gnomonic`, `LambertConformalConic`, `PolarStereographic` and `TransverseMercator`)
 - [x] Geocodes (`GARS`, `Geohash`, `Georef`, `MGRS` and `OSGB`)
 - [x] Coordinate conversions (`UTMUPS`, `Geocentric` and `LocalCartesian`)
 - [x] Coordinate parsing/formatting (`DMS` and `GeoCoords`)
 - [x] Geodesic (`Geodesic`, `GeodesicLine`, `GeodesicExact` and `GeodesicLineExact`)
 - [x] Rhumb (`Rhumb` and `RhumbLine`)
 - [x] PolygonArea (`PolygonArea<T>`, `PolygonArea`, `PolygonAreaExact` and `PolygonAreaRhumb`)
 - [x] Geoid (`Geoid`)
 - [x] GravityModel (`GravityCircle`, `NormalGravity` and `GravityModel`)
 - [x] MagneticModel (`MagneticModel`, `MagneticCircle`)
 - [x] Auxilary classes (`MathEx`, `Ellipoid`, `EllipticFunction` and `SphericalHarmonic`)

`Geodesic` and `GeodesicExact` are tested with the [test set for geodesic](https://zenodo.org/record/32156#.YCFzsFBLQ_0).

`TransverseMercator` and `TransverseMercatorExact` are tested with data generated by 64-bit `TransverseMercatorProj` utility ran on Windows.

Managed implemetation of C mathematical functions in `MathEx` are tested with data generated by 64-bit Windows Universal C Runtime.

## Installing
### Stable ![Release](https://buildstats.info/nuget/GeographicLib.NET?includePreReleases=false)
Stable releases of GeographicLib.NET are hosted on NuGet.
You can install them using the following command:
```
dotnet add package GeographicLib.NET
```

### Preview ![Preview](https://buildstats.info/nuget/GeographicLib.NET?includePreReleases=true)
Preview versions of GeographicLib.NET are hosted on NuGet pre-release channel.
You can install them using the following command:

```
dotnet add package GeographicLib.NET --prerelease
```

## Mathematical Functions
GeographicLib uses several C mathematical functions that are not available in all versions of .NET. These functions are:
 - remquo
 - hypot
 - log1p
 - expm1
 - frexp
 - log2 (available since .NET 5.0)
 - fma (available since .NET 5.0)
 - scalbn (available since .NET 5.0)
 - copysign (available since .NET 5.0)
 - atanh (available since .NET Standard 2.1)
 - asinh (available since .NET Standard 2.1)
 - cbrt (available since .NET Standard 2.1)

GeographicLib.NET provides managed implementations of these functions (ported from [musl libc](https://musl.libc.org/)).

`GeographicLib.MathEx` class will use implementations provided by .NET runtime whenenver possible, and will fallback to managed implementations when not available in .NET runtime. 

You can also force `GeographicLib.MathEx` to fallback to platform dependent implementations provided by system C runtime libraries,
rather than managed implementations, by setting `GeographicLib.MathEx.UseManagedCMath` property to `false`.
These functions provide better performance, but may produce different results on different platform in some edge cases.

## Documentation
GeographicLib.NET includes a detailed XML documentation for all public APIs.
Since the API surface of GeographicLib.NET is highly compatible with the original GeographicLib,
you can also refer the original documentation [here](https://geographiclib.sourceforge.io/html/index.html) for usage and explanation.

## Change Log
GeographicLib.NET adopts changes made in original GeographicLib and aligns its version number with original GeographicLib releases.

Bellow is a list of stable releases of GeographicLib.NET and changes made in .NET side in each release.
For changes adopted from original GeographicLib, please refer the its change log [here](https://geographiclib.sourceforge.io/html/changes.html).

### Version 1.52.0 (unreleased)
- **BREAKING**
  - Modify overloads of `Forward` and `Reverse` in `AlbersEqualArea`, `AzimuthalEquidistant`, `CassiniSoldner` and `LambertConformalConic` to return coordinates as tuples.

- **NEW**
  - Add constructor overloads that accept `IEllipsoid` as parameter for `AlbersEqualArea` and `LambertConformalConic`.
  - Add managed implementation of `log2`.

### Version 1.51.0 (released 2021/03/14)
Initial stable release.