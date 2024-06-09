# Overview

`sgp4` is a Golang translation of the C++ astrodynamics algorithm publsihed by [CelesTrak](https://celestrak.org/software/vallado/cpp.zip), written by David Vallado on 28 Jun 2005, and based on the methodology originally published through the AIAA Spacetrack Report #3.

# Use

`sgp4` takes one of two inputs. First is a struct derived from the Orbit Mean-Elements Message (OMM), the second is by a Two Line Element (TLE).

## Initialize the Satrec Coefficients

These functions take an OMM or TLE input and return a pointer to the sgp4 Satrec coefficients that can be used to propogate the orbital elements.

### OMM

```go
func OMMinit(o *OMM) *satrec
```

Where the `OMM` struct is defined as and map directly to the space-track.org OMM datatype.

```go
type OMM struct {
    NORAD_CAT_ID      string  // "44071"
    EPOCH             string  // "2023-11-30T14:51:58.008096"
    MEAN_MOTION       float64 // 1.00271227
    MEAN_MOTION_DOT   float64 // 0.00000055
    MEAN_MOTION_DDOT  float64 // 0.0
    BSTAR             float64 // 0.0
    INCLINATION       float64 // 0.0155
    RA_OF_ASC_NODE    float64 // 354.7718
    ECCENTRICITY      float64 // 0.00003150
    ARG_OF_PERICENTER float64 // 254.5248
    MEAN_ANOMALY      float64 // 103.2389
}
```

### TLE

```go
func TLEinit(tle1, tle2 string) *satrec
```

## Propogate

`Propogate` takes a pointer to an initialized `satrec` struct and a time at propogation. Note that error increases the further the time diverges from the EPOCH.

```go
func Propogate(s *satrec, tprop time.Time) (eci, v Vector)
```

Where `eci` is the current position in the Earth Centered Inertial frame and the corresponding `v` velocity vectors.

## Conversions

### ECI to LLA

Converts Earth Centered Inertial to Latitude, Longitude, and Altitude at the identified `t` time. 

```go
func ECItoLLA(eci Vector, t time.Time) LLA
```

### ECI to Look Angles

Given a satellite position `eci`, observer position `obs` (in LLA), and time.

```go
func ECItoLookAngles(eci Vector, obs LLA, t time.Time) (as, rg, el float64)
```

# Testing

Unit tests were built by taking individual function inputs and outputs from the initial C++ libary. 