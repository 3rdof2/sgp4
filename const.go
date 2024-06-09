package sgp4

import "math"

// GLOBAL VARIABLES

const (
	x2o3   float64 = 2.0 / 3.0
	pi     float64 = math.Pi
	twopi  float64 = 2.0 * pi
	degRad float64 = pi / 180.0
	radDeg float64 = 180 / pi
	temp4  float64 = 1.5e-12
	xpdotp float64 = 1440.0 / (2.0 * pi)
)

var ss float64 = 78.0/g.radius + 1.0

// GRAVITATIONAL CONSTANTS WGS84

type Grav struct{ mu, radius, xke, tumin, j2, j3, j4, j3oj2 float64 }

var g Grav = Grav{
	// WGS84 gravitational constants
	mu:     398600.5,
	radius: 6378.137,
	xke:    0.07436685316871385,
	tumin:  13.446851082,
	j2:     0.00108262998905,
	j3:     -0.00000253215306,
	j4:     -0.00000161098761,
	j3oj2:  -0.00233889055874,
}
