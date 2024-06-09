package sgp4

import (
	"math"
	"time"
)

// Convert earth centered intertial coordinates into lat, long, and altitude.
// Takes ECI output and Golang time.Time input from the SGP4/SDP4 propogator. Returns
// lat and long (in degrees).
func ECItoLLA(eci Vector, t time.Time) LLA {
	// Reference: https://celestrak.org/columns/v02n03/
	var l LLA
	a := 6378.137     // earth semi-major axis
	b := 6356.7523142 // earth semi-minor axis
	f := (a - b) / a  // flatenning
	e2 := ((2 * f) - f*f)

	gmst := Gmst(t)

	sqx2y2 := math.Sqrt(eci.x*eci.x + eci.y*eci.y)

	// spherical earth coordinates
	long := math.Atan2(eci.y, eci.x) - gmst
	lat := math.Atan2(eci.z, sqx2y2)

	// oblate earth fix
	c := 0.0
	for i := 0; i < 20; i++ {
		c = 1 / math.Sqrt(1-e2*(math.Sin(lat)*math.Sin(lat)))
		lat = math.Atan2(eci.z+(a*c*e2*math.Sin(lat)), sqx2y2)
	}

	l.Alt = (sqx2y2 / math.Cos(lat)) - (a * c)
	// v = math.Sqrt(398600.4418 / (l.alt + 6378.137))

	l.Long = math.Mod(long/pi*180, 360)
	if l.Long > 180 {
		l.Long = 360 - l.Long
	} else if l.Long < -180 {
		l.Long = 360 + l.Long
	}

	l.Lat = lat / math.Pi * 180

	return l
}

// Calculates gmst from a Golang time.Time value.
func Gmst(t time.Time) float64 {
	// Based on the 1982 IAU precession model. As described by David Hammen
	// https://astronomy.stackexchange.com/questions/21002/how-to-find-greenwich-mean-sideral-time
	jday := JDay(t)
	_, ut := math.Modf(jday + 0.5)
	jday = jday - ut
	tu := (jday - 2451545.0) / 36525.0
	gmst := 24110.54841 + tu*(8640184.812866+tu*(0.093104-tu*6.2e-6))
	gmst = math.Mod(gmst+86400.0*1.00273790934*ut, 86400.0)
	return 2 * math.Pi * gmst / 86400.0
}

func JDay(t time.Time) float64 {
	j := 367.0*float64(t.Year()) -
		math.Floor((7*(float64(t.Year())+math.Floor((float64(int(t.Month()))+9)/12.0)))*0.25) +
		math.Floor(275*float64(int(t.Month()))/9.0) +
		float64(t.Day()) +
		1721013.5 +
		((float64(t.Second())/60.0+float64(t.Minute()))/60.0+float64(t.Hour()))/24.0

	return j
}

// Calculate look angles for given satellite position in ECI and observer position in lat, long (deg, decimal deg) and alt
//
//	in meters. Returns azmuth in degrees, range in km, and elevation in degrees.
//
// Reference https://celestrak.com/columns/v02n02/
func ECItoLookAngles(eci Vector, obs LLA, t time.Time) (az, rg, el float64) {
	g := Gmst(t)
	re := 6378.137

	obs.Alt /= 1000 // convert observer alt from meters to km

	thetaLat := math.Mod(g+obs.Lat, twopi)
	thetaLong := math.Mod(g+obs.Long, twopi)

	r := (re + obs.Alt) * math.Cos(obs.Lat)
	obsx := r * math.Cos(thetaLong)
	obsy := r * math.Sin(thetaLong)
	obsz := (re + obs.Alt) * math.Sin(obs.Lat)

	rx := eci.x - obsx
	ry := eci.y - obsy
	rz := eci.z - obsz

	tops := math.Sin(obs.Lat)*math.Cos(thetaLat)*rx + math.Sin(obs.Lat)*math.Sin(thetaLat)*ry - math.Cos(obs.Lat)*rz
	tope := -math.Sin(thetaLat)*rx + math.Cos(thetaLat)*ry
	topz := math.Cos(obs.Lat)*math.Cos(thetaLat)*rx + math.Cos(obs.Lat)*math.Sin(thetaLat)*ry + math.Sin(obs.Lat)*rz

	az = math.Atan(-tope / tops)

	if tops > 0 {
		az = az + pi
	}
	if az < 0 {
		az = az + twopi
	}

	az *= radDeg // convert azmuth radians to degrees

	rg = math.Sqrt(rx*rx + ry*ry + rz*rz)
	el = math.Asin(topz / rg)
	el *= radDeg // convert elevation radians to degrees

	return az, rg, el
}
