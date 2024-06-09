package sgp4

import (
	"fmt"
	"log"
	"math"
	"strconv"
	"time"
)

func OMMinit(o *OMM) *satrec {
	var s satrec

	s.satnum = o.NORAD_CAT_ID
	s.epoch, _ = time.Parse("2006-01-02T15:04:05.999999", o.EPOCH)
	s.ndot = o.MEAN_MOTION_DOT / (xpdotp * 1440.0)
	s.nddot = o.MEAN_MOTION_DDOT / (xpdotp * 1440.0 * 1440.0)
	s.bstar = o.BSTAR
	s.inclo = o.INCLINATION * degRad
	s.nodeo = o.RA_OF_ASC_NODE * degRad
	s.ecco = o.ECCENTRICITY
	s.argpo = o.ARG_OF_PERICENTER * degRad
	s.mo = o.MEAN_ANOMALY * degRad
	s.no = o.MEAN_MOTION / (xpdotp)

	sgp4init(&s)

	return &s
}

func TLEinit(tle1, tle2 string) *satrec {
	var s satrec

	s.satnum = tle1[2:7]
	s.inclo = parseFt(tle2[9:16])
	s.inclo *= degRad
	s.nodeo = parseFt(tle2[18:25])
	s.ecco, _ = strconv.ParseFloat("0."+tle2[26:33], 64)
	s.argpo = parseFt(tle2[34:42])
	s.argpo *= degRad
	s.mo = parseFt(tle2[43:51])
	s.mo *= degRad
	s.no = parseFt(tle2[52:63])
	s.no /= xpdotp
	s.nddot = parseM(tle1[45:52])
	s.nddot /= (xpdotp * 1400.0 * 1440.0)
	s.bstar = parseM(tle1[54:61])
	s.ndot, _ = strconv.ParseFloat(tle1[34:43], 64)
	if tle1[33:34] == "-" {
		s.ndot *= -1
	}
	s.ndot /= (xpdotp * 1440.0)

	str, _ := strconv.ParseFloat(tle1[23:32], 64)
	hr := math.Floor(24 * str)
	min := math.Floor((24*str - hr) * 60)
	sec := ((24*str-hr)*60 - min) * 60
	x := fmt.Sprintf("20%v-%vT%v:%v:%v", tle1[18:20], tle1[20:23], hr, min, math.Floor(sec*1000000)/1000000)
	s.epoch, _ = time.Parse("2006-002T15:04:05.999999", x)

	sgp4init(&s)

	return &s
}

func parseFt(x string) float64 {
	y := x
	for i := 0; i < len(x); i++ {
		if x[i:i+1] == " " {
			y = x[i+1:]
		} else {
			break
		}
	}
	z, err := strconv.ParseFloat(y, 64)
	if err != nil {
		log.Fatal("Error parsing float, ", err)
	}
	return z
}

func parseM(x string) float64 {
	a := len(x)
	n, err := strconv.Atoi(x[a-1:])
	if err != nil {
		log.Fatal("Error parsing 00000-0, ", err)
	}
	x = "0." + x[:a-2]
	y, err := strconv.ParseFloat(x, 64)
	if err != nil {
		log.Fatal("Error parsing 00000-0, ", err)
	}
	return y * math.Pow10(-n)
}
