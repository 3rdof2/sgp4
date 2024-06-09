package sgp4

import "math"

func initl(s *satrec) (ainv, ao, con42, cosio, cosio2, eccsq, omeosq, posq, rp, rteosq, sinio float64) {
	// calculate auxillary epoch quantities
	eccsq = s.ecco * s.ecco
	omeosq = 1.0 - eccsq
	rteosq = math.Sqrt(omeosq)
	cosio = math.Cos(s.inclo)
	sinio = math.Sin(s.inclo)
	cosio2 = cosio * cosio

	// adjust mean motion to un-kozai
	ak := math.Pow(g.xke/s.no, x2o3)
	d1 := 0.75 * g.j2 * (3.0*cosio2 - 1.0) / (rteosq * omeosq)
	del := d1 / (ak * ak)
	adel := ak * (1.0 - del*del - del*(1.0/3.0+134.0*del*del/81.0))
	del = d1 / (adel * adel)
	s.no = s.no / (1.0 + del)

	ao = math.Pow(g.xke/s.no, x2o3)
	po := ao * omeosq
	con42 = 1.0 - 5.0*cosio2
	s.con41 = -con42 - cosio2 - cosio2
	ainv = 1.0 / ao
	posq = po * po
	rp = ao * (1.0 - s.ecco)

	s.method = "n"

	s.epoch_jday = JDay(s.epoch)

	s.epochF = s.epoch_jday - 2433281.5 // fractional jday conversion

	// find gst at the time of observation
	tut1 := (s.epoch_jday - 2451545.0) / 36525.0
	s.gsto = -6.2e-6*tut1*tut1*tut1 + 0.093104*tut1*tut1 + (876600.0*3600+8640184.812866)*tut1 + 67310.54841
	s.gsto = math.Mod((s.gsto * degRad / 240.0), twopi)

	// TODO
	// Test case for gsto < 0.0

	if s.gsto < 0.0 {
		s.gsto += twopi
	}

	return ainv, ao, con42, cosio, cosio2, eccsq, omeosq, posq, rp, rteosq, sinio
}
