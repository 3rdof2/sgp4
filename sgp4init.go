package sgp4

import "math"

func sgp4init(s *satrec) {
	var c cof

	s.isimp = 0
	s.method = "n"
	s.operationmode = "i"
	s.init = "y"

	qzms2t := (120.0 - 78.0) / g.radius // could be renamed as qzms24
	qzms2t = qzms2t * qzms2t * qzms2t * qzms2t

	// SGP4 INITIALIZATION
	_, ao, con42, cosio, cosio2, eccsq, omeosq, posq, rp, rteosq, sinio := initl(s)
	c.eccsq = eccsq

	if rp < 220.0/g.radius+1.0 {
		s.isimp = 1
	} // use case doesn't meet these conditions; deep space only

	sfour := ss
	qzms24 := qzms2t

	// TODO
	// Update satrec to include perige from the OMM measurement or include the equation below if taken from TLE
	// Update test functions to leverage the TLE variant to avoid accuracy conflicts if OMM perige different than TLE version

	perige := (rp - 1.0) * g.radius

	if perige < 156.0 {
		// you would have to be obsurdly close to earth for perige to be less than 156 km
		sfour = perige - 78.0
		if perige < 98.0 {
			sfour = 20.0
		}
		qzms24 = (120.0 - sfour) / g.radius
		qzms24 = qzms24 * qzms24 * qzms24 * qzms24
		sfour = sfour/g.radius + 1.0
	}

	pinvsq := 1.0 / posq
	tsi := 1.0 / (ao - sfour)
	s.eta = ao * s.ecco * tsi

	etasq := s.eta * s.eta
	eeta := s.ecco * s.eta
	psisq := math.Abs(1.0 - etasq)
	coef := qzms24 * tsi * tsi * tsi * tsi
	coef1 := coef / math.Pow(psisq, 3.5)
	cc2 := coef1 * s.no * (ao*(1.0+1.5*etasq+eeta*(4.0+etasq)) + 0.375*g.j2*tsi/psisq*s.con41*
		(8.0+3.0*etasq*(8.0+etasq)))
	s.cc1 = s.bstar * cc2

	cc3 := 0.0

	// TODO
	// Implement test case with eccentricity less thaqn 1.0e-4

	if s.ecco > 1.0e-4 {
		cc3 = -2.0 * coef * tsi * g.j3oj2 * s.no * sinio / s.ecco
	}
	s.x1mth2 = 1.0 - cosio2

	s.cc4 = 2.0 * s.no * coef1 * ao * omeosq * (s.eta*(2.0+0.5*etasq) + s.ecco*
		(0.5+2.0*etasq) - g.j2*tsi/(ao*psisq)*(-3.0*s.con41*(1.0-2.0*eeta+etasq*
		(1.5-0.5*eeta))+0.75*s.x1mth2*(2.0*etasq-eeta*(1.0+etasq))*math.Cos(2.0*s.argpo)))

	s.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75*(etasq+eeta) + eeta*etasq)

	cosio4 := cosio2 * cosio2

	temp1 := 1.5 * g.j2 * pinvsq * s.no
	temp2 := 0.5 * temp1 * g.j2 * pinvsq
	temp3 := -0.46875 * g.j4 * pinvsq * pinvsq * s.no
	s.mdot = s.no + 0.5*temp1*rteosq*s.con41 + 0.0625*temp2*rteosq*(13.0-78.0*cosio2+137.0*cosio4)
	s.argpdot = -0.5*temp1*con42 + 0.0625*temp2*(7.0-114.0*cosio2+395.0*cosio4) + temp3*(3.0-36.0*cosio2+49.0*cosio4)

	xhdot1 := -temp1 * cosio
	s.nodedot = xhdot1 + (0.5*temp2*(4.0-19.0*cosio2)+2.0*temp3*(3.0-7.0*cosio2))*cosio
	c.xpidot = s.argpdot + s.nodedot

	s.omgcof = s.bstar * cc3 * math.Cos(s.argpo)

	s.xmcof = 0.0
	if s.ecco > 1.0e-4 {
		s.xmcof = -x2o3 * coef * s.bstar / eeta
	}

	s.nodecf = 3.5 * omeosq * xhdot1 * s.cc1

	s.t2cof = 1.5 * s.cc1

	s.xlcof = 0.0
	if math.Abs(cosio+1.0) > 1.5e-12 {
		s.xlcof = -0.25 * g.j3oj2 * sinio * (3.0 + 5.0*cosio) / (1.0 + cosio)
	} else {
		s.xlcof = -0.25 * g.j3oj2 * sinio * (3.0 + 5.0*cosio) / temp4
	}

	s.aycof = -0.5 * g.j3oj2 * sinio

	s.delmo = 1.0 + s.eta*math.Cos(s.mo)
	s.delmo = s.delmo * s.delmo * s.delmo

	s.sinmao = math.Sin(s.mo)

	s.x7thm1 = 7.0*cosio2 - 1.0

	// DEEPSPACE
	if (twopi / s.no) >= 225.0 {
		s.method = "d"
		s.isimp = 1
		tc := 0.0
		inclm := s.inclo

		dscom(s, &c, tc)
		dpper(s, s.init, s.operationmode, s.ecco, s.inclo, s.nodeo, s.argpo, s.mo)
		dsinit(s, &c, tc, inclm)
	}

	// TODO
	// Test case for LEO use case where isimp == 1

	// NOT DEEPSPACE
	if s.isimp != 1 {
		cc1sq := s.cc1 * s.cc1
		s.d2 = 4.0 * ao * tsi * cc1sq
		temp := s.d2 * tsi * s.cc1 / 3.0
		s.d3 = (17.0*ao + sfour) * temp
		s.d4 = 0.5 * temp * ao * tsi * (221.0*ao + 31.0*sfour) * s.cc1
		s.t3cof = s.d2 + 2.0*cc1sq
		s.t4cof = 0.25 * (3.0*s.d3 + s.cc1*(12.0*s.d2+10.0*cc1sq))
		s.t5cof = 0.2 * (3.0*s.d4 + 12.0*s.cc1*s.d3 + 6.0*s.d2*s.d2 + 15.0*cc1sq*(2.0*s.d2+cc1sq))
	}
}
