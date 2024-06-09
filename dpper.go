package sgp4

import "math"

func dpper(s *satrec, init, opsmode string, ep, inclp, nodep, argpp, mp float64) (ep_o, inclp_o, nodep_o, argpp_o, mp_o float64) {
	zns := 1.19459e-5
	zes := 0.01675
	znl := 1.5835218e-4
	zel := 0.05490

	// Time varying periodics
	zm := s.zmos + zns*s.t

	if init == "y" {
		zm = s.zmos
	}

	zf := zm + 2.0*zes*math.Sin(zm)
	sinzf := math.Sin(zf)
	f2 := 0.5*sinzf*sinzf - 0.25
	f3 := -0.5 * sinzf * math.Cos(zf)
	ses := s.se2*f2 + s.se3*f3
	sis := s.si2*f2 + s.si3*f3
	sls := s.sl2*f2 + s.sl3*f3 + s.sl4*sinzf
	sghs := s.sgh2*f2 + s.sgh3*f3 + s.sgh4*sinzf
	shs := s.sh2*f2 + s.sh3*f3
	zm = s.zmol + znl*s.t

	if init == "y" {
		zm = s.zmol
	}

	zf = zm + 2.0*zel*math.Sin(zm)
	sinzf = math.Sin(zf)
	f2 = 0.5*sinzf*sinzf - 0.25
	f3 = -0.5 * sinzf * math.Cos(zf)
	sel := s.ee2*f2 + s.e3*f3
	sil := s.xi2*f2 + s.xi3*f3
	sll := s.xl2*f2 + s.xl3*f3 + s.xl4*sinzf
	sghl := s.xgh2*f2 + s.xgh3*f3 + s.xgh4*sinzf
	shll := s.xh2*f2 + s.xh3*f3
	pe := ses + sel
	pinc := sis + sil
	pl := sls + sll
	pgh := sghs + sghl
	ph := shs + shll

	if init == "n" {
		pe -= s.peo
		pinc -= s.pinco
		pl -= s.plo
		pgh -= s.pgho
		ph -= s.pho

		inclp += pinc
		ep += pe

		sinip := math.Sin(inclp)
		cosip := math.Cos(inclp)

		// TODO
		// Test case for inclo >= 0.2

		if inclp >= 0.2 {
			ph /= sinip
			pgh -= cosip * ph
			argpp += pgh
			nodep += ph
			mp += pl
		} else {
			// Apply periodics with lyddane modification
			sinop := math.Sin(nodep)
			cosop := math.Cos(nodep)
			alfdp := sinip * sinop
			betdp := sinip * cosop
			dalf := ph*cosop + pinc*cosip*sinop
			dbet := -ph*sinop + pinc*cosip*cosop
			alfdp = alfdp + dalf
			betdp = betdp + dbet
			nodep = math.Mod(nodep, twopi)
			//  sgp4fix for afspc written intrinsic functions
			// nodep used without a trigonometric function ahead
			// if ((nodep < 0.0) && (opsmode == 'a')) { o.nodep += twopi }
			xls := mp + argpp + pl + pgh + (cosip-pinc*sinip)*nodep
			xnoh := nodep
			nodep = math.Atan2(alfdp, betdp)
			//  sgp4fix for afspc written intrinsic functions
			// nodep used without a trigonometric function ahead
			// if ((nodep < 0.0) && (opsmode == 'a')) { nodep += twopi }
			if math.Abs(xnoh-nodep) > pi {
				if nodep < xnoh {
					nodep += twopi
				} else {
					nodep -= twopi
				}
			}
			mp += pl
			argpp = xls - mp - cosip*nodep
		}
	}
	return ep, inclp, nodep, argpp, mp
}
