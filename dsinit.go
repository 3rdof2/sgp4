package sgp4

import "math"

func dsinit(s *satrec, c *cof, tc, inclm float64) (em, argpm, mm, nm, inclm_o, nodem float64) {

	q22 := 1.7891679e-6
	q31 := 2.1460748e-6
	q33 := 2.2123015e-7
	root22 := 1.7891679e-6
	root44 := 7.3636953e-9
	root54 := 2.1765803e-9
	rptim := 4.37526908801129966e-3
	root32 := 3.7393792e-7
	root52 := 1.1428639e-7
	znl := 1.5835218e-4
	zns := 1.19459e-5

	// Deep space initialization
	s.irez = 0
	if s.no < 0.0052359877 && s.no > 0.0034906585 {
		s.irez = 1
	}
	if s.no >= 8.26e-3 && s.no <= 9.24e-3 && s.ecco >= 0.5 {
		s.irez = 2
	}

	// Solar Terms
	ses := c.ss1 * zns * c.ss5
	sis := c.ss2 * zns * (c.sz11 + c.sz13)
	sls := -zns * c.ss3 * (c.sz1 + c.sz3 - 14.0 - 6.0*c.eccsq)
	sghs := c.ss4 * zns * (c.sz31 + c.sz33 - 6.0)
	shs := -zns * c.ss2 * (c.sz21 + c.sz23)

	if s.inclo < 5.2359877e-2 || s.inclo > twopi/2-5.2359877e-2 {
		shs = 0.0
	}
	if c.sinim != 0.0 {
		shs = shs / c.sinim
	}

	sgs := sghs - c.cosim*shs

	// Lunar Terms
	s.dedt = ses + c.s1*znl*c.s5
	s.didt = sis + c.s2*znl*(c.z11+c.z13)
	s.dmdt = sls - znl*c.s3*(c.z1+c.z3-14.0-6.0*c.emsq)
	sghl := c.s4 * znl * (c.z31 + c.z33 - 6.0)
	shll := -znl * c.s2 * (c.z21 + c.z23)

	// sgp4fix for 180 deg incl
	if s.inclo < 5.2359877e-2 || s.inclo > twopi/2-5.2359877e-2 {
		shll = 0.0
	}

	s.domdt = sgs + sghl
	s.dnodt = shs

	if c.sinim != 0.0 {
		s.domdt -= c.cosim / c.sinim * shll
		s.dnodt += shll / c.sinim
	}

	// Deep Space resonsnace effects
	s.dndt = 0.0
	theta := math.Mod(s.gsto+tc*rptim, twopi)
	em = s.ecco + s.dedt*s.t
	inclm += s.didt * s.t
	argpm = s.domdt * s.t
	nodem = s.dnodt * s.t
	mm = s.dmdt * s.t

	// INITIALIZE THE RESONANCE TERMS
	var g201, g211, g310, g322, g410, g422, g520, g533, g521, g532, f220, f221, f321, f322, f441, f442 float64
	var f522, f523, f542, f543 float64

	if s.irez != 0 {
		aonv := math.Pow(s.no/g.xke, x2o3)

		// TODO
		// Implement test case where irez == 2
		// Conditions required for irez == 2
		//	s.no >= 8.26e-3 &&
		//	s.no <= 9.24e-3 &&
		//	s.ecco >= 0.5

		// Geopotential resonance for 12 hour orbits
		if s.irez == 2 {
			cosisq := c.cosim * c.cosim
			emo := em
			em = s.ecco
			emsqo := c.emsq
			c.emsq = c.eccsq
			eoc := em * c.emsq
			g201 = -0.306 - (em-0.64)*0.440

			if em <= 0.65 {
				g211 = 3.616 - 13.2470*em + 16.2900*c.emsq
				g310 = -19.302 + 117.3900*em - 228.4190*c.emsq + 156.5910*eoc
				g322 = -18.9068 + 109.7927*em - 214.6334*c.emsq + 146.5816*eoc
				g410 = -41.122 + 242.6940*em - 471.0940*c.emsq + 313.9530*eoc
				g422 = -146.407 + 841.8800*em - 1629.014*c.emsq + 1083.4350*eoc
				g520 = -532.114 + 3017.977*em - 5740.032*c.emsq + 3708.2760*eoc
			} else {
				g211 = -72.099 + 331.819*em - 508.738*c.emsq + 266.724*eoc
				g310 = -346.844 + 1582.851*em - 2415.925*c.emsq + 1246.113*eoc
				g322 = -342.585 + 1554.908*em - 2366.899*c.emsq + 1215.972*eoc
				g410 = -1052.797 + 4758.686*em - 7193.992*c.emsq + 3651.957*eoc
				g422 = -3581.690 + 16178.110*em - 24462.770*c.emsq + 12422.520*eoc
				if em > 0.715 {
					g520 = -5149.66 + 29936.92*em - 54087.36*c.emsq + 31324.56*eoc
				} else {
					g520 = 1464.74 - 4664.75*em + 3763.64*c.emsq
				}
			}

			if em < 0.7 {
				g533 = -919.22770 + 4988.6100*em - 9064.7700*c.emsq + 5542.21*eoc
				g521 = -822.71072 + 4568.6173*em - 8491.4146*c.emsq + 5337.524*eoc
				g532 = -853.66600 + 4690.2500*em - 8624.7700*c.emsq + 5341.4*eoc
			} else {
				g533 = -37995.780 + 161616.52*em - 229838.20*c.emsq + 109377.94*eoc
				g521 = -51752.104 + 218913.95*em - 309468.16*c.emsq + 146349.42*eoc
				g532 = -40023.880 + 170470.89*em - 242699.48*c.emsq + 115605.82*eoc
			}

			sini2 := c.sinim * c.sinim
			f220 = 0.75 * (1.0 + 2.0*c.cosim + cosisq)
			f221 = 1.5 * sini2
			f321 = 1.875 * c.sinim * (1.0 - 2.0*c.cosim - 3.0*cosisq)
			f322 = -1.875 * c.sinim * (1.0 + 2.0*c.cosim - 3.0*cosisq)
			f441 = 35.0 * sini2 * f220
			f442 = 39.3750 * sini2 * sini2
			f522 = 9.84375 * c.sinim * (sini2*(1.0-2.0*c.cosim-5.0*cosisq) + 0.33333333*(-2.0+4.0*c.cosim+6.0*cosisq))
			f523 = c.sinim * (4.92187512*sini2*(-2.0-4.0*c.cosim+10.0*cosisq) + 6.56250012*(1.0+2.0*c.cosim-3.0*cosisq))
			f542 = 29.53125 * c.sinim * (2.0 - 8.0*c.cosim + cosisq*(-12.0+8.0*c.cosim+10.0*cosisq))
			f543 = 29.53125 * c.sinim * (-2.0 - 8.0*c.cosim + cosisq*(12.0+8.0*c.cosim-10.0*cosisq))
			xno2 := s.no * s.no
			ainv2 := aonv * aonv
			temp1 := 3.0 * xno2 * ainv2
			temp := temp1 * root22
			s.d2201 = temp * f220 * g201
			s.d2211 = temp * f221 * g211
			temp1 = temp1 * aonv
			temp = temp1 * root32
			s.d3210 = temp * f321 * g310
			s.d3222 = temp * f322 * g322
			temp1 = temp1 * aonv
			temp = 2.0 * temp1 * root44
			s.d4410 = temp * f441 * g410
			s.d4422 = temp * f442 * g422
			temp1 = temp1 * aonv
			temp = temp1 * root52
			s.d5220 = temp * f522 * g520
			s.d5232 = temp * f523 * g532
			temp = 2.0 * temp1 * root54
			s.d5421 = temp * f542 * g521
			s.d5433 = temp * f543 * g533
			s.xlamo = math.Mod(s.mo+s.nodeo+s.nodeo-theta-theta, twopi)
			s.xfact = s.mdot + s.dmdt + 2.0*(s.nodedot+s.dnodt-rptim) - s.no
			em = emo
			c.emsq = emsqo
		}

		if s.irez == 1 {
			g200 := 1.0 + c.emsq*(-2.5+0.8125*c.emsq)
			g310 = 1.0 + 2.0*c.emsq
			g300 := 1.0 + c.emsq*(-6.0+6.60937*c.emsq)
			f220 = 0.75 * (1.0 + c.cosim) * (1.0 + c.cosim)
			f311 := 0.9375*c.sinim*c.sinim*(1.0+3.0*c.cosim) - 0.75*(1.0+c.cosim)
			f330 := 1.0 + c.cosim
			f330 = 1.875 * f330 * f330 * f330
			s.del1 = 3.0 * s.no * s.no * aonv * aonv
			s.del2 = 2.0 * s.del1 * f220 * g200 * q22
			s.del3 = 3.0 * s.del1 * f330 * g300 * q33 * aonv
			s.del1 = s.del1 * f311 * g310 * q31 * aonv
			s.xlamo = math.Mod(s.mo+s.nodeo+s.argpo-theta, twopi) // << ERROR
			s.xfact = s.mdot + c.xpidot - rptim + s.dmdt + s.domdt + s.dnodt - s.no
		}
		s.xli = s.xlamo
		s.xni = s.no
		s.atime = 0.0
		nm = s.no + s.dndt
	}
	return em, argpm, mm, nm, inclm, nodem
}
