package sgp4

import "math"

func dspace(s *satrec, tc, em, argpm, mm, inclm, nodem, nm float64) (em_o, argpm_o, mm_o, inclm_o, nodem_o, nm_o float64) {

	var delt, xndt, xldot, xnddt, xomi, x2omi, x2li, xl float64

	fasx2 := 0.13130908
	fasx4 := 2.8843198
	fasx6 := 0.37448087
	g22 := 5.7686396
	g32 := 0.95240898
	g44 := 1.8014998
	g52 := 1.0508330
	g54 := 4.4108898
	rptim := 4.37526908801129966e-3
	stepp := 720.0
	stepn := -720.0
	step2 := 259200.0

	// Deep Space Resonance Effects
	s.dndt = 0.0
	theta := math.Mod((s.gsto + tc*rptim), twopi)
	em = em + s.dedt*s.t

	inclm = inclm + s.didt*s.t
	argpm = argpm + s.domdt*s.t
	nodem = nodem + s.dnodt*s.t
	mm = mm + s.dmdt*s.t

	// sgp4fix for negative inclinations
	// the following if statement should be commented out
	if inclm < 0.0 {
		inclm = -inclm
		argpm = argpm - pi
		nodem = nodem + pi
	}

	// update resonance: numerical (euler-maclaurin) integration
	// epoch restart --------------------------------------------
	// sgp4fix for propagator problems
	// the following integration works for neative time steps and periods
	// the specific changes are unknown because the original code was so convoluted

	// sgp4fix take out atime = 0.0 and fix for faster operation
	ft := 0.0

	if s.irez != 0.0 {
		// spg4fix streamline check
		if s.atime == 0.0 || s.t*s.atime <= 0.0 || math.Abs(s.t) < math.Abs(s.atime) {
			s.atime = 0.0
			s.xni = s.no
			s.xli = s.xlamo
		}
		// sgp4fix move check outside loop
		if s.t > 0.0 {
			delt = stepp
		} else {
			delt = stepn
		}

		j := 0

		for j == 0 {
			if s.irez != 2 {
				// near-synchronous resonance terms
				xndt = s.del1*math.Sin(s.xli-fasx2) +
					s.del2*math.Sin(2.0*(s.xli-fasx4)) +
					s.del3*math.Sin(3.0*(s.xli-fasx6))
				xldot = s.xni + s.xfact
				xnddt = s.del1*math.Cos(s.xli-fasx2) +
					2.0*s.del2*math.Cos(2.0*(s.xli-fasx4)) +
					3.0*s.del3*math.Cos(3.0*(s.xli-fasx6))
				xnddt = xnddt * xldot
			} else {
				// near-half-day resonance terms
				xomi = s.argpo + s.argpdot*s.atime
				x2omi = xomi + xomi
				x2li = s.xli + s.xli
				xndt = (s.d2201*math.Sin(x2omi+s.xli-g22) +
					s.d2211*math.Sin(s.xli-g22) +
					s.d3210*math.Sin(xomi+s.xli-g32) +
					s.d3222*math.Sin(-xomi+s.xli-g32) +
					s.d4410*math.Sin(x2omi+x2li-g44) +
					s.d4422*math.Sin(x2li-g44) +
					s.d5220*math.Sin(xomi+s.xli-g52) +
					s.d5232*math.Sin(-xomi+s.xli-g52) +
					s.d5421*math.Sin(xomi+x2li-g54) +
					s.d5433*math.Sin(-xomi+x2li-g54))
				xldot = s.xni + s.xfact
				xnddt = (s.d2201*math.Cos(x2omi+s.xli-g22) +
					s.d2211*math.Cos(s.xli-g22) +
					s.d3210*math.Cos(xomi+s.xli-g32) +
					s.d3222*math.Cos(-xomi+s.xli-g32) +
					s.d5220*math.Cos(xomi+s.xli-g52) +
					s.d5232*math.Cos(-xomi+s.xli-g52) +
					2.0*(s.d4410*math.Cos(x2omi+x2li-g44)+
						s.d4422*math.Cos(x2li-g44)+
						s.d5421*math.Cos(xomi+x2li-g54)+
						s.d5433*math.Cos(-xomi+x2li-g54)))
				xnddt = xnddt * xldot
			}
			// integrator
			if math.Abs(s.t-s.atime) >= stepp {
				j = 0
			} else {
				ft = s.t - s.atime
				j = 1
			}

			if j == 0 {
				s.xli = s.xli + xldot*delt + xndt*step2
				s.xni = s.xni + xndt*delt + xnddt*step2
				s.atime += delt
			}
		}

		nm = s.xni + xndt*ft + xnddt*ft*ft*0.5
		xl = s.xli + xldot*ft + xndt*ft*ft*0.5

		if s.irez != 1 {
			mm = xl - 2.0*nodem + 2.0*theta
			s.dndt = nm - s.no
		} else {
			mm = xl - nodem - argpm + theta
			s.dndt = nm - s.no
		}

		nm = s.no + s.dndt
	}

	return em, argpm, mm, inclm, nodem, nm
}
