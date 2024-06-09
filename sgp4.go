package sgp4

import (
	"fmt"
	"math"
	"time"
)

func Propogate(s *satrec, tprop time.Time) (p, v Vector) {
	j := JDay(tprop)
	tsince := (j - s.epoch_jday) * 1440
	return sgp4(s, tsince)
}

func sgp4(s *satrec, tsince float64) (p, v Vector) {
	var C cof

	mrt := 0.0
	vkmpersec := g.radius * g.xke / 60.0

	s.t = tsince

	xmdf := s.mo + s.mdot*tsince
	argpdf := s.argpo + s.argpdot*tsince
	nodedf := s.nodeo + s.nodedot*tsince
	argpm := argpdf
	mm := xmdf
	t2 := tsince * tsince
	nodem := nodedf + s.nodecf*t2
	tempa := 1.0 - s.cc1*tsince
	tempe := s.bstar * s.cc4 * tsince
	templ := s.t2cof * t2

	if s.isimp != 1 {
		delomg := s.omgcof * tsince
		delmtemp := 1.0 + s.eta*math.Cos(xmdf)
		delm := s.xmcof * (delmtemp*delmtemp*delmtemp - s.delmo)
		temp := delomg + delm
		mm = xmdf + temp
		argpm = argpdf - temp
		t3 := t2 * tsince
		t4 := t3 * tsince
		tempa = tempa - s.d2*t2 - s.d3*t3 - s.d4*t4
		tempe = tempe + s.bstar*s.cc5*(math.Sin(mm)-s.sinmao)
		templ = templ + s.t3cof*t3 + t4*(s.t4cof+tsince*s.t5cof)
	}

	nm := s.no
	em := s.ecco
	inclm := s.inclo

	if s.method == "d" {
		em, argpm, mm, inclm, nodem, nm = dspace(s, tsince, em, argpm, mm, inclm, nodem, nm)
	}

	if nm < 0.0 {
		fmt.Println("Error. Mean motion is less than zero")
	}

	am := math.Pow((g.xke/nm), x2o3) * tempa * tempa
	nm = g.xke / math.Pow(am, 1.5)
	em = em - tempe

	if em >= 1.0 || em < -0.001 {
		fmt.Println("Error. Mean eccentricity not within range 0 - 1")
	}

	if em < 1.0e-6 {
		em = 1.0e-6
	}

	mm = mm + s.no*templ
	xlm := mm + argpm + nodem
	C.emsq = em * em

	nodem = math.Mod(nodem, twopi)
	argpm = math.Mod(argpm, twopi)
	xlm = math.Mod(xlm, twopi)
	mm = math.Mod((xlm - argpm - nodem), twopi)

	C.sinim = math.Sin(inclm)
	C.cosim = math.Cos(inclm)

	ep := em
	xincp := inclm
	argpp := argpm
	nodep := nodem
	mp := mm
	sinip := C.sinim
	cosip := C.cosim

	if s.method == "d" {
		ep, xincp, nodep, argpp, mp = dpper(s, "n", s.operationmode, ep, xincp, nodep, argpp, mp)

		// TODO
		// replace dpper with dpper2 throughout the code base
		// for sgp4init use | dpper2(&G, init, opsmode, G.ecco, G.inclo, G.nodeo, G.argpo, G.mo)
		// need to adjust so that these are pass by reference not by copy. This is of little value here as we aren't passing complex variables by copy.

		if xincp < 0.0 {
			xincp = -xincp
			nodep += pi
			argpp -= pi
		}

		if ep < 0.0 || ep > 1.0 {
			fmt.Println("Perturbed eccentricity not within range 0 - 1")
		}

		sinip = math.Sin(xincp)
		cosip = math.Cos(xincp)
		s.aycof = -0.5 * g.j3oj2 * sinip

		if math.Abs(cosip+1.0) > 1.5e-12 {
			s.xlcof = -0.25 * g.j3oj2 * sinip * (3.0 + 5.0*cosip) / (1.0 + cosip)
		} else {
			s.xlcof = -0.25 * g.j3oj2 * sinip * (3.0 + 5.0*cosip) / temp4
		}
	}

	axnl := ep * math.Cos(argpp)
	temp := 1.0 / (am * (1.0 - ep*ep))
	aynl := ep*math.Sin(argpp) + temp*s.aycof
	xl := mp + argpp + nodep + temp*s.xlcof*axnl

	u := math.Mod((xl - nodep), twopi)
	eo1 := u
	tem5 := 9999.9
	ktr := 1

	var sineo1, coseo1 float64

	for math.Abs(tem5) >= 1.0e-12 && ktr <= 10 {
		sineo1 = math.Sin(eo1)
		coseo1 = math.Cos(eo1)
		tem5 = 1.0 - coseo1*axnl - sineo1*aynl
		tem5 = (u - aynl*coseo1 + axnl*sineo1 - eo1) / tem5
		if math.Abs(tem5) >= 0.95 {
			if tem5 > 0.0 {
				tem5 = 0.95
			} else {
				tem5 = -0.95
			}
		}
		eo1 = eo1 + tem5
		ktr = ktr + 1
	}

	ecose := axnl*coseo1 + aynl*sineo1
	esine := axnl*sineo1 - aynl*coseo1
	el2 := axnl*axnl + aynl*aynl
	pl := am * (1.0 - el2)

	if pl < 0.0 {
		fmt.Println("semilatus rectum is less than zero")
	}

	rl := am * (1.0 - ecose)
	rdotl := math.Sqrt(am) * esine / rl
	rvdotl := math.Sqrt(pl) / rl
	betal := math.Sqrt(1.0 - el2)
	temp = esine / (1.0 + betal)
	sinu := am / rl * (sineo1 - aynl - axnl*temp)
	cosu := am / rl * (coseo1 - axnl + aynl*temp)
	su := math.Atan2(sinu, cosu)
	sin2u := (cosu + cosu) * sinu
	cos2u := 1.0 - 2.0*sinu*sinu
	temp = 1.0 / pl
	temp1 := 0.5 * g.j2 * temp
	temp2 := temp1 * temp

	if s.method == "d" {
		cosisq := cosip * cosip
		s.con41 = 3.0*cosisq - 1.0
		s.x1mth2 = 1.0 - cosisq
		s.x7thm1 = 7.0*cosisq - 1.0
	}

	mrt = rl*(1.0-1.5*temp2*betal*s.con41) + 0.5*temp1*s.x1mth2*cos2u
	su = su - 0.25*temp2*s.x7thm1*sin2u
	xnode := nodep + 1.5*temp2*cosip*sin2u
	xinc := xincp + 1.5*temp2*cosip*sinip*cos2u
	mvt := rdotl - nm*temp1*s.x1mth2*sin2u/g.xke
	rvdot := rvdotl + nm*temp1*(s.x1mth2*cos2u+1.5*s.con41)/g.xke

	sinsu := math.Sin(su)
	cossu := math.Cos(su)
	snod := math.Sin(xnode)
	cnod := math.Cos(xnode)
	sini := math.Sin(xinc)
	cosi := math.Cos(xinc)

	xmx := -snod * cosi
	xmy := cnod * cosi
	ux := xmx*sinsu + cnod*cossu
	uy := xmy*sinsu + snod*cossu
	uz := sini * sinsu
	vx := xmx*cossu - cnod*sinsu
	vy := xmy*cossu - snod*sinsu
	vz := sini * cossu

	p.x = mrt * g.radius * ux
	p.y = mrt * g.radius * uy
	p.z = mrt * g.radius * uz

	v.x = (mvt*ux + rvdot*vx) * vkmpersec
	v.y = (mvt*uy + rvdot*vy) * vkmpersec
	v.z = (mvt*uz + rvdot*vz) * vkmpersec

	if mrt < 1.0 {
		fmt.Println("mrt is less than 1.0 indicating the satellite has decayed")
	}

	return
}
