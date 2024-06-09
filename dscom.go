package sgp4

import "math"

func dscom(s *satrec, c *cof, tc float64) {

	zes := 0.01675
	zel := 0.05490
	c1ss := 2.9864797e-6
	c1l := 4.7968065e-7
	zsinis := 0.39785416
	zcosis := 0.91744867
	zcosgs := 0.1945905
	zsings := -0.98088458

	snodm := math.Sin(s.nodeo)
	cnodm := math.Cos(s.nodeo)
	sinomm := math.Sin(s.argpo)
	cosomm := math.Cos(s.argpo)
	c.sinim = math.Sin(s.inclo)
	c.cosim = math.Cos(s.inclo)
	c.emsq = s.ecco * s.ecco
	betasq := 1.0 - c.emsq
	rtemsq := math.Sqrt(betasq)

	day := s.epochF + 18261.5 + tc/1440.0
	xnodce := math.Mod(4.5236020-9.2422029e-4*day, twopi)
	stem := math.Sin(xnodce)
	ctem := math.Cos(xnodce)
	zcosil := 0.91375164 - 0.03568096*ctem
	zsinil := math.Sqrt(1.0 - zcosil*zcosil)
	zsinhl := 0.089683511 * stem / zsinil
	zcoshl := math.Sqrt(1.0 - zsinhl*zsinhl)
	gam := 5.8351514 + 0.0019443680*day
	zx := 0.39785416 * stem / zsinil
	zy := zcoshl*ctem + 0.91744867*zsinhl*stem
	zx = math.Atan2(zx, zy)
	zx = gam + zx - xnodce
	zcosgl := math.Cos(zx)
	zsingl := math.Sin(zx)

	// DO SOLAR TERMS
	zcosg := zcosgs
	zsing := zsings
	zcosi := zcosis
	zsini := zsinis
	zcosh := cnodm
	zsinh := snodm
	cc := c1ss
	xnoi := 1.0 / s.no

	var s6, s7, ss6, ss7, sz2, sz12, sz22, sz32, z2, z12, z22, z32 float64

	for i := 1; i <= 2; i++ {
		a1 := zcosg*zcosh + zsing*zcosi*zsinh
		a3 := -zsing*zcosh + zcosg*zcosi*zsinh
		a7 := -zcosg*zsinh + zsing*zcosi*zcosh
		a8 := zsing * zsini
		a9 := zsing*zsinh + zcosg*zcosi*zcosh
		a10 := zcosg * zsini
		a2 := c.cosim*a7 + c.sinim*a8
		a4 := c.cosim*a9 + c.sinim*a10
		a5 := -c.sinim*a7 + c.cosim*a8
		a6 := -c.sinim*a9 + c.cosim*a10

		x1 := a1*cosomm + a2*sinomm
		x2 := a3*cosomm + a4*sinomm
		x3 := -a1*sinomm + a2*cosomm
		x4 := -a3*sinomm + a4*cosomm
		x5 := a5 * sinomm
		x6 := a6 * sinomm
		x7 := a5 * cosomm
		x8 := a6 * cosomm

		c.z31 = 12.0*x1*x1 - 3.0*x3*x3
		z32 = 24.0*x1*x2 - 6.0*x3*x4
		c.z33 = 12.0*x2*x2 - 3.0*x4*x4
		c.z1 = 3.0*(a1*a1+a2*a2) + c.z31*c.emsq
		z2 = 6.0*(a1*a3+a2*a4) + z32*c.emsq
		c.z3 = 3.0*(a3*a3+a4*a4) + c.z33*c.emsq
		c.z11 = -6.0*a1*a5 + c.emsq*(-24.0*x1*x7-6.0*x3*x5)
		z12 = -6.0*(a1*a6+a3*a5) + c.emsq*(-24.0*(x2*x7+x1*x8)-6.0*(x3*x6+x4*x5))
		c.z13 = -6.0*a3*a6 + c.emsq*(-24.0*x2*x8-6.0*x4*x6)
		c.z21 = 6.0*a2*a5 + c.emsq*(24.0*x1*x5-6.0*x3*x7)
		z22 = 6.0*(a4*a5+a2*a6) + c.emsq*(24.0*(x2*x5+x1*x6)-6.0*(x4*x7+x3*x8))
		c.z23 = 6.0*a4*a6 + c.emsq*(24.0*x2*x6-6.0*x4*x8)
		c.z1 = c.z1 + c.z1 + betasq*c.z31
		z2 = z2 + z2 + betasq*z32
		c.z3 = c.z3 + c.z3 + betasq*c.z33
		c.s3 = cc * xnoi
		c.s2 = -0.5 * c.s3 / rtemsq
		c.s4 = c.s3 * rtemsq
		c.s1 = -15.0 * s.ecco * c.s4
		c.s5 = x1*x3 + x2*x4
		s6 = x2*x3 + x1*x4
		s7 = x2*x4 - x1*x3

		if i == 1 {
			c.ss1 = c.s1
			c.ss2 = c.s2
			c.ss3 = c.s3
			c.ss4 = c.s4
			c.ss5 = c.s5
			ss6 = s6
			ss7 = s7
			c.sz1 = c.z1
			sz2 = z2
			c.sz3 = c.z3
			c.sz11 = c.z11
			sz12 = z12
			c.sz13 = c.z13
			c.sz21 = c.z21
			sz22 = z22
			c.sz23 = c.z23
			c.sz31 = c.z31
			sz32 = z32
			c.sz33 = c.z33
			zcosg = zcosgl
			zsing = zsingl
			zcosi = zcosil
			zsini = zsinil
			zcosh = zcoshl*cnodm + zsinhl*snodm
			zsinh = snodm*zcoshl - cnodm*zsinhl
			cc = c1l
		}
	}

	s.zmol = math.Mod(4.7199672+0.22997150*day-gam, twopi)
	s.zmos = math.Mod(6.2565837+0.017201977*day, twopi)

	s.se2 = 2.0 * c.ss1 * ss6
	s.se3 = 2.0 * c.ss1 * ss7
	s.si2 = 2.0 * c.ss2 * sz12
	s.si3 = 2.0 * c.ss2 * (c.sz13 - c.sz11)
	s.sl2 = -2.0 * c.ss3 * sz2
	s.sl3 = -2.0 * c.ss3 * (c.sz3 - c.sz1)
	s.sl4 = -2.0 * c.ss3 * (-21.0 - 9.0*c.emsq) * zes
	s.sgh2 = 2.0 * c.ss4 * sz32
	s.sgh3 = 2.0 * c.ss4 * (c.sz33 - c.sz31)
	s.sgh4 = -18.0 * c.ss4 * zes
	s.sh2 = -2.0 * c.ss2 * sz22
	s.sh3 = -2.0 * c.ss2 * (c.sz23 - c.sz21)

	// LUNAR TERMS
	s.ee2 = 2.0 * c.s1 * s6
	s.e3 = 2.0 * c.s1 * s7
	s.xi2 = 2.0 * c.s2 * z12
	s.xi3 = 2.0 * c.s2 * (c.z13 - c.z11)
	s.xl2 = -2.0 * c.s3 * z2
	s.xl3 = -2.0 * c.s3 * (c.z3 - c.z1)
	s.xl4 = -2.0 * c.s3 * (-21.0 - 9.0*c.emsq) * zel
	s.xgh2 = 2.0 * c.s4 * z32
	s.xgh3 = 2.0 * c.s4 * (c.z33 - c.z31)
	s.xgh4 = -18.0 * c.s4 * zel
	s.xh2 = -2.0 * c.s2 * z22
	s.xh3 = -2.0 * c.s2 * (c.z23 - c.z21)
}
