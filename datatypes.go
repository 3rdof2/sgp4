package sgp4

import "time"

type satrec struct {
	// OMM INPUT
	satnum                                                string
	ndot, nddot, bstar, inclo, nodeo, ecco, argpo, mo, no float64
	epoch                                                 time.Time
	epoch_jday, epochF                                    float64
	// FLAGS
	method, operationmode, init string
	// NEAR SPACE COEFFICIENT
	argpdot, aycof, cc1, cc4, cc5, con41                                        float64
	d2, d3, d4, delmo, eta, gsto, mdot, nodecf, nodedot                         float64
	omgcof, sinmao, t, t2cof, t3cof, t4cof, t5cof, x1mth2, x7thm1, xlcof, xmcof float64
	isimp                                                                       int
	// DEEP SPACE COEFFICIENT
	atime, e3, ee2, d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421 float64
	d5433, del1, del2, del3, dedt, didt, domdt, dmdt, dndt, dnodt                 float64
	peo, pinco, pgho, plo, pho                                                    float64
	se2, se3, sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4                 float64
	xl2, xl3, xl4, xlamo, xli, xfact, xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xni   float64
	zmol, zmos                                                                    float64
	irez                                                                          int
}

type cof struct {
	sinim, cosim, emsq, eccsq, xpidot            float64
	s1, s2, s3, s4, s5                           float64
	ss1, ss2, ss3, ss4, ss5                      float64
	sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33 float64
	z1, z3, z11, z13, z21, z23, z31, z33         float64
}

type OMM struct {
	NORAD_CAT_ID      string  // satnum
	EPOCH             string  // epoch
	MEAN_MOTION_DOT   float64 // ndot
	MEAN_MOTION_DDOT  float64 // nddot
	BSTAR             float64 // bstar
	INCLINATION       float64 // inclo
	RA_OF_ASC_NODE    float64 // nodeo
	ECCENTRICITY      float64 // ecco
	ARG_OF_PERICENTER float64 // argpo
	MEAN_ANOMALY      float64 // mo
	MEAN_MOTION       float64 // no
}

type Vector struct{ x, y, z float64 }

type LLA struct {
	Lat  float64
	Long float64
	Alt  float64
}
