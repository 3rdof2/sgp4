package sgp4

import (
	"fmt"
	"time"
)

func Convert(x interface{}) ([]float64, []string, []int, time.Time) {
	switch x := x.(type) {
	case *satrec:
		// convert the satrec
		aString := []string{x.satnum, x.method, x.operationmode, x.init}
		aFloat := []float64{x.ndot, x.nddot, x.bstar, x.inclo, x.nodeo, x.ecco, x.argpo, x.mo, x.no, x.epoch_jday, x.epochF,
			x.argpdot, x.aycof, x.cc1, x.cc4, x.cc5, x.con41, x.d2, x.d3, x.d4, x.delmo, x.eta,
			x.gsto, x.mdot, x.nodecf, x.nodedot, x.omgcof, x.sinmao, x.t, x.t2cof, x.t3cof, x.t4cof, x.t5cof,
			x.x1mth2, x.x7thm1, x.xlcof, x.xmcof, x.atime, x.e3, x.ee2, x.d2201, x.d2211, x.d3210, x.d3222, x.d4410,
			x.d4422, x.d5220, x.d5232, x.d5421, x.d5433, x.del1, x.del2, x.del3, x.dedt, x.didt, x.domdt, x.dmdt,
			x.dndt, x.dnodt, x.peo, x.pinco, x.pgho, x.plo, x.pho, x.se2, x.se3, x.sgh2, x.sgh3, x.sgh4, x.sh2, x.sh3,
			x.si2, x.si3, x.sl2, x.sl3, x.sl4, x.xl2, x.xl3, x.xl4, x.xlamo, x.xli, x.xfact, x.xgh2, x.xgh3, x.xgh4,
			x.xh2, x.xh3, x.xi2, x.xi3, x.xni, x.zmol, x.zmos}
		return aFloat, aString, []int{x.isimp, x.irez}, x.epoch
	case cof:
		// convert the cof
		aFloat := []float64{x.sinim, x.cosim, x.emsq, x.eccsq, x.xpidot, x.s1, x.s2, x.s3, x.s4, x.s5, x.ss1, x.ss2,
			x.ss3, x.ss4, x.ss5, x.sz1, x.sz3, x.sz11, x.sz13, x.sz21, x.sz23, x.sz31, x.sz33, x.z1, x.z3, x.z11,
			x.z13, x.z21, x.z23, x.z31, x.z33}
		return aFloat, nil, nil, time.Time{}
	case Vector:
		aFloat := []float64{x.x, x.y, x.z}
		return aFloat, nil, nil, time.Time{}
	case LLA:
		aFloat := []float64{x.Alt, x.Lat, x.Long}
		return aFloat, nil, nil, time.Time{}
	}
	// If an inappropriate data type was given, the function will return a nil answer
	return nil, nil, nil, time.Time{}
}

var sat []string = []string{"ndot", "nddot", "bstar", "inclo", "nodeo", "ecco", "argpo", "mo", "no", "epoch_jday",
	"epochF", "argpdot", "aycof", "cc1", "cc4", "cc5", "con41", "d2", "d3", "d4", "delmo", "eta",
	"gsto", "mdot", "nodecf", "nodedot", "omgcof", "sinmao", "t", "t2cof", "t3cof", "t4cof", "t5cof",
	"x1mth2", "x7thm1", "xlcof", "xmcof", "atime", "e3", "ee2", "d2201", "d2211", "d3210", "d3222", "d4410",
	"d4422", "d5220", "d5232", "d5421", "d5433", "del1", "del2", "del3", "dedt", "didt", "domdt", "dmdt",
	"dndt", "dnodt", "peo", "pinco", "pgho", "plo", "pho", "se2", "se3", "sgh2", "sgh3", "sgh4", "sh2", "sh3",
	"si2", "si3", "sl2", "sl3", "sl4", "xl2", "xl3", "xl4", "xlamo", "xli", "xfact", "xgh2", "xgh3", "xgh4",
	"xh2", "xh3", "xi2", "xi3", "xni", "zmol", "zmos"}
var satString []string = []string{"satnum", "method", "operationsmode", "init"}
var satInt []string = []string{"isimp", "irez"}

func Check(a, b interface{}) ([]string, error) {

	if fmt.Sprintf("%T", a) != fmt.Sprintf("%T", b) {
		return nil, fmt.Errorf("a is not the same type as b.\n")
	}

	f, s, i, t := Convert(a)
	F, S, I, T := Convert(b)

	var resp []string

	for j := 0; j < len(f); j++ {
		if f[j] != F[j] {
			e := (f[j] - F[j]) / f[j]
			if e > 0.000001 {
				resp = append(resp, fmt.Sprintf("%v, a: %v, b: %v, err: %v", sat[j], f[j], F[j], e))
			}
		}
	}
	for j := 0; j < len(s); j++ {
		if s[j] != S[j] {
			resp = append(resp, fmt.Sprintf("%v : %v != %v", satString[j], s[j], S[j]))
		}
	}
	for j := 0; j < len(i); j++ {
		if i[j] != I[j] {
			resp = append(resp, fmt.Sprintf("%v, %v != %v", satInt[j], i[j], I[j]))
		}
	}
	if t != T {
		resp = append(resp, fmt.Sprintf("a: %v, b: %v", t, T))
	}

	return resp, nil
}
