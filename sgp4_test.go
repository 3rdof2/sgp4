package sgp4

import (
	"math"
	"reflect"
	"testing"
	"time"
)

func TestTLEinit(t *testing.T) {
	tle1 := "1 36867U 82025BE  24004.89010392  .00003744  00000-0  18028-2 0  9998"
	tle2 := "2 36867  82.5273 228.2287 0022169 119.6504 240.6886 14.18411700142970"

	var o OMM

	o.NORAD_CAT_ID = "36867"
	o.EPOCH = "2024-01-04T21:21:44.978688"
	o.MEAN_MOTION = 14.18411700
	o.ECCENTRICITY = 0.00221690
	o.MEAN_MOTION_DOT = 0.00003744
	o.MEAN_MOTION_DDOT = 0.0000000000000
	o.BSTAR = 0.00180280000000
	o.INCLINATION = 82.5273
	o.RA_OF_ASC_NODE = 228.2287
	o.ARG_OF_PERICENTER = 119.6504
	o.MEAN_ANOMALY = 240.6886

	s := OMMinit(&o)
	a := TLEinit(tle1, tle2)

	str, err := Check(s, a)
	if err != nil {
		t.Errorf("%v", err)
	}
	for _, v := range str {
		t.Errorf("%v", v)
	}

}

func TestSgp4init(t *testing.T) {
	var s satrec
	s.satnum = "44071"
	s.epoch, _ = time.Parse("2006-01-02T15:04:05.999999", "2023-11-30T14:51:58.008096")
	s.ndot = 1.6665470288140301e-12
	s.nddot = 0
	s.bstar = 0
	s.inclo = 0.0002705260340591211
	s.nodeo = 6.191936003226819
	s.ecco = 3.15e-05
	s.argpo = 4.442295787980063
	s.mo = 1.801858721137174
	s.no = 0.00437515764041159

	var b satrec
	b = satrec{
		satnum:        "44071",
		epoch_jday:    2.4602791194212963e+06,
		ndot:          1.6665470288140301e-12,
		nddot:         0,
		bstar:         0,
		inclo:         0.0002705260340591211,
		nodeo:         6.191936003226819,
		ecco:          3.15e-05,
		argpo:         4.442295787980063,
		mo:            1.801858721137174,
		no:            0.004374995068569369,
		method:        "d",
		operationmode: "i",
		init:          "y",
		gsto:          5.100098369122165,
		isimp:         1,
		con41:         1.9999997804470002,
		cc5:           2.5301454032485523e-11,
		d4:            0,
		argpdot:       3.251834409309059e-07,
		t:             0,
		t4cof:         0,
		x7thm1:        5.999999487709667,
		xlcof:         6.327307763430076e-07,
		cc1:           0,
		d2:            0,
		delmo:         0.9999744459150449,
		omgcof:        -0,
		t2cof:         0,
		t5cof:         0,
		mdot:          0.004375157637389939,
		xmcof:         0,
		aycof:         3.1636538961856574e-07,
		cc4:           2.9404363505917423e-17,
		d3:            0,
		eta:           3.719502845151615e-05,
		sinmao:        0.9734236435400782,
		t3cof:         0,
		x1mth2:        7.318433326020113e-08,
		nodedot:       -1.625887093095899e-07,
		nodecf:        -0,
		irez:          1,
		d3210:         0,
		d4422:         0,
		d5421:         0,
		del1:          -6.397477159775755e-13,
		didt:          -4.3587945820787964e-09,
		domdt:         2.715064341617845e-08,
		peo:           0,
		pinco:         0,
		se3:           5.326430990719636e-07,
		sgh4:          -0.0002058113472108935,
		si2:           -0.0014469831366815606,
		sl3:           0.00623549779429091,
		xfact:         -0.004375014923455862,
		xgh4:          -0.0001083476462637093,
		xi2:           0.00021383874647060809,
		xl3:           -0.0025975660572708273,
		zmol:          2.2017398948110127,
		xli:           1.0528068360423042,
		d2201:         0,
		d3222:         0,
		d5220:         0,
		d5433:         0,
		del2:          1.4104020776079415e-11,
		dmdt:          -9.814964016877237e-08,
		e3:            7.644510976836127e-09,
		pgho:          0,
		plo:           0,
		sgh2:          -0.01651601381446915,
		sh2:           -0.0007051137583276835,
		si3:           0.0007448223807063925,
		sl4:           0.0004802264772678869,
		xgh2:          -8.477119050457478e-05,
		xh2:           -0.00019163490436296422,
		xi3:           0.0002249974031644886,
		xl4:           0.0002528111748482556,
		zmos:          5.687936935123361,
		xni:           0.004374995068569369,
		d2211:         0,
		d4410:         0,
		d5232:         0,
		dedt:          9.025642197376453e-13,
		del3:          1.978470220518558e-12,
		dnodt:         0,
		ee2:           9.138261842041732e-08,
		pho:           0,
		se2:           -2.603293268959062e-07,
		sgh3:          -0.007432491692699144,
		sh3:           -0.0013187293886349073,
		sl2:           0.01602180512327805,
		xgh3:          0.0027676424659859498,
		xh3:           0.0001976717311208955,
		xl2:           -0.0001558820432923891,
		xlamo:         1.0528068360423042,
		atime:         0,
	}

	sgp4init(&s)

	values := reflect.ValueOf(s)
	values_b := reflect.ValueOf(b)
	types_b := values.Type()

	for i := 0; i < values.NumField(); i++ {
		if i == 0 || i == 13 || i == 14 || i == 15 { // satnum method, operationsmode, init
			b := values_b.Field(i).String()
			a := values.Field(i).String()
			if a != b {
				t.Errorf("%v %v %v", types_b.Field(i).Name, b, a)
			}
		} else if i == 10 || i == 12 || i == 42 || i == 98 { // epoch, epochF, isimp, irez
		} else if values_b.Field(i).Float() != values.Field(i).Float() {
			b := values_b.Field(i).Float()
			a := values.Field(i).Float()
			e := math.Abs((a-b)/b) * 100
			if e >= 0.001 {
				t.Errorf("%v | %v n: %v err: %v", types_b.Field(i).Name, b, a, e)
			}
		}
	}
}

func TestInitialize(t *testing.T) {
	decoder := []string{"ainv", "G.no", "ao", "G.con41", "con42", "cosio", "cosio2", "eccsq", "omeosq", "posq", "rp", "rteosq", "sinio", "G.ßgsto"}
	var G satrec
	G_ans := []float64{}
	G_base := []float64{0.15126279075433127, 0.004374995068569369, 6.611011174745008, 1.9999997804470002, -3.9999996340783337,
		0.9999999634078327, 0.9999999268156667, 9.9225e-10, 0.99999999900775, 43.70546866586987, 6.610802927893004, 0.999999999503875,
		0.00027052603075940977, 5.100098369122165}
	G.ecco = 3.15e-05
	G.inclo = 0.0002705260340591211
	G.no = 0.00437515764041159
	G.epoch, _ = time.Parse("2006-01-02T15:04:05.999999", "2023-11-30T14:51:58.008096")
	ainv, ao, con42, cosio, cosio2, eccsq, omeosq, posq, rp, rteosq, sinio := initl(&G)

	G_ans = append(G_ans, ainv, G.no, ao, G.con41, con42, cosio, cosio2, eccsq, omeosq, posq, rp, rteosq, sinio, G.gsto)

	for i, v := range G_ans {
		e := math.Abs((v-G_base[i])/v) * 100
		if e >= 0.001 {
			t.Errorf("%v, b:%v, a:%v, error: %v", decoder[i], G_base[i], v, e*100)
		}
	}
}

func TestDscom(t *testing.T) {
	decoder := []string{"c.sinim", "c.cosim", "G.e3", "G.ee2", "G.ecco", "c.emsq", "G.se2", "G.se3", "G.sgh2", "G.sgh3", "G.sgh4", "G.sh2", "G.sh3", "G.si2", "G.si3", "G.sl2", "G.sl3", "G.sl4", "c.s1", "c.s2", "c.s3", "c.s4", "c.s5", "c.ss1", "c.ss2", "c.ss3", "c.ss4", "c.ss5", "c.sz1", "c.sz3", "c.sz11", "c.sz13", "c.sz21", "c.sz23", "c.sz31", "c.sz33", "G.xgh2", "G.xgh3", "G.xgh4", " G.xh2", " G.xh3", "G.xi2", "G.xi3", "G.xl2", "G.xl3", "G.xl4", "G.no", "c.z1", "c.z3", "c.z11", "c.z13", "c.z21", "c.z23", "c.z31", "c.z33", "G.zmol", "G.zmos"}

	var G satrec
	var c cof
	G_ans := []float64{}
	G_base := []float64{0.000270526, 0.999999963, 7.64451097683612e-09, 9.13826184204173e-08, 3.15e-05, 9.9225e-10, -2.60329326895906e-07, 5.32643099071963e-07, -0.016516014, -0.007432492, -0.000205811, -0.000705114, -0.001318729, -0.001446983, 0.000744822, 0.016021805, 0.006235498, 0.000480226, -5.18055685687134e-08, -0.0000548207075345159, 0.000109641, 0.000109641, -0.085463061, -3.22540171002146e-07, -0.000341312, 0.000682625, 0.000682625, -0.052287615, 11.50906211, 6.941766362, 0.645351725, -0.445763758, 2.055846393, 0.123995134, 6.421715157, 0.977660989, -0.0000847711905045748, 0.002767642, -0.000108348, -0.000191635, 0.000197672, 0.000213839, 0.000224997, -0.000155882, -0.002597566, 0.000252811, 0.004374995, 2.85847208, 14.704206, 1.230243367, -0.821877203, 0.332495329, 2.13538825, -2.857397595, 9.763939273, 2.201739895, 5.687936935}
	tc := 0.0
	G.nodeo = 6.191936003226819
	G.argpo = 4.442295787980063
	G.inclo = 0.0002705260340591211
	G.ecco = 3.15e-05
	G.epochF = 26997.619421296287
	G.no = 0.004374995068569369

	dscom(&G, &c, tc)

	G_ans = append(G_ans, c.sinim, c.cosim, G.e3, G.ee2, G.ecco, c.emsq, G.se2, G.se3, G.sgh2, G.sgh3, G.sgh4, G.sh2, G.sh3, G.si2, G.si3, G.sl2, G.sl3, G.sl4, c.s1, c.s2, c.s3, c.s4, c.s5, c.ss1, c.ss2, c.ss3, c.ss4, c.ss5, c.sz1, c.sz3, c.sz11, c.sz13, c.sz21, c.sz23, c.sz31, c.sz33, G.xgh2, G.xgh3, G.xgh4, G.xh2, G.xh3, G.xi2, G.xi3, G.xl2, G.xl3, G.xl4, G.no, c.z1, c.z3, c.z11, c.z13, c.z21, c.z23, c.z31, c.z33, G.zmol, G.zmos)

	for i, v := range G_ans {
		e := math.Abs((v - G_base[i]) / v)
		if e >= 0.00001 {
			t.Errorf("GEO value: %v, e: %v, this: %v, baseline: %v", decoder[i], e*10, v, G_base[i])
		}
	}
}

func TestDpper2(t *testing.T) {
	decoder := []string{"ep", "inclp", "nodep", "argpp", "mp"}

	// ================================= GEO TEST CASE 1 - DPPER =================================
	var G satrec
	G_ans1 := []float64{}
	G_base1 := []float64{3.15e-05, 0.0002705260340591211, 6.191936003226819, 4.442295787980063, 1.801858721137174}

	G.t = 0.0
	G.inclo = 0.0002705260340591211
	init := "y"
	G.ecco = 3.15e-05
	G.nodeo = 6.191936003226819
	G.argpo = 4.442295787980063
	G.mo = 1.801858721137174
	opsmode := "i"
	G.e3 = 7.644510976836127e-09
	G.ee2 = 9.138261842041732e-08
	G.peo = 0
	G.pgho = 0
	G.pho = 0
	G.pinco = 0
	G.plo = 0
	G.se2 = -2.603293268959062e-07
	G.se3 = 5.326430990719636e-07
	G.sgh2 = -0.01651601381446915
	G.sgh3 = -0.007432491692699144
	G.sgh4 = -0.0002058113472108935
	G.sh2 = -0.0007051137583276835
	G.sh3 = -0.0013187293886349073
	G.si2 = -0.0014469831366815606
	G.si3 = 0.0007448223807063925
	G.sl2 = 0.01602180512327805
	G.sl3 = 0.00623549779429091
	G.sl4 = 0.0004802264772678869
	G.xgh2 = -8.477119050457478e-05
	G.xgh3 = 0.0027676424659859498
	G.xgh4 = -0.0001083476462637093
	G.xh2 = -0.00019163490436296422
	G.xh3 = 0.0001976717311208955
	G.xi2 = 0.00021383874647060809
	G.xi3 = 0.0002249974031644886
	G.xl2 = -0.0001558820432923891
	G.xl3 = -0.0025975660572708273
	G.xl4 = 0.0002528111748482556
	G.zmol = 2.2017398948110127
	G.zmos = 5.687936935123361

	dpper(&G, init, opsmode, G.ecco, G.inclo, G.nodeo, G.argpo, G.mo)
	G_ans1 = append(G_ans1, G.ecco, G.inclo, G.nodeo, G.argpo, G.mo)

	for i, v := range G_ans1 {
		e := math.Abs((v - G_base1[i]) / v)
		if e >= 0.000001 {
			t.Errorf("GEO1 value: %v, e: %v, this: %v, baseline: %v", decoder[i], e*10, v, G_base1[i])
		}
	}

	// ================================= GEO TEST CASE 2 - DPPER =================================
	G_ans2 := []float64{}
	G_base2 := []float64{3.1598213173349556e-05, 0.00040788478363964895, 6.436349906567948, 4.210453941719098, -1.5551617572068033}

	G.t = 45188.03333334625
	init = "n"
	ep := 3.154078510204699e-05
	inclp := 7.35606791909354e-05
	nodep := 6.1845889392109115
	argpp := 4.458217072328013
	mp := -1.551309399144703

	ep, inclp, nodep, argpp, mp = dpper(&G, init, opsmode, ep, inclp, nodep, argpp, mp)
	G_ans2 = append(G_ans2, ep, inclp, nodep, argpp, mp)

	for i, v := range G_ans2 {
		e := math.Abs((v - G_base2[i]) / v)
		if e >= 0.000001 {
			t.Errorf("GEO2 value: %v, e: %v, this: %v, baseline: %v", decoder[i], e*10, v, G_base2[i])
		}
	}
}

func TestDsinit(t *testing.T) {
	decoder := []string{"em", "argpm", "inclm", "mm", "nm", "nodem", "atime", "d2201", "d2211", "d3210", "d3222", "d4410",
		"d4422", "d5220", "d5232", "d5421", "d5433", "dedt", "didt", "dmdt", "dnodt", "domdt", "del1", "del2", "del3",
		"xfact", "xlamo", "xli", "xni"}

	// ================================= GEO TEST CASE 1 - DPPER =================================
	var G satrec
	var c cof
	G_ans := []float64{}
	G_base := []float64{3.15e-05, 0, 0.0002705260340591211, 0, 0.004374995068569369, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9.025642197376453e-13,
		-4.3587945820787964e-09, -9.814964016877237e-08, 0, 2.715064341617845e-08, -6.397477159775755e-13, 1.4104020776079415e-11,
		1.978470220518558e-12, -0.004375014923455862, 1.0528068360423042, 1.0528068360423042, 0.004374995068569369}

	c.cosim = 0.9999999634078327
	c.emsq = 9.9225e-10
	G.argpo = 4.442295787980063
	c.s1 = -5.180556856871346e-08
	c.s2 = -5.4820707534515914e-05
	c.s3 = 0.00010964141501463598
	c.s4 = 0.00010964141496024014
	c.s5 = -0.085463060852343
	c.sinim = 0.00027052603075940977
	c.ss1 = -3.225401710021465e-07
	c.ss2 = -0.00034131235060548895
	c.ss3 = 0.0006826247008723107
	c.ss4 = 0.0006826247005336434
	c.ss5 = -0.05228761520727895
	c.sz1 = 11.509062114870048
	c.sz3 = 6.941766361930908
	c.sz11 = 0.6453517254794157
	c.sz13 = -0.4457637577210542
	c.sz21 = 2.0558463933535784
	c.sz23 = 0.12399513409577811
	c.sz31 = 6.4217151569697375
	c.sz33 = 0.9776609886333322
	G.t = 0.0
	tc := 0.0
	G.gsto = 5.100098369122165
	G.mo = 1.801858721137174
	G.mdot = 0.004375157637389939
	G.no = 0.004374995068569369
	G.nodeo = 6.191936003226819
	G.nodedot = -1.625887093095899e-07
	G.inclo = 0.0002705260340591211
	G.argpo = 4.44229578798006
	inclm := 0.000270526034059121
	c.xpidot = 1.6259473162131596e-07
	c.z1 = 2.8584720803167136
	c.z3 = 14.704205998125415
	c.z11 = 1.2302433670870876
	c.z13 = -0.8218772027069723
	c.z21 = 0.3324953289171527
	c.z23 = 2.1353882503231607
	c.z31 = -2.8573975953790516
	c.z33 = 9.763939273492966
	G.ecco = 3.15e-05
	c.eccsq = 9.9225e-10

	G.irez = 0
	G.atime = 0

	em, argpm, mm, nm, inclm, nodem := dsinit(&G, &c, tc, inclm)

	G_ans = append(G_ans, em, argpm, inclm, mm, nm, nodem, G.atime, G.d2201, G.d2211, G.d3210, G.d3222, G.d4410, G.d4422, G.d5220, G.d5232,
		G.d5421, G.d5433, G.dedt, G.didt, G.dmdt, G.dnodt, G.domdt, G.del1, G.del2, G.del3, G.xfact, G.xlamo, G.xli, G.xni)

	for i, v := range G_ans {
		e := math.Abs((v - G_base[i]) / v)
		if e >= 0.00001 {
			t.Errorf("GEO value: %v, e: %v, this: %v, baseline: %v", decoder[i], e*10, v, G_base[i])
		}
	}
	if G.irez != 1 {
		t.Errorf("GEO irez is %v, was supposed to be %v.", G.irez, 1)
	}
}

func TestDspace(t *testing.T) {
	decoder := []string{"atime", "em", "argpm", "inclm", "xli", "mm", "xni", "nodem", "dndt", "nm"}

	var G satrec
	G_ans := []float64{}
	G_base := []float64{0, 3.15e-05, 4.442295787980063, 0.0002705260340591211, 1.0528068360423042, -4.481326586042414, 0.004374995068569369, 6.191936003226819, 0, 0.004374995068569369}

	G.irez = 1
	G.d2201 = 0
	G.d2211 = 0
	G.d3210 = 0
	G.d3222 = 0
	G.d4410 = 0
	G.d4422 = 0
	G.d5220 = 0
	G.d5232 = 0
	G.d5421 = 0
	G.d5433 = 0
	G.dedt = 9.025642197376453e-13
	G.del1 = -6.397477159775755e-13
	G.del2 = 1.4104020776079415e-11
	G.del3 = 1.978470220518558e-12
	G.didt = -4.3587945820787964e-09
	G.dmdt = -9.814964016877237e-08
	G.dnodt = 0
	G.domdt = 2.715064341617845e-08
	G.argpo = 4.442295787980063
	G.argpdot = 3.251834409309059e-07
	G.t = 0
	tc := 0.0
	G.gsto = 5.100098369122165
	G.xfact = -0.004375014923455862
	G.xlamo = 1.0528068360423042
	G.no = 0.004374995068569369
	G.atime = 0
	em := 3.15e-05
	argpm := 4.442295787980063
	inclm := 0.0002705260340591211
	G.xli = 1.0528068360423042
	mm := 1.801858721137174
	G.xni = 0.004374995068569369
	nodem := 6.191936003226819
	nm := 0.004374995068569369

	em, argpm, mm, inclm, nodem, nm = dspace(&G, tc, em, argpm, mm, inclm, nodem, nm)
	G_ans = append(G_ans, G.atime, em, argpm, inclm, G.xli, mm, G.xni, nodem, G.dndt, nm)

	for i, v := range G_ans {
		e := math.Abs((v - G_base[i]) / v)
		if e >= 0.0000001 {
			t.Errorf("GEO value: %v, e: %v, this: %v, baseline: %v", decoder[i], e*10, v, G_base[i])
		}
	}

	// TEST CASE 2
	G_ans2 := []float64{}
	G_base2 := []float64{44640, 3.154078510204699e-05, 4.458217072328013, 7.35606791909354e-05, 1.0601459078371092, -7.834494706324289, 0.004375362079161073, 6.1845889392109115, 3.714383965836704e-07, 0.004375366506965953}

	G.t = 45188.03333334625
	tc = 45188.03333334625

	em = 3.15e-05
	argpm = 4.456990188148301 // 4.442295787980063
	inclm = 0.0002705260340591211
	mm = 199.50662787815813    // 1.801858721137174
	nodem = 6.1845889392109115 // 6.191936003226819
	nm = 0.004374995068569369

	em, argpm, mm, inclm, nodem, nm = dspace(&G, tc, em, argpm, mm, inclm, nodem, nm)
	G_ans2 = append(G_ans2, G.atime, em, argpm, inclm, G.xli, mm, G.xni, nodem, G.dndt, nm)

	for i, v := range G_ans2 {
		e := math.Abs((v - G_base2[i]) / v)
		if e >= 0.0000001 {
			t.Errorf("GEO2 value: %v, e: %v, this: %v, baseline: %v", decoder[i], e*10, v, G_base2[i])
		}
	}
}

func TestSgp4(t *testing.T) {
	var G satrec = satrec{
		epoch_jday:    2.4602791194212963e+06,
		ndot:          1.6665470288140301e-12,
		nddot:         0,
		bstar:         0,
		inclo:         0.0002705260340591211,
		nodeo:         6.191936003226819,
		ecco:          3.15e-05,
		argpo:         4.442295787980063,
		mo:            1.801858721137174,
		no:            0.004374995068569369,
		method:        "d",
		operationmode: "i",
		init:          "n",
		gsto:          5.100098369122165,
		isimp:         1,
		con41:         1.9999988083621068,
		cc5:           2.5301454032485523e-11,
		d4:            0,
		argpdot:       3.251834409309059e-07,
		t:             0,
		t4cof:         0,
		x7thm1:        5.999997219511582,
		xlcof:         1.4740812262608144e-06,
		cc1:           0,
		d2:            0,
		delmo:         0.9999744459150449,
		omgcof:        -0,
		t2cof:         0,
		t5cof:         0,
		mdot:          0.004375157637389939,
		xmcof:         0,
		aycof:         7.370406314280264e-07,
		cc4:           2.9404363505917423e-17,
		d3:            0,
		eta:           3.719502845151615e-05,
		sinmao:        0.9734236435400782,
		t3cof:         0,
		x1mth2:        3.972126310886883e-07,
		nodedot:       -1.625887093095899e-07,
		nodecf:        -0,
		irez:          1,
		d3210:         0,
		d4422:         0,
		d5421:         0,
		del1:          -6.397477159775755e-13,
		didt:          -4.3587945820787964e-09,
		domdt:         2.715064341617845e-08,
		peo:           0,
		pinco:         0,
		se3:           5.326430990719636e-07,
		sgh4:          -0.0002058113472108935,
		si2:           -0.0014469831366815606,
		sl3:           0.00623549779429091,
		xfact:         -0.004375014923455862,
		xgh4:          -0.0001083476462637093,
		xi2:           0.00021383874647060809,
		xl3:           -0.0025975660572708273,
		zmol:          2.2017398948110127,
		xli:           1.0528068360423042,
		d2201:         0,
		d3222:         0,
		d5220:         0,
		d5433:         0,
		del2:          1.4104020776079415e-11,
		dmdt:          -9.814964016877237e-08,
		e3:            7.644510976836127e-09,
		pgho:          0,
		plo:           0,
		sgh2:          -0.01651601381446915,
		sh2:           -0.0007051137583276835,
		si3:           0.0007448223807063925,
		sl4:           0.0004802264772678869,
		xgh2:          -8.477119050457478e-05,
		xh2:           -0.00019163490436296422,
		xi3:           0.0002249974031644886,
		xl4:           0.0002528111748482556,
		zmos:          5.687936935123361,
		xni:           0.004374995068569369,
		d2211:         0,
		d4410:         0,
		d5232:         0,
		dedt:          9.025642197376453e-13,
		del3:          1.978470220518558e-12,
		dnodt:         0,
		ee2:           9.138261842041732e-08,
		pho:           0,
		se2:           -2.603293268959062e-07,
		sgh3:          -0.007432491692699144,
		sh3:           -0.0013187293886349073,
		sl2:           0.01602180512327805,
		xgh3:          0.0027676424659859498,
		xh3:           0.0001976717311208955,
		xl2:           -0.0001558820432923891,
		xlamo:         1.0528068360423042,
		atime:         0,
	}

	tprop, err := time.Parse("2006-01-02T15:04:05.99", "2024-01-01T00:00:00.00")
	if err != nil {
		t.Errorf("Error parsing input tprop: %v", err)
	}

	p, v := Propogate(&G, tprop)

	t.Errorf("%+v, %+v", p, v)
}
