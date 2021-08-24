package main

import (
	"math"
	"strconv"
)

func runManySims(name, name2 string, param, param2 float64, maxGen, runType int) (averages, []float64) {

	step := 1.0

	SetParam()

	currParam := strconv.FormatFloat(param, 'f', -1, 64)
	importInputParam("../Params/InputParams.csv", true, []string{name, currParam})
	currParam = strconv.FormatFloat(param2, 'f', -1, 64)
	importInputParam("../Params/InputParams.csv", true, []string{name2, currParam})

	var n3 nucleosomes

	maxGen = nGen
	postGen = int(math.Ceil(40.0 / float64(nGen)))

	maxJ := nGen * int(tCC)

	tot := initResTot(nSim)

	for i := 0; i < nSim; i++ {

		if runType == 1 {
			maxJ = nGen * 1 * int(24)
			n3 = lhp1_1sim()
		}

		if runType == 5 {

			maxJ = nGen*1*int(24*7) + nGen*1*int(24)*postGen
			n3 = test_loop_sim()

		}

		tot.me3[i] = make([]float64, maxJ+1)
		tot.me0[i] = make([]float64, maxJ+1)
		tot.NR[i] = make([]float64, maxJ+1)
		tot.body[i] = make([]float64, maxJ+1)
		tot.protein[i] = make([]float64, maxJ+1)
		tot.ON[i] = make([]float64, maxJ+1)

		n3 = WriteTimeAll(n3, step, maxJ)

		tot.trans[i] = make([]float64, len(n3.trans))

		tot.me3[i] = n3.me3
		tot.me0[i] = n3.me0
		tot.NR[i] = n3.NR
		tot.body[i] = n3.body
		tot.protein[i] = n3.protein
		tot.ON[i] = FLCONMany(ONtres, n3.NR, n3.NRme2)

		tot.trans[i] = n3.trans
	}
	a := averageResTot(tot, nSim, maxJ+1)
	a = averageResTrans(tot, nSim, len(n3.trans), a)

	return a, n3.time
}

func lhp1_1sim() nucleosomes {

	noReplication = false
	alpha = vernAlpha
	var startTime, genTime float64
	startTime = 0
	timeRec = 0.0
	timeRecWarm = 0.0
	genTime = 0

	preV = false
	postV = true
	n3 := nucleosomes{}

	tCC = tCCwarm
	timeRec = timeRecWarm

	kme = lowKme
	NRkme = highKme

	gp01 = 0

	SetParam()

	if noProtein {
		n3 = initSystem(N, 0, 0, 3, 0)
		n3.Resetme3(3, 0, true)

	} else {
		n3 = initSystem(N, 0, 0, 3, 1)
		n3.Resetme3(3, 1, true)
	}

	currLen := math.Ceil(float64(nGen+extraGens) * (24.0 / tCCwarm))

	n3.trans = make([]float64, int(currLen+10))

	genTime = tCCwarm * rng.Float64()

	startTime = 0

	noReplication = true
	n3, startTime, genTime = GStep(n3, tCC, extraGens, startTime, 0, genTime)
	alpha = NValpha
	noReplication = false

	n3, startTime, genTime = GStep(n3, tCC, nGen, startTime, extraGens, genTime)

	return n3
}

func test_loop_sim() nucleosomes {

	alpha = NValpha

	preGens := 10

	ccDiff := 7.0
	noReplication = false
	var startTime, genTime float64
	startTime = 0
	timeRec = 0.0
	timeRecWarm = 0.0
	genTime = 0
	genFirst := 0.0

	preV = true
	postV = false
	n3 := nucleosomes{}

	tCC = tCCwarm
	timeRec = timeRecWarm

	kme = highKme
	NRkme = highKme

	if lhp1 {
		kme = lowKme
	}

	gp01 = 0

	SetParam()

	n3 = initSystem(N, 0, 0, 0, 0)

	if !Lov1 {
		n3.initme3(0.0, true)
	}

	n3.trans = make([]float64, preGens+2)

	startTime = 0.0

	genFirst = rng.Float64() * tCCwarm

	n3, startTime, genTime = GStep(n3, tCC, preGens, startTime, 0, genFirst)

	genTime = genTime * ccDiff
	startTime = startTime - (float64(preGens) * 24)

	n3 = resetSystem(n3)

	preV = false

	alpha = vernAlpha

	tCC = tCCwarm * ccDiff
	timeRec = timeRecWarm * ccDiff

	NRkme = highKme
	kme = lowKme

	if lhp1 {
		kme = lowKme
	}

	SetParam()

	currLen := math.Ceil((float64(nGen*(postGen+2)) * (24.0 / tCCwarm)))

	n3.trans = make([]float64, int(currLen+10))

	n3, startTime, genTime = GStep(n3, tCC, nGen, startTime, 0, genTime)
	genTime = genTime / ccDiff

	tCC = tCCwarm
	timeRec = timeRecWarm

	postV = true

	kme = highKme
	NRkme = highKme

	if lhp1 {
		kme = lowKme
	}

	SetParam()

	noReplication = true
	n3, startTime, _ = GStep(n3, tCC, extraGens, startTime, nGen, 0)
	noReplication = false

	if noWarmRep {
		noReplication = true
	}

	alpha = NValpha

	if PRE {

		n3, startTime, genTime = GStep(n3, tCC, 11, startTime, nGen+extraGens, genTime)
		kme = PREkme
		NRkme = 0.0

		SetParam()

		n3.Resetme3(0, 0, true)
		n3, startTime, genTime = GStep(n3, tCC, (nGen * postGen), startTime, nGen+11, genTime)

	} else {

		n3, startTime, genTime = GStep(n3, tCC, (nGen * postGen), startTime, nGen+extraGens, genTime)

	}

	return n3

}
