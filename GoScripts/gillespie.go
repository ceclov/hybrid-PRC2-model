package main

import (
	"fmt"
	"math"
)

func GStep(n nucleosomes, tmax float64, nGen int, currTime float64, genStart int, genTime float64) (nucleosomes, float64, float64) {

	Btime := 0.0

	thisTime := 0.0
	currCellCycle := currTime

	adjustGtime := 0.0
	adjustEtime := 0.0

	n.Pme3end = make([]float64, nGen)

	totTime := float64(nGen)*(tCC/tCCwarm)*24 + currTime

	if genTime > 0 {
		adjustGtime = genTime
		adjustEtime = genTime

	}

	for {

		if currTime > totTime {

			currTime = totTime
			genTime = math.Mod((totTime - (tCC - adjustEtime)), tCC)
			n.thisTime = currTime

			break
		}

		n.loopRates()

		n.transcription = n.Rtf()

		thisTime, _ = DirectGillespie(n, false)

		genTime = genTime + thisTime
		currTime = currTime + thisTime

		n.thisTime = currTime

		if genTime > tmax {

			if !noReplication {
				if EqRep {
					n.RepEq()

				} else {
					n.Rep()

				}
			}

			n.nCellCycle++

			genTime = 0.0

			currCellCycle = currCellCycle + tmax - adjustGtime

			adjustGtime = 0

			currTime = currCellCycle

			Btime = genTime - timeRec + float64(n.nCellCycle)
		}

		if genTime > timeRec {

			n.time = append(n.time, currTime)

			Btime = genTime - timeRec + float64(n.nCellCycle)

			n.Btime = append(n.Btime, Btime)

			n.me0 = append(n.me0, meanState(n, 0))

			n.me3 = append(n.me3, meanState(n, 3))

			n.me2 = append(n.me2, meanState(n, 2))

			n.protein = append(n.protein, meanStateP(n, 1))

			n.NR = append(n.NR, meanStateNuc(n, 3, true))
			n.NRme2 = append(n.NRme2, meanStateNuc(n, 2, true))
			n.body = append(n.body, meanStateNuc(n, 3, false))
			n.ON = append(n.ON, ONOFFNuc(n, true))

		}

	}

	if len(n.me3) == 0 {

		currTime = tmax * (float64(nGen))
		n.time = append(n.time, currTime)
		n.me3 = append(n.me3, meanState(n, 3))
		n.me0 = append(n.me0, meanState(n, 0))
		n.protein = append(n.protein, meanStateP(n, 1))
		n.NR = append(n.NR, meanStateNuc(n, 3, true))
		n.NRme2 = append(n.NRme2, meanStateNuc(n, 2, true))
		n.body = append(n.body, meanStateNuc(n, 3, false))
		n.ON = append(n.ON, ONOFFNuc(n, true))

	}

	return n, currTime, genTime

}

func DirectGillespie(n nucleosomes, fast bool) (float64, bool) {
	var time float64
	var d []reaction
	var dr reaction
	var r2 float64

	d = n.Props()
	d = cumReac(d)

	var err bool

	r := rng.Float64()

	if len(d) == 0 {
		err = true
		panic("no reactions!!")

	}

	a0 := d[len(d)-1].rate

	if a0 == 0 {
		fmt.Println(d[0:6])
		panic("a0 is zero")

	}

	time = (1 / a0) * math.Log(1/r)

	r2 = rng.Float64()

	dr = FindReac(d, r2*a0)

	if dr.prop == 0 {
		zeroProp++
		panic("zero propensities")
	}

	r2n(n, dr)

	return (time / 3600.0), err
}

func FindReac(d []reaction, item float64) reaction {
	var dr reaction
	found := false

	for i := 0; i < len(d); i++ {
		if d[i].rate > item || d[i].rate == item {
			dr = d[i]
			found = true
			break
		}

	}
	if !found {
		panic("no reaction found")
	}

	return dr
}
