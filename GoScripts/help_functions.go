package main

import (
	"fmt"
	"math/rand"
	"strconv"
)

func initResTot(m int) results {

	var r results

	r.me3 = make([][]float64, m)
	r.me0 = make([][]float64, m)
	r.NR = make([][]float64, m)
	r.protein = make([][]float64, m)
	r.body = make([][]float64, m)

	r.ON = make([][]float64, m)
	r.trans = make([][]float64, m)

	return r

}

func averageResTot(r results, m, n int) averages {

	var a averages

	a.me3 = make([]float64, n)
	a.me0 = make([]float64, n)
	a.body = make([]float64, n)
	a.protein = make([]float64, n)
	a.NR = make([]float64, n)

	a.ONtot = make([]float64, n)

	for i := 0; i < n; i++ {

		a.me3[i] = meanFloat(getColumn(r.me3, m, i))
		a.me0[i] = meanFloat(getColumn(r.me0, m, i))
		a.body[i] = meanFloat(getColumn(r.body, m, i))
		a.protein[i] = meanFloat(getColumn(r.protein, m, i))
		a.NR[i] = meanFloat(getColumn(r.NR, m, i))

		a.ONtot[i] = meanFloat(averageON(getColumn(r.ON, m, i)))

	}
	return a
}

func averageResTrans(r results, m, n int, a averages) averages {
	a.trans = make([]float64, n)

	for i := 0; i < n; i++ {

		a.trans[i] = meanFloat(getColumn(r.trans, m, i))

	}
	return a
}

func averageON(v []float64) []float64 {

	l := len(v)

	ncells := int(float64(l) / 2.0)

	out := make([]float64, ncells)
	ind := 0

	var cell1, cell2 float64

	for i := 0; i < ncells; i++ {

		ind = i * 2

		cell1 = v[ind]
		cell2 = v[ind+1]

		if cell1 == 1.0 || cell2 == 1.0 {
			out[i] = 1.0
		} else {
			out[i] = 0.0
		}

	}

	return out
}

func initSystem(N int, state, pstate, nucState, protState int) nucleosomes {

	var n nucleosomes
	n.values = make([]nucleosome, N)
	n.nucleation = make([]protein, NPnuc)
	switch {

	case (N < 0):
		panic("Length of system to small")
	case N == 0:
		panic("System size is zero")

	case state < 0:
		panic("Init state to low")
	case state > 3:
		panic("init state to high")
	}

	for i := 0; i < N; i++ {
		n.values[i].h1 = state
		n.values[i].h2 = state
		n.values[i].p1 = pstate

		n.values[i].p2 = pstate
		n.values[i].lab = "old"
	}

	for i := 0; i < NPnuc; i++ {
		n.nucleation[i].p = protState

	}

	for i := 0; i < Nnuc; i++ {
		n.values[i].h1 = nucState
		n.values[i].h2 = nucState

	}

	n.transcription = 0.0

	n.me0 = []float64{}

	n.me3 = []float64{}

	n.me2 = []float64{}
	n.ON = []float64{}

	n.me1 = []float64{}

	n.time = []float64{}

	n.Btime = []float64{}
	n.NR = []float64{}
	n.NRme2 = []float64{}
	n.body = []float64{}
	n.nCellCycle = 0

	n.trans = make([]float64, nGen)

	return n

}

func (n nucleosomes) initme3(frac float64, nuc bool) {
	var l, start int

	if nuc {
		start = 0
		l = Nnuc
	} else {
		start = Nnuc
		l = N
	}

	for i := start; i < l; i++ {

		if rand.Float64() < frac {
			n.values[i].h1 = 3
			n.values[i].p1 = 1
		} else {
			n.values[i].h1 = 0
			n.values[i].p1 = 0
		}

		if rand.Float64() < frac {
			n.values[i].h2 = 3
			n.values[i].p2 = 1
		} else {
			n.values[i].h2 = 0
			n.values[i].p2 = 0
		}

	}

}

func (n nucleosomes) Resetme3(nucState, protState int, nuc bool) {
	var l, start int

	if nuc {
		start = 0
		l = Nnuc
	} else {
		start = Nnuc
		l = N
	}

	for i := start; i < l; i++ {

		n.values[i].h1 = nucState
		n.values[i].h2 = nucState

	}

	for i := 0; i < NPnuc; i++ {
		n.nucleation[i].p = protState

	}

}

func resetSystem(n nucleosomes) nucleosomes {

	n.me3 = []float64{}

	n.time = []float64{}
	n.protein = []float64{}

	n.NR = []float64{}
	n.NRme2 = []float64{}
	n.body = []float64{}
	n.ON = []float64{}
	n.trans = make([]float64, nGen)

	n.nCellCycle = 0

	return n
}

func printFloat64(v []float64) {

	for i := 0; i < len(v); i++ {
		fmt.Print(v[i], "\t")
	}
	fmt.Println()
}

func StringFloat64(v []float64) string {
	var s, t string

	for i := 0; i < len(v); i++ {
		t = strconv.FormatFloat(v[i], 'f', -1, 64)
		s = s + t + "\t"
	}
	s = s + "\n"
	return s
}

func StringFloat64short(v float64) string {
	s := strconv.FormatFloat(v, 'f', -1, 64) + "\t"

	return s
}

func meanState(n nucleosomes, state int) float64 {

	L := len(n.values)

	if L == 0 {

		return 0
	}
	sum := 0
	for i := 0; i < L; i++ {
		switch n.values[i].h1 {
		case state:
			sum++
		}

		switch n.values[i].h2 {
		case state:
			sum++
		}

	}

	return float64(sum) / float64(2*L)

}

func meanStateNuc(n nucleosomes, state int, NR bool) float64 {
	var divLen float64

	L := len(n.values)

	start := Nnuc
	if NR {

		start = 0
		L = Nnuc
	}
	if L == 0 {

		return 0
	}
	sum := 0
	for i := start; i < L; i++ {
		switch n.values[i].h1 {
		case state:
			sum++
		}

		switch n.values[i].h2 {
		case state:
			sum++
		}
		divLen++
	}

	if divLen == 0 {
		return 0
	}

	return float64(sum) / (2 * divLen)

}

func ONOFFNuc(n nucleosomes, NR bool) float64 {
	var divLen, thisON float64

	L := len(n.values)

	start := Nnuc
	if NR {

		start = 0
		L = Nnuc
	}
	if L == 0 {

		return 0
	}
	sum := 0.0
	for i := start; i < L; i++ {

		sum = Kdelta(n.values[i].h1, 3) + Kdelta(n.values[i].h2, 3)
		sum = sum + Kdelta(n.values[i].h1, 2) + Kdelta(n.values[i].h2, 2)
		divLen++
	}

	if divLen == 0 {
		return 0
	}

	thisfrac := sum / (2 * divLen)

	if thisfrac > PT {
		thisON = 0.0
	} else {
		thisON = 1.0
	}

	return thisON

}

func countPNucState(n nucleosomes, state int) float64 {
	L := len(n.nucleation)

	if L == 0 {
		return 0
	}

	sum := 0

	for i := 0; i < L; i++ {

		switch n.nucleation[i].p {
		case state:
			sum++
		}
	}
	return float64(sum)
}

func meanStateP(n nucleosomes, state int) float64 {
	L := len(n.nucleation)

	if L == 0 {
		return 0
	}

	sum := 0

	for i := 0; i < L; i++ {

		switch n.nucleation[i].p {
		case state:
			sum++
		}
	}
	return float64(sum) / float64(L)

}

func FLCON(tres float64, input float64) float64 {

	var fracTres, out float64
	fracTres = PT

	if input > fracTres {
		out = 0.0
	} else {
		out = 1.0
	}
	return (out)
}

func FLCONMany(tres float64, input1, input2 []float64) []float64 {
	res := make([]float64, len(input1)+1)
	var input float64

	for i := 0; i < len(input1); i++ {
		input = input1[i] + input2[i]
		res[i] = FLCON(tres, input)
	}

	return (res)
}
