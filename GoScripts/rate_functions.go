package main

import (
	"bufio"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"io"
	"io/ioutil"
	"math"
	"os"
	"strconv"
	"strings"
)

func Kdelta(x int, y int) float64 {
	switch x {
	case y:
		return 1
	}
	return 0

}

func antiKdelta(x float64, y float64) float64 {
	switch x {
	case y:
		return 0
	}
	return 1

}

func (n nucleosomes) loopRates() {

	var currNuc nucleosome
	var ind int
	var currEi []float64
	var currWeeks float64

	nuc := n.values
	prot := n.nucleation

	currPfrac := 0.0
	NumP := countPNucState(n, 1)

	currWeeks = n.thisTime / (7 * 24)

	if NumP > float64(vin3Lim) {

		currPfrac = 1

		gp01 = (currWeeks * normVin3) / (1 + (currWeeks / wzero))

	} else {

		gp01 = (currWeeks * normVin3) / (1 + (currWeeks / wzero))
	}

	if preV || postV || noProtein {
		gp01 = 0.0

	}

	SetParam()

	for i := 0; i < Nnuc; i++ {

		ind = i

		currNuc = nuc[ind]

		currEi = AllEi(nuc, i)

		nuc[ind].me1 = NRHelpRm(currNuc.h1, currEi[0], currPfrac)
		nuc[ind].me2 = NRHelpRm(currNuc.h2, currEi[1], currPfrac)
		nuc[ind].de1 = helpDm(currNuc.h1)
		nuc[ind].de2 = helpDm(currNuc.h2)

	}

	for i := Nnuc; i < len(nuc); i++ {

		ind = i

		currEi = AllEi(nuc, ind)

		currNuc = nuc[ind]

		nuc[ind].me1 = BodyHelpRm(currNuc.h1, currEi[0])
		nuc[ind].me2 = BodyHelpRm(currNuc.h2, currEi[1])
		nuc[ind].de1 = helpDm(currNuc.h1)
		nuc[ind].de2 = helpDm(currNuc.h2)

	}

	for i := 0; i < NPnuc; i++ {

		prot[i].add = (gp01 + kp01*currPfrac) * antiKdelta(float64(NPnuc), NumP)
		prot[i].deP = gdep * (Kdelta(prot[i].p, 1))
	}

}

func BodyHelpRm(h int, currEi float64) float64 {

	res := Kdelta(h, 0) * (gme01 + kme01*currEi)
	res = res + Kdelta(h, 1)*(gme12+kme12*currEi)
	res = res + Kdelta(h, 2)*(gme23+kme*currEi)

	return beta * res
}

func NRHelpRm(h int, currEi, currPfrac float64) float64 {
	res := Kdelta(h, 0) * (NRgme01*NVnoise + NRkme01*currEi)
	res = res + Kdelta(h, 1)*(NRgme12*NVnoise+NRkme12*currEi)
	res = res + Kdelta(h, 2)*(NRgme23*NVnoise+NRkme*currEi)

	res = res + currPfrac*ProtMe3*(Kdelta(h, 0)+Kdelta(h, 1)+Kdelta(h, 2))

	return beta * res
}

func helpDm(h int) float64 {

	return gdem * (Kdelta(h, 1) + Kdelta(h, 2) + Kdelta(h, 3))
}

func (n nucleosomes) Rtf() float64 {

	nuc := n.values
	L := len(nuc)

	bodySum := 0.0
	NRSum := 0.0
	sum := float64(0)
	var res float64

	for i := 0; i < Nnuc; i++ {
		NRSum = NRSum + Kdelta(n.values[i].h1, 2) + Kdelta(n.values[i].h2, 2) + Kdelta(n.values[i].h1, 3) + Kdelta(n.values[i].h2, 3) //how should we do it here?

	}

	for i := Nnuc; i < L; i++ {
		bodySum = bodySum + Kdelta(n.values[i].h1, 2) + Kdelta(n.values[i].h2, 2) + Kdelta(n.values[i].h1, 3) + Kdelta(n.values[i].h2, 3)

	}

	if run_loop {

		if NRSum > bodySum {
			sum = NRSum / float64(Nnuc*2)

		} else if NRSum < bodySum {
			bodyL := L - Nnuc
			sum = bodySum / float64(bodyL*2)

		} else {
			bodyL := L - Nnuc
			sum = bodySum / float64(bodyL*2)

		}
	} else {
		sum = (NRSum + bodySum) / float64(L*2)
	}

	if usePT {
		switch {
		case sum < PT:
			res = alpha * (fmax - (sum/PT)*(fmax-fmin))
		case sum > PT:
			res = alpha * fmin
		case sum == PT:
			res = alpha * fmin
		}
	} else {

		res = alpha * (fmax - sum*(fmax-fmin))
	}

	switch {
	case res > transMax:
		res = transMax

	case res == transMax:
		res = transMax
	}

	return res

}

func (n nucleosomes) DeM() {

	for i := 0; i < len(n.values); i++ {
		switch {
		case rng.Float64() < pdem: //check that this is correct
			n.values[i].h1 = deMeth(n.values[i].h1)

		}

		switch {
		case rng.Float64() < pdem:
			n.values[i].h2 = deMeth(n.values[i].h2)

		}

	}

}

func (n nucleosomes) Ex() {

	nuc := n.values

	for i := 0; i < len(nuc); i++ {
		switch {
		case rng.Float64() < pex: //check that this is correct
			nuc[i].h1 = 0
			nuc[i].h2 = 0

		}

		switch {
		case rng.Float64() < pex:
			nuc[i].h1 = 0
			nuc[i].h2 = 0

		}

	}

}

func (n nucleosomes) Rep() {

	nuc := n.values

	for i := 0; i < len(nuc); i++ {
		switch {

		case rng.Float64() < 0.5:
			nuc[i].h1 = 0
			nuc[i].h2 = 0

		}

	}

	for i := 0; i < NPnuc; i++ {
		if rng.Float64() < 0.5 {
			n.nucleation[i].p = 0
		}
	}

	NumP := countPNucState(n, 1)

	for i := 0; i < int(NumP); i++ {

		n.nucleation[i].p = 1
	}

	for i := int(NumP); i < NPnuc; i++ {
		n.nucleation[i].p = 0
	}

}

func (n nucleosomes) RepEq() {
	var start int
	nuc := n.values

	choice := rng.Intn(2)

	switch choice {

	case 0:
		start = 0

	case 1:
		start = 1

	}

	l := len(n.values)

	for {
		if start > (l - 1) {

			break
		}

		nuc[start].h1 = 0
		nuc[start].h2 = 0

		if rng.Float64() < (1-RepEff) && !(start == l) {

			nuc[start+1].h1 = 0
			nuc[start+1].h2 = 0
		}

		start = start + 2

	}

}

func (n nucleosomes) n2r() []reaction {

	nCol := 4
	numR := 1 + nCol*len(n.values)
	d := make([]reaction, numR)
	var ind int
	nuc := n.values
	var currRate float64

	for i := 0; i < len(nuc); i++ {
		ind = (i * nCol)

		currRate = (-1 / nuc[i].me1) * math.Log(rng.Float64())
		d[ind].rate = currRate
		d[ind].action = "meth"
		d[ind].nucleosome = i
		d[ind].histone = "h1"

		currRate = (-1 / nuc[i].me2) * math.Log(rng.Float64())
		d[ind+1].rate = currRate
		d[ind+1].action = "meth"
		d[ind+1].nucleosome = i
		d[ind+1].histone = "h2"

		currRate = (-1 / nuc[i].de1) * math.Log(rng.Float64())
		d[ind+2].rate = currRate
		d[ind+2].action = "demeth"
		d[ind+2].nucleosome = i
		d[ind+2].histone = "h1"

		currRate = (-1 / nuc[i].de2) * math.Log(rng.Float64())
		d[ind+3].rate = currRate
		d[ind+3].action = "demeth"
		d[ind+3].nucleosome = i
		d[ind+3].histone = "h2"

	}
	currRate = (-1 / n.transcription) * math.Log(rng.Float64())
	d[numR-1].rate = currRate
	d[numR-1].action = "transcription"

	return d
}

func (n nucleosomes) Props() []reaction {

	nCol := 4

	l := len(n.values)

	numR := 1 + nCol*l

	d := make([]reaction, numR)
	dp := make([]reaction, (NPnuc + 1))
	var ind int
	nuc := n.values

	prot := n.nucleation

	for i := 0; i < l; i++ {
		ind = (i * nCol)

		d[ind].prop = nuc[i].me1
		d[ind].action = "meth"
		d[ind].nucleosome = i
		d[ind].histone = "h1"

		d[ind+1].prop = nuc[i].me2
		d[ind+1].action = "meth"
		d[ind+1].nucleosome = i
		d[ind+1].histone = "h2"

		d[ind+2].prop = nuc[i].de1
		d[ind+2].action = "demeth"
		d[ind+2].nucleosome = i
		d[ind+2].histone = "h1"

		d[ind+3].prop = nuc[i].de2
		d[ind+3].action = "demeth"
		d[ind+3].nucleosome = i
		d[ind+3].histone = "h2"

	}

	d[numR-1].prop = n.transcription
	d[numR-1].action = "transcription"

	ind = 0

	for i := 0; i < NPnuc; i++ {
		ind = i

		dp[ind].prop = prot[i].deP
		dp[ind].action = "dePn"
		dp[ind].nucprot = i

	}

	dp[NPnuc].prop = prot[0].add
	dp[NPnuc].action = "add"
	dp[NPnuc].nucprot = ind

	return append(d, dp...)
}

func cumReac(d []reaction) []reaction {
	var res []reaction
	currRate := 0.0

	for i := 0; i < len(d); i++ {

		currRate = currRate + d[i].prop

		d[i].rate = currRate
		res = append(res, d[i])

	}
	return res
}

func r2n(n nucleosomes, r reaction) {
	currAction := r.action

	NumP := countPNucState(n, 1)

	switch currAction {
	case "meth":
		if r.histone == "h1" {
			n.values[r.nucleosome].h1 = addMeth(n.values[r.nucleosome].h1)
		} else {
			n.values[r.nucleosome].h2 = addMeth(n.values[r.nucleosome].h2)
		}

	case "demeth":
		if r.histone == "h1" {
			n.values[r.nucleosome].h1 = deMeth(n.values[r.nucleosome].h1)
		} else {
			n.values[r.nucleosome].h2 = deMeth(n.values[r.nucleosome].h2)
		}

	case "transcription":
		n.DeM()
		n.Ex()

		n.trans[n.nCellCycle]++

	case "P":
		if r.protein == "p2" {
			n.values[r.nucleosome].p2 = addP(n.values[r.nucleosome].p2)
		} else {
			n.values[r.nucleosome].p1 = addP(n.values[r.nucleosome].p1)
		}

	case "deP":
		if r.protein == "p2" {
			n.values[r.nucleosome].p2 = deP(n.values[r.nucleosome].p2)
		} else {
			n.values[r.nucleosome].p1 = deP(n.values[r.nucleosome].p1)
		}

	case "add":

		n.nucleation = addtoOligo(n.nucleation, int(NumP)+1)

	case "dePn":

		n.nucleation[r.nucprot].p = deP(n.nucleation[r.nucprot].p)
	}
}

func deP(value int) int {

	switch value {
	case 1:
		return value - 1

	}

	return value
}

func addP(value int) int {

	switch value {
	case 0:
		return value + 1
	}

	return value
}

func addtoOligo(v []protein, l int) []protein {
	L := len(v)
	for i := 0; i < l; i++ {
		v[i].p = 1
	}

	for i := l; i < L; i++ {
		v[i].p = 0
	}
	return v
}

func deMeth(value int) int {

	switch value {
	case 1, 2, 3:
		return value - 1

	}
	return value
}

func addMeth(value int) int {

	switch value {
	case 0, 1, 2:
		return value + 1

	}
	return value
}

func importParam(name string) []float64 {

	var input []float64
	raw, err := ioutil.ReadFile(name)
	if err != nil {
		fmt.Println(err.Error())

	}
	json.Unmarshal(raw, &input)

	return input
}

func BM(n nucleosomes) float64 {

	sumOFF := 0.0
	sumON := 0.0
	var ON, OFF []float64
	var res float64

	ON = ONOFF(n.me2, n.me3, true)
	OFF = ONOFF(n.me2, n.me3, false)

	sumOFF = timeAverage(n.Btime, OFF)
	sumON = timeAverage(n.Btime, ON)

	res = 4 * sumOFF * sumON

	return res

}

func ONOFF(me2, me3 []float64, ON bool) []float64 {
	l := len(me2)
	res := make([]float64, l)
	var sum float64

	for i := 0; i < l; i++ {
		sum = me2[i] + me3[i]

		res[i] = 0

		if ON {
			if sum < 0.25 {
				res[i] = 1.0
			}
		} else {
			if sum > 0.75 {
				res[i] = 1.0
			}
		}

	}

	return res
}

func timeAverage(time, states []float64) float64 {

	var timeDiff float64
	K := len(states)
	res := 0.0
	nom := time[len(time)-1] - time[0]

	for i := 0; i < (K - 1); i++ {
		timeDiff = time[i+1] - time[i]
		switch {
		case timeDiff < 0:
			fmt.Println(timeDiff, time[i], time[i+1])
			panic("In timeAverage timeDiff is less than 0.0")

		}

		res = (states[i] * (timeDiff) / nom) + res
	}

	return res
}

func meanFloat(v []float64) float64 {

	l := len(v)
	sum := 0.0

	if l == 0 {
		panic("vector zero")
	}

	for i := 0; i < l; i++ {
		sum += v[i]
	}

	return sum / float64(l)
}

func WriteTimeAll(n nucleosomes, step float64, l int) nucleosomes {

	tout := 0.0

	currT := 0.0

	output := nucleosomes{}

	output.time = make([]float64, l+1)

	output.me3 = make([]float64, l+1)

	output.me0 = make([]float64, l+1)

	output.protein = make([]float64, l+1)
	output.NR = make([]float64, l+1)
	output.NRme2 = make([]float64, l+1)
	output.body = make([]float64, l+1)
	output.ON = make([]float64, l+1)
	output.trans = make([]float64, len(n.trans))
	output.trans = n.trans

	L := len(n.time)

	if len(n.me3) == 0 || len(n.protein) == 0 {

		panic("length of me3 or protein is zero")
	}

	output.time[0] = tout

	output.me3[0] = n.me3[0]

	output.ON[0] = n.ON[0]

	output.me0[0] = n.me0[0]

	output.protein[0] = n.protein[0]

	output.body[0] = n.body[0]

	output.NR[0] = n.NR[0]
	output.NRme2[0] = n.NRme2[0]

	j := 0

	maxJ := float64(l)

	if maxJ-n.time[len(n.time)-1] > 1 {
		fmt.Println("curr MAXJ", maxJ)
		fmt.Println("currMax time", n.time[len(n.time)-1])
		panic("maxJ does not correspond")
	}

	for i := 1; i < L; i++ {

		for {

			if tout > n.time[i] || tout > maxJ {

				break

			} else {

				currT = tout

				output.me3[j] = n.me3[i]

				output.me0[j] = n.me0[i]
				output.ON[j] = n.ON[i]

				output.protein[j] = n.protein[i]

				output.time[j] = currT

				output.body[j] = n.body[i]

				output.NR[j] = n.NR[i]
				output.NRme2[j] = n.NRme2[i]

				tout += step

				j++

			}

		}

	}

	if !(currT > maxJ || currT == maxJ) {

		for {

			if tout > maxJ {

				break

			}

			currT = tout

			output.time[j] = currT

			output.me3[j] = n.me3[L-1]
			output.ON[j] = n.ON[L-1]

			output.me0[j] = n.me0[L-1]
			output.body[j] = n.body[L-1]

			output.NR[j] = n.NR[L-1]
			output.NRme2[j] = n.NRme2[L-1]
			output.protein[j] = n.protein[L-1]

			j++

			tout += step

		}

	}

	return output

}

func initRes(m int) results {

	var r results

	r.me3 = make([][]float64, m)
	r.protein = make([][]float64, m)

	return r

}

func getColumn(mat [][]float64, m, n int) []float64 {

	res := make([]float64, m)

	for k := 0; k < m; k++ {

		res[k] = mat[k][n]

	}

	return res

}

func ImportParamCsv(name string) []float64 {
	f, _ := os.Open(name)

	var input []float64
	r := csv.NewReader(bufio.NewReader(f))

	for {
		record, err := r.Read()

		if err == io.EOF {
			break
		}
		v, err := strconv.ParseFloat(strings.Join(record, ""), 64)
		input = append(input, float64(v))

	}
	return input
}

func printStrings(S []string, name string) {
	f, _ := os.Create(name)
	defer f.Close()

	for i := 0; i < len(S); i++ {

		f.WriteString(S[i])

	}

	f.Close()

}

func getParam(start, end float64, l int, log bool) []float64 {
	res := make([]float64, l)
	step := (end - start) / float64(l)

	if log {
		if start == 0 || start < 0 {
			panic("start is equal or less than zero")

		}
		step = (math.Log(end) - math.Log(start)) / float64(l)
	}
	currRes := start
	res[0] = currRes
	for i := 1; i < l; i++ {
		if log {
			currRes = currRes * math.Exp(step)
		} else {
			currRes = currRes + step
		}
		res[i] = currRes
	}

	return res
}

func AllEihelp(item int) float64 {
	return pme2*(Kdelta(item, 2)) + (Kdelta(item, 3))
}

func AllEi(n []nucleosome, ind int) []float64 {

	var neighbors []int
	Ein := make([]float64, 2)
	sum := float64(0)
	L := len(n)

	switch ind {

	case (L - 1):
		neighbors = append(neighbors, ind-1)

	case 0:
		neighbors = append(neighbors, ind+1)
	default:
		neighbors = append(neighbors, ind-1, ind+1)
	}

	for i := 0; i < len(neighbors); i++ {
		sum = sum + AllEihelp(n[neighbors[i]].h1) + AllEihelp(n[neighbors[i]].h2)
	}

	Ein[0] = sum + AllEihelp(n[ind].h2)
	Ein[1] = sum + AllEihelp(n[ind].h1)

	return Ein
}
