package main

import (
	"fmt"
	"os"
)

func printOutput(outString, transString, name string, params []string) {

	params = append(params, "gdem", "gme23", "NRgme23", "zeroProp", "fmax", "postV", "EqRep", "noReplication", "NRtoB", "BtoNR", "PT", "pex", "ProtMe3", "divider", "wzero")
	PrintParams(params, "../Output/"+name+"_"+extraName+"_NucleationParams.txt")

	resFile, _ := os.Create(("../Output/" + name + "_" + extraName + ".txt"))

	defer resFile.Close()

	resFile.WriteString(outString)

	transFile, _ := os.Create(("../Output/" + name + "_" + extraName + "_Trans.txt"))

	defer transFile.Close()

	transFile.WriteString(transString)

}

func TestLoop() (string, string) {

	ONtres = 1.0

	run_loop = true

	par1 := []float64{2, 4, 6, 8}

	par2 := []float64{0}

	name1 := "nGen"

	name2 := "dummy"

	Lov1 = false
	PRE = false
	lhp1 = false
	noProtein = false
	noWarmRep = false

	timeRec = 0.0
	timeRecWarm = 0.0

	tCC = 22.0

	tCCwarm = 22.0

	fmax = 7.5e-4
	fmin = fmax / 25

	pdem = 0.16595869

	pex = 0.08317638

	NVnoise = 1

	kp01 = 0.05

	gdep = 0.001
	PREkme = 1e-6

	normVin3 = 1.8e-4
	ONtres = 1.0

	vin3Lim = 3
	ProtMe3 = 0.01
	NPnuc = 17
	extraGens = 1

	nSim = 4000

	run_type := 5
	//run_type = 1
	fig_type := "S7"
	wzero = 3

	highKme = 8e-6 * 5

	lowKme = (8e-6 * 4) / 1.9

	extraName = "TestLoop" + fig_type

	var res averages
	var time []float64

	outStr := ""
	transStr := ""

	gens := 8
	nGen = 8

	if fig_type == "2C" {
		par1 = []float64{8}

	}

	if fig_type == "4A" {
		PRE = true
		nGen = 8
		gens = 8
		name1 = "PREkme"
		par1 = []float64{0.0, 1.2 * 8e-6, 1.5 * 8e-6}

	}

	if fig_type == "4B" || fig_type == "S7" {

		Lov1 = true

		nGen = 12
		gens = 12
		par1 = []float64{2, 4, 8, 12}

		name1 = "nGen"

		NValpha = 2.0

		postGen = 40

		normVin3 = 2.1e-4
		vin3Lim = 3

		ProtMe3 = 5e-2

		extraGens = 1
	}

	if fig_type == "S7" {
		noWarmRep = true
	}

	if run_type == 1 {

		nSim = 4000
		divider = 100
		lowKme = 0
		nGen = 20
		gens = 20

		EqRep = false
		noProtein = false

		ONtres = 1.0
		PT = (1.0 / 3.0)

		if fig_type == "1C" {
			noProtein = true
			par1 = []float64{3, 4, 5, 6, 7, 8, 9, 10}
			ONtres = 0.5
			PT = (1.0 / 6.0)
			name1 = "Nnuc"

		}

		if fig_type == "2D" {
			name1 = "NPnuc"
			par1 = []float64{10, 12, 14, 15, 16, 17, 18, 20}

		}

		if fig_type == "1S2A" {
			noProtein = true
			par1 = []float64{3, 4, 5, 6, 7, 8, 9, 10}
			name1 = "Nnuc"

		}

		if fig_type == "1S2B" {
			EqRep = true
			RepEff = 1.0

			ONtres = 0.5
			PT = (1.0 / 6.0)

			noProtein = true
			par1 = []float64{3, 4, 5, 6, 7, 8, 9, 10}
			name1 = "Nnuc"

		}

		if fig_type == "1S2C" {
			name1 = "NPnuc"
			par1 = []float64{5, 6, 7, 8, 9, 10}

			gdep = 0.0008
			kp01 = 0.05

			EqRep = true
			RepEff = 1.0
		}

	}

	fmt.Println(nSim, name1, name2, fig_type)

	for i := 0; i < len(par1); i++ {

		for j := 0; j < len(par2); j++ {
			res, time = runManySims(name1, name2, par1[i], par2[j], gens, run_type)

			outStr = outStr + name1 + ":" + StringFloat64short(par1[i]) + name2 + ":" + StringFloat64short(par2[j]) + "\n"
			outStr = outStr + "time" + "\t" + StringFloat64(time) + "me3" + "\t" + StringFloat64(res.me3) + "NR" + "\t" + StringFloat64(res.NR) + "body" + "\t" + StringFloat64(res.body) + "protein" + "\t" + StringFloat64(res.protein) + "FLCON" + "\t" + StringFloat64(res.ONtot)
			transStr = transStr + name1 + ":" + StringFloat64short(par1[i]) + name2 + ":" + StringFloat64short(par2[j]) + "\n" + "Trans" + "\t" + StringFloat64(res.trans)
		}
	}
	return outStr, transStr
}
