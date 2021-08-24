package main

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"strconv"
)

func importInputParam(name string, notFile bool, s []string) []string {
	//Not all parameters are used.
	var thisVar string
	var inputs []string

	f, _ := os.Open(name)

	r := csv.NewReader(bufio.NewReader(f))

	for {

		record, err := r.Read()

		if notFile {
			record = s
		}

		if err == io.EOF {
			break
		}

		thisVar = record[0]
		inputs = append(inputs, thisVar)

		switch thisVar {

		case "Lov1":
			b, err := strconv.ParseBool(record[1])
			if err == nil {
				Lov1 = b
			}

		case "PRE":
			b, err := strconv.ParseBool(record[1])
			if err == nil {
				PRE = b
			}

		case "nSim":
			b, err := strconv.ParseInt(record[1], 10, 0)
			if err == nil {
				nSim = int(b)

			}

		case "nGen":
			b, err := strconv.ParseInt(record[1], 10, 0)
			if err == nil {
				nGen = int(b)

			}

		case "postGen":

			b, err := strconv.ParseInt(record[1], 10, 0)
			if err == nil {
				postGen = int(b)

			}

		case "Nnuc":
			b, err := strconv.ParseInt(record[1], 10, 0)
			if err == nil {
				Nnuc = int(b)

			}

		case "N":
			b, err := strconv.ParseInt(record[1], 10, 0)
			if err == nil {
				N = int(b)

			}

		case "Ntot":
			b, err := strconv.ParseInt(record[1], 10, 0)
			if err == nil {
				Ntot = int(b)

			}

		case "NPnuc":
			b, err := strconv.ParseInt(record[1], 10, 0)
			if err == nil {
				NPnuc = int(b)

			}

		case "timeRec":

			b, err := strconv.ParseFloat(record[1], 64)

			if err == nil {
				timeRec = b

			}

		case "tCC":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				tCC = b

			}

		case "ONtres":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				ONtres = b

			}

		case "kme":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				kme = b

			}

		case "lowKme":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				lowKme = b

			}

		case "highKme":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				highKme = b

			}

		case "PREkme":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				PREkme = b

			}

		case "pdem":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				pdem = b

			}

		case "pex":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				pex = b
			}
		case "usePT":
			b, err := strconv.ParseBool(record[1])
			if err == nil {
				usePT = b
			}

		case "kp01":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				kp01 = b

			}

		case "gdep":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				gdep = b

			}

		case "gp01":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				gp01 = b

			}
		case "Switchgp01":

			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				Switchgp01 = b

			}

		case "alpha":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				alpha = b

			}

		case "normVin3":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				normVin3 = b

			}

		case "divider":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				divider = b

			}

		case "vin3Lim":
			b, err := strconv.ParseInt(record[1], 10, 0)
			if err == nil {
				vin3Lim = int(b)

			}

		case "extraGens":
			b, err := strconv.ParseInt(record[1], 10, 0)
			if err == nil {
				extraGens = int(b)

			}

		case "NVnoise":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				NVnoise = b
			}

		case "vernAlpha":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				vernAlpha = b

			}

		case "wzero":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				wzero = b

			}

		case "ProtMe3":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				ProtMe3 = b
			}

		case "NValpha":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				NValpha = b
			}

		case "RepEff":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				RepEff = b
			}

		case "fmax":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				fmax = b
			}

		case "fmin":
			b, err := strconv.ParseFloat(record[1], 64)
			if err == nil {
				fmin = b
			}

		}
	}
	return inputs
}

func PrintParams(params []string, name string) {

	f, err := os.Create(name)

	if err != nil {
		fmt.Println(err)
		fmt.Println(name)
		panic("something wrong with filename")
	}
	defer f.Close()

	var thisVar, t string

	str := "VarName" + "\t" + "VarVariable" + "\n"
	n3, err := f.WriteString(str)

	for i := 0; i < len(params); i++ {

		thisVar = params[i]

		switch thisVar {

		case "Lov1":
			t = strconv.FormatBool(Lov1)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "PRE":
			t = strconv.FormatBool(PRE)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "nSim":
			t = strconv.Itoa(nSim)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "nGen":
			t = strconv.Itoa(nGen)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "postGen":
			t = strconv.Itoa(postGen)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "Nnuc":
			t = strconv.Itoa(Nnuc)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "N":
			t = strconv.Itoa(N)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "Ntot":
			t = strconv.Itoa(Ntot)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "NPnuc":
			t = strconv.Itoa(NPnuc)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "timeRec":
			t := strconv.FormatFloat(timeRec, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "ONtres":
			t := strconv.FormatFloat(ONtres, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "tCC":
			t := strconv.FormatFloat(tCC, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "divider":
			t := strconv.FormatFloat(divider, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "kme":
			t := strconv.FormatFloat(kme, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "highKme":
			t := strconv.FormatFloat(highKme, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "PREkme":
			t := strconv.FormatFloat(PREkme, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "lowKme":
			t := strconv.FormatFloat(lowKme, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "pdem":
			t := strconv.FormatFloat(pdem, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "pex":
			t := strconv.FormatFloat(pex, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "usePT":
			t = strconv.FormatBool(usePT)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "kp01":
			t := strconv.FormatFloat(kp01, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "gdep":
			t := strconv.FormatFloat(gdep, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "gp01":
			t := strconv.FormatFloat(gp01, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "Switchgp01":
			t := strconv.FormatFloat(Switchgp01, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "alpha":
			t := strconv.FormatFloat(alpha, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "normVin3":
			t := strconv.FormatFloat(normVin3, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "vin3Lim":
			t = strconv.Itoa(vin3Lim)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "NVnoise":
			t := strconv.FormatFloat(NVnoise, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "ProtMe3":
			t := strconv.FormatFloat(ProtMe3, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "NValpha":
			t := strconv.FormatFloat(NValpha, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "vernAlpha":
			t := strconv.FormatFloat(vernAlpha, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "gdem":
			t := strconv.FormatFloat(gdem, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "gme23":
			t := strconv.FormatFloat(gme23, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "NRgme23":
			t := strconv.FormatFloat(NRgme23, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "zeroProp":
			t = strconv.Itoa(zeroProp)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "fmax":
			t := strconv.FormatFloat(fmax, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "fmin":
			t := strconv.FormatFloat(fmin, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "RepEff":

			t := strconv.FormatFloat(RepEff, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "postV":
			t = strconv.FormatBool(postV)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "EqRep":
			t = strconv.FormatBool(EqRep)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)
		case "noReplication":
			t = strconv.FormatBool(noReplication)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "PT":
			t := strconv.FormatFloat(PT, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "wzero":
			t := strconv.FormatFloat(wzero, 'f', -1, 64)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		case "extraGens":
			t = strconv.Itoa(extraGens)
			str = thisVar + "\t" + t + "\n"
			n3, err = f.WriteString(str)

		}

		if err != nil {
			fmt.Println(n3)
			panic("error in print params")
		}

	}
	f.Close()

}
