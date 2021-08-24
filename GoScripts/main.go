package main

import (
	"fmt"
	"math/rand"
	"os"
	"time"

	"github.com/seehuhn/mt19937"
)

//env GOOS=linux GOARCH=amd64 go build package-import-path

var pdem float64
var fmax float64
var fmin float64
var pex float64
var gdem float64
var PT = 1.0 / 3.0
var alpha float64
var vernAlpha float64
var NValpha float64
var usePT bool
var transMax = 1 / 60.0

var ONtres float64

const beta = 1

const pme2 = 0.1

var kme float64
var NRkme float64
var divider float64

var kme01 float64
var kme12 float64
var gme01 float64
var gme12 float64
var gme23 float64

var NRkme01 float64
var NRkme12 float64
var NRgme01 float64
var NRgme12 float64
var NRgme23 float64

var NVnoise float64
var highKme float64
var lowKme float64
var PREkme float64

var N int

var Ntot int
var Nnuc int
var NPnuc int

var tCC float64
var tCCwarm float64
var timeRec float64
var timeRecWarm float64
var nGen int
var postGen int
var nSim int

var count0 = false
var count2 = true
var count3 = true
var count1 = false
var printBtime = false

var rng = rand.New(mt19937.New())

var gdep float64
var gp01 float64
var kp01 float64
var ProtMe3 float64
var wzero float64 = 3

var PRE = false
var Lov1 = false

var vin3Lim int
var Switchgp01 float64
var normVin3 float64

var extraName = ""

var cold = true
var postV = false
var noProtein = false

var preV = false
var lhp1 = false
var noReplication = false
var Switch = false
var extraGens int

var run_loop = false

var printDummy = true
var zeroProp = 0
var EqRep = false
var noWarmRep = false

var stats []float64

var RepEff = 1.0

func main() {

	if len(os.Args) < 2 {
		fmt.Println("Missing parameter, provide file name!")
		return
	}
	data := os.Args[1]

	outName := string(data)

	now := time.Now()
	params := importInputParam("../Params/InputParams.csv", false, []string{})

	SetParam()

	outString, outTrans := TestLoop()

	printOutput(outString, outTrans, outName, params)

	end := time.Now()

	fmt.Println(end.Sub(now))
}

func SetParam() {

	kme01 = 9.0 * kme
	kme12 = 6.0 * kme
	gme01 = kme01 / divider
	gme12 = kme12 / divider
	gme23 = kme / divider

	NRkme01 = 9.0 * NRkme
	NRkme12 = 6.0 * NRkme
	NRgme01 = NRkme01 / divider
	NRgme12 = NRkme12 / divider
	NRgme23 = NRkme / divider

	gdem = fmin * pdem
	usePT = true

}
