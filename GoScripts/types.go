package main

type nucleosome struct {
	h1  int `validate:"nonzero"`
	h2  int `validate:"number,min=2,max=3"`
	p1  int
	p2  int
	me1 float64
	me2 float64

	de1 float64
	de2 float64

	P1   float64
	P2   float64
	deP1 float64
	deP2 float64
	lab  string
}

type protein struct {
	p   int
	add float64
	deP float64
}

type reaction struct {
	rate       float64
	prop       float64
	part       string
	nucleosome int
	histone    string
	action     string
	protein    string
	nucprot    int
}

type nucleosomes struct {
	nucleation    []protein
	values        []nucleosome
	transcription float64
	time          []float64
	Btime         []float64
	me3           []float64
	body          []float64
	NR            []float64
	NRme2         []float64

	me2 []float64

	me1 []float64

	me0 []float64

	Pme3end []float64

	nCellCycle int

	protein []float64

	trans []float64

	ON       []float64
	thisTime float64
}

type rates interface {
	Rates()
	Rtf() float64
	Pn2r() float64
}

type conversion interface {
	n2r()
	//n2rNuc()
}

type transcription interface {
	DeM()
	Ex()
}

type results struct {
	me3     [][]float64
	me0     [][]float64
	protein [][]float64
	NR      [][]float64
	body    [][]float64
	ON      [][]float64
	trans   [][]float64
}

type averages struct {
	me3     []float64
	me0     []float64
	time    []float64
	protein []float64
	ON      float64
	ONtot   []float64
	NR      []float64
	body    []float64

	trans []float64
}
