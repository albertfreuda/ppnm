CFLAGS = -Wall -std=gnu11 -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY:default
default: out.higgs.png

out.higgs.png: out.txt Makefile
	echo '\
	set terminal png;\
	set output $@;\
	set xlabel "$$E_{cm} [GeV]$$";\
	set ylabel "$$\sigma$$";\
	set key outside;\
	set title "Fit to CERN data on Higgs boson";\
	plot \
	"$<" index 1 using 1:2:5 with yErrorBars title "Data"\
	,"$<" index 1 using 1:3 with line title "Fit"\
	,"$<" index 1 using 1:4 with line title "Initial parameter guess"\
	' | pyxplot

out.txt: main
	./$< > $@

main: main.c newton_min.c simplex.c

simplex_test: simplex_test.c simplex.c

.PHONEY:clean
clean:
	$(RM) main out*
