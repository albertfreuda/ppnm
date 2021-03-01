CC = gcc
CFLAGS = -Wall -std=gnu11 -pthread 
CFLAGS += -fopenmp
LDLIBS = -lm -lgomp

.PHONEY: default
default: out.error.png

out.error.png:out.data.txt Makefile
	echo '\
	set term png ;\
	set out "$@" ;\
	set title "Monte Carlo calculation of $$\pi$$" ;\
	set logscale x ;\
	set logscale y ;\
	set xlabel "Number of points $$N$$" ;\
	set ylabel "Absolute Error" ;\
	a=1 ;\
	f(x) = a/sqrt(x) ;\
	fit f(x) "$<" using 1:3 via a ;\
	plot [][1e-5:8e-3]\
	"$<" using 1:3 with points title "Error"\
	,f(x) with lines title "$$1/\sqrt{N}$$"\
	'|pyxplot

ns = 1 1.5 2 4 8
out.data.txt: main Makefile
	>$@
	for n in $(ns); do ./$< "$$n"e5 >> $@; done
	for n in $(ns); do ./$< "$$n"e6 >> $@; done
	for n in $(ns); do ./$< "$$n"e7 >> $@; done

.PHONEY:clean
clean: 
	$(RM) out* main mainB *.png
