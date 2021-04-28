CFLAGS = -Wall -std=gnu11 -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY:default
default: out.txt

#out.mc_comparison.png: out.txt Makefile\
	echo '\
	set terminal png;\
	set output $@;\
	set xlabel "Number of points";\
	set ylabel "Error estimate";\
	set key bottom left;\
	set logscale x;\
	set logscale y;\
	set title "Error estimates from Monte Carlo Integration";\
	plot \
	"$<" using 1:2 with linespoints title "Plain Monte Carlo"\
	,"$<" using 1:3 with linespoints title "Low Discrepancy Monte Carlo"\
	' | pyxplot

out.txt: main
	./$< > $@

main: main.c newton.c runge_kutta.c
