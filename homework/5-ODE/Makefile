CFLAGS = -Wall -std=gnu11 -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY:default
default: out.threebody.png out.SIR.png out.harm_osc.png 

out.threebody.png: out.txt Makefile	
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Stable solution to the 3-body problem";\
	plot \
	"$<" index 2 using 2:4 with points title "Body 1"\
	,"$<" index 2 using 6:8 with points title "Body 2"\
	,"$<" index 2 using 10:12 with points title "Body 3"\
	'| pyxplot

out.SIR.png: out.txt Makefile	
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Solution to the SIR model";\
	plot \
	"$<" index 1 using 1:2 with linespoints title "Suscpetible"\
	,"$<" index 1 using 1:3 with linespoints title "Infected"\
	,"$<" index 1 using 1:4 with linespoints title "Removed"\
	'| pyxplot

out.harm_osc.png:out.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Solution to the harmonic oscillator";\
	plot \
	"$<" index 0 using 1:2 with linespoints title "ODE solution with adaptive step-size"\
	'| pyxplot

out.txt: main
	./$< > $@

main: main.c runge_kutta.c

.PHONEY:clean
clean:
	$(RM) out* main 
