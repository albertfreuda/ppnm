CFLAGS = -O1 -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY:default
defaul:plot.cinterp.png plot.cinterp_integ.png plot.cinterp_deriv.png

plot.cinterp.png: out.cinterp.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key bottom left;\
	set tics out;\
	set xlabel "$$x$$";\
	set ylabel "$$y$$";\
	set title "Interpolation";\
	plot \
	"$<" index 0 using 1:2 with line title "Cubic interpolation"\
	,"$<" index 0 using 1:5 with line title "Real function $$y=e^{-x^2}$$"\
	,"$<" index 1 using 1:2 with line title "Cubic interpolation from GSL"\
	,"$<" index 2 using 1:2 with points title "Data points"\
	'| pyxplot

plot.cinterp_integ.png: out.cinterp.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key bottom left;\
	set tics out;\
	set xlabel "$$x$$";\
	set ylabel "$$y$$";\
	set title "Interpolation (integrals)";\
	plot \
	"$<" index 0 using 1:3 with line title "Integral of interpolation"\
	,"$<" index 0 using 1:6 with line title "Real integral (erf from math.h)"\
	,"$<" index 1 using 1:3 with line title "Cubic interpolation integral from GSL"\
	'| pyxplot

plot.cinterp_deriv.png: out.cinterp.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key bottom left;\
	set tics out;\
	set xlabel "$$x$$";\
	set ylabel "$$y$$";\
	set title "Interpolation (derivatives)";\
	plot \
	"$<" index 0 using 1:4 with line title "Derivative of interpolation"\
	,"$<" index 0 using 1:7 with line title "Real derivative $$y=-2xe^{-x^2}$$"\
	,"$<" index 1 using 1:4 with line title "Cubic interpolation derivative from GSL"\
	'| pyxplot

out.cinterp.txt: main
	./$< > $@

main: main.c cinterp.c 

.PHONEY:clean
clean:
	$(RM) main out* plot*
