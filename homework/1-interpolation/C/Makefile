CFLAGS = -O1 -Wall -std=gnu11
LDLIBS = -lm

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
	"$<" using 1:2 with line title "Quadratic interpolation"\
	,"$<" using 1:5 with line title "Real function $$y=e^{-x^2}$$"\
	,"$<" using 1:3 with line title "Integral of interpolation"\
	,"$<" using 1:6 with line title "Real integral (erf from math.h)"\
	,"$<" using 1:4 with line title "Derivative of interpolation"\
	,"$<" using 1:7 with line title "Real derivative $$y=-2xe^{-x^2}$$"\
	'| pyxplot

out.cinterp.txt: main
	./$< > $@

main: main.c cinterp.c 

.PHONEY:clean
clean:
	$(RM) main out*
