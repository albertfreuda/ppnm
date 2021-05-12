CFLAGS = -O1 -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm 
LDLIBS += $(shell gsl-config --libs)

.PHONEY:default
default: plot.linterp.png plot.linterp_integ.png

plot.linterp.png: out.linterp.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "$$x$$";\
	set ylabel "$$y$$";\
	set title "Interpolation";\
	plot \
	"$<" index 1 using 1:2 with line title "Linear interpolation"\
	,"$<" index 1 using 1:6 with line title "Linear interpolation (gsl)"\
	,"$<" index 1 using 1:4 with line title "Real function $$y=x^2$$"\
	,"$<" index 0 using 1:2 with points title "Data"\
	'| pyxplot

plot.linterp_integ.png: out.linterp.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "$$x$$";\
	set ylabel "$$y$$";\
	set title "Interpolation (integrals)";\
	plot \
	"$<" index 1 using 1:3 with line title "Integral of linear interpolation"\
	,"$<" index 1 using 1:7 with line title "Integral of linear interpolation (gsl)"\
	,"$<" index 1 using 1:5 with line title "Real integral $$y=x^3/3$$"\
	'| pyxplot

out.linterp.txt: main
	./$< > $@

main: main.c linterp.c linterp_int.c

.PHONEY:clean
clean:
	$(RM) main out* *.png
