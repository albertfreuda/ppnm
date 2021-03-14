CFLAGS = -O1 -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm 
LDLIBS += $(shell gsl-config --libs)

.PHONEY:default
defaul: plot.linterp.png plot.linterp_gsl.png

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
	"$<" using 1:2 with line title "Linear interpolation"\
	,"$<" using 1:3 with line title "Integral of linear interpolation"\
	,"$<" using 1:4 with line title "Real function $$y=x^2$$"\
	,"$<" using 1:5 with line title "Real integral $$y=x^3/3$$"\
	'| pyxplot

plot.linterp_gsl.png: out.linterp.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "$$x$$";\
	set ylabel "$$y$$";\
	set title "Interpolation (gsl)";\
	plot \
	"$<" using 1:6 with line title "Linear interpolation (gsl)"\
	,"$<" using 1:7 with line title "Integral of linear interpolation (gsl)"\
	,"$<" using 1:4 with line title "Real function $$y=x^2$$"\
	,"$<" using 1:5 with line title "Real integral $$y=x^3/3$$"\
	'| pyxplot

out.linterp.txt: main
	./$< > $@

main: main.c linterp.c linterp_int.c

.PHONEY:clean
clean:
	$(RM) main out*
