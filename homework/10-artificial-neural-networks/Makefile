CFLAGS = -Wall -std=gnu11 -O
CFLAGS = $(shell gsl-config --cflags)
LDLIBS += -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY:default
default: out.spline.png

out.spline.png:out.txt Makefile
	echo '\
	set terminal png;\
	set output $@;\
	set xlabel "$$x$$";\
	set ylabel "$$y$$";\
	set key outside;\
	set title "Deeplearning interpolation";\
	f(x) = x*exp(-x);\
	plot \
	"$<" index 0 using 1:2 with points title "Data"\
	,"$<" index 1 using 1:2 with line title "Interpolant"\
	,"$<" index 1 using 1:3 with line title "Interpolant integral"\
	,"$<" index 1 using 1:4 with line title "Interpolant derivative"\
	' | pyxplot

out.txt:main
	./$< > $@

main: main.c ann.c newton.c amoeba.c

.PHONEY:clean
clean: 
	$(RM) out* main

