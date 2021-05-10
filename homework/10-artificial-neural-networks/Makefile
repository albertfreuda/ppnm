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
	set key top left;\
	set title "Deeplearning interpolation";\
	f(x) = x*exp(-x);\
	plot \
	"$<" index 0 using 1:2 with points title "Data"\
	,"$<" index 1 using 1:2 with line title "Interpolant"\
	' | pyxplot

out.txt:main
	./$< > $@

main: main.c ann.c newton.c amoeba.c

.PHONEY:clean
clean: 
	$(RM) out* main

