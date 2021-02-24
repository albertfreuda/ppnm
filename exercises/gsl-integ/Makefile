CC = gcc
CFLAGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: out.txt bessel.png erf.png

bessel.png: out.bessel.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Bessel function of first kind";\
	plot \
	"$<" using 1:2 with line title "n=0"\
	,"$<" using 1:3 with line title "n=1"\
	,"$<" using 1:4 with line title "n=2"\
	,"$<" using 1:5 with line title "n=3"\
	'| pyxplot

erf.png: out.erf.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key bottom right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	"$<" with line title "erf from integral form"\
	'| pyxplot

out.erf.txt: funB
	./$< > $@

out.bessel.txt: funC
	./$< > $@

out.txt: funA
	./$< > $@

test:
	echo $(CFLAGS)
	echo $(LDLIBS)

.PHONEY: clean
clean:
	$(RM) funA funB funC *.png out*
