CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

erf.pyxplot.png: data.txt tab_val_erf.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key bottom right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	"tab_val_erf.txt" with points title "Tabulated values"\
	,"$<" using 1:2 with line title "erf from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-erf"\
	,"$<" using 1:4 with line title "homemade Erf"\
	'| pyxplot


gamma.pyxplot.png: data2.txt tab_val_gamma.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Gamma function";\
	plot \
	"tab_val_gamma.txt" with points title "Tabulated values"\
	,"$<" using 1:2 with line title "tgamma from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-gamma"\
	,"$<" using 1:4 with line title "Stirlings app."\
	'| pyxplot

lngamma.pyxplot.png: dataB.txt tab_val_lngamma.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "lnGamma function";\
	plot \
	"tab_val_lngamma.txt" with points title "Tabulated values"\
	,"$<" using 1:2 with line title "gamma from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-lngamma"\
	,"$<" using 1:4 with line title "Stirlings app."\
	'| pyxplot

data.txt: main
	./$< > $@

main: main.o Erf.o

data2.txt: main2
	./$< > $@

main2: main2.o Gamma.o

dataB.txt: mainB
	./$< > $@

mainB: mainB.o Gamma.o

.PHONEY:clean
clean: 
	$(RM) mainB out* main data2.txt data.txt *.o *.pn *.png main2

test:
	echo $(CFLAGS)
	echo $(LDLIBS)
