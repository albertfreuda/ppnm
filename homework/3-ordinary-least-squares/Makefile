CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: fit.png 

fit.png:out.txt Makefile
	echo '\
	set terminal png;\
	set output $@;\
	set key top right;\
	set tics out;\
	set xlabel "Time [days]";\
	set ylabel "Activity";\
	set title "Decay";\
	f(x) =exp(4.95867-0.171379*x);\
	g(x) =exp(4.95867-0.0298375-(0.171379+0.00353553)*x);\
	h(x) =exp(4.95867+0.0298375-(0.171379-0.00353553)*x);\
	q(x) =exp(4.95867-(0.19254)*x);\
	plot \
	"$<" index 0 with yerrorbars title "Data"\
	,f(x) with color red line title "Fit"\
	,g(x) with color black line lt 2 title "Lower bound"\
	,h(x) with color black line lt 2 title "Upper bound"\
	,q(x) with color red line lt 2 title "Newest estimate"\
	'| pyxplot

out.txt: main
	./$< > $@

main: main.c functions.c ls_fit.c 

.PHONEY:clean
clean: 
	$(RM) out* main *.png
