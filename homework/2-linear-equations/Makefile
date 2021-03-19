CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY:default
default: timing.png out.txt

timing.png:out.time.txt
	echo '\
	set terminal png;\
	set output $@;\
	set key top left;\
	set tics out;\
	set xlabel "$$N$$";\
	set ylabel "Time";\
	set title "Timing of QR-decomposition";\
	a=1 ;\
	f(x)=a*x**3 ;\
	fit f(x) "$<" via a;\
	plot \
	"$>" with linespoints title "CPU times"\
	,f(x) with lines title "$$\propto N^3$$"\
	'|pyxplot

out.time.txt:QRtime Makefile
	cat /dev/null > $@
	for n in $$(seq 500 50 1000); do \time --format "$$n %U" --append --output $@ ./$< $$n; done

out.txt: main
	./$< > $@

QRtime:QRtime.c functions.c

main: main.c functions.c

.PHONEY:clean
clean: 
	$(RM) out* main *.png
