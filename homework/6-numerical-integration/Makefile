CFLAGS = -Wall -std=gnu11 -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY:default
default: out.evaluations.png

out.evaluations.png: out.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "$$\pi_{real}$$ = 3.1415926535897932384626433...";\
	set logscale x;\
	set logscale y;\
	set key outside;\
	set xlabel "Tolerance";\
	set ylabel "Number of evaluations";\
	plot \
	"$<" index 1 using 1:2 with points title "$$\pi$$        = 3.141592653589793115997963"\
	,"$<" index 2 using 1:2 with points title "$$\pi_{CC}$$  = 3.141592653589793115997963"\
	,"$<" index 3 using 1:2 with points title "$$\pi_{gsl}$$ = 3.141592653589792671908754"\
	' | pyxplot

out.txt: main
	./$< > $@

main: main.c quad.c

.PHONEY:clean
clean:
	$(RM) out* main 
