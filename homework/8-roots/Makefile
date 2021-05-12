CFLAGS = -Wall -std=gnu11 -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY:default
default: out.txt out.s_wave.png out.better_boundary.png

out.s_wave.png: out.s_wave.txt Makefile
	echo '\
	set terminal png;\
	set output $@;\
	set xlabel "$$r$$";\
	set ylabel "$$\psi_r$$";\
	set key top right;\
	set title "B: S-wave hydrogen solution";\
	f(x) = x*exp(-x);\
	plot \
	"$<" using 1:2 with line title "ODE S-wave"\
	, f(x) with line title "Real S-wave"\
	' | pyxplot

out.better_boundary.png: out.better_boundary.txt Makefile
	echo '\
	set terminal png;\
	set output $@;\
	set xlabel "$$r_{max}$$";\
	set ylabel "Energy";\
	set key top right;\
	set title "Convergence of energy";\
	f(x) = -0.5;\
	plot \
	"$<" using 1:2 with line title "$$\psi_r(r_{max}) = 0$$"\
	,"$<" using 1:3 with line title "$$\psi_r(r_{max}) = re^{-kr}$$"\
	, f(x) with line title "Real energy"\
	' | pyxplot

out.better_boundary.txt: partC
	./$< > $@

out.txt: main
	./$< > $@

main: main.c newton.c runge_kutta.c linsolve.c

partC: partC.c newton.c runge_kutta.c linsolve.c

.PHONEY:clean
clean: 
	$(RM) main partC out*
