CFLAGS = -Wall -std=gnu11 -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY:default
default: out.txt out.energies.png out.modes.png out.scaling.png

out.energies.png: out.qm.txt Makefile
	echo '\
	set terminal png;\
	set output $@;\
	set title "Eigen-energies of particle in a box";\
	plot \
	"$<" index 0 using 1:2 with points title "Analytic"\
	,"$<" index 0 using 1:3 with points title "Numeric"\
	'| pyxplot

out.modes.png: out.qm.txt Makefile
	echo '\
	set terminal png;\
	set output $@;\
	set key bottom right;\
	set title "Eigen-modes of particle in a box";\
	plot \
	"$<" index 1 using 1:2 with line title "Ground State"\
	,"$<" index 1 using 1:3 with line title "1st"\
	,"$<" index 1 using 1:4 with line title "2nd"\
	'| pyxplot

out.scaling.png: out.scaling.txt out.gsl_scaling.txt out.optimized.txt Makefile
	echo '\
	set terminal png;\
	set output $@;\
	set key top left;\
	set tics out;\
	set xlabel "$$N$$";\
	set ylabel "Time";\
	set title "Timing of Jacobi eigenvalue decomposition";\
	a=1;\
	f(x) = a*x**3;\
	fit f(x) "$<" via a;\
	plot \
	"$<" with linespoints title "CPU times (own routine)"\
	,f(x) with lines title "$$\propto N^3$$"\
	,"out.gsl_scaling.txt" with linespoints title "CPU times (gsl)"\
	,"out.optimized.txt" with linespoints title "CPU times (optimized)"\
	'| pyxplot

out.txt: main
	./$< > $@

out.qm.txt: particle_in_a_box
	./$< > $@

out.scaling.txt:scaling Makefile
	cat /dev/null > $@
	for n in $$(seq 100 10 200);do \time --format "$$n %U" --append --output $@ ./$< $$n; done
	
out.gsl_scaling.txt:gsl_scaling Makefile
	cat /dev/null > $@
	for n in $$(seq 100 10 200);do \time --format "$$n %U" --append --output $@ ./$< $$n; done

out.optimized.txt:optimized Makefile
	cat /dev/null > $@
	for n in $$(seq 100 10 200);do \time --format "$$n %U" --append --output $@ ./$< $$n; done

main: main.c functions.c

scaling: scaling.c functions.c
gsl_scaling: gsl_scaling.c functions.c
optimized: optimized.c functions.c

particle_in_a_box: particle_in_a_box.c functions.c

.PHONEY:clean
clean:
	$(RM) out* main particle_in_a_box scaling gsl_scaling
