CFLAGS = -Wall -O1 -std=gnu11
LDLIBS = -lm

exponential_implementation.pdf:exponential_implementation.tex true_exp_plot.pdf simple_exp_plot.pdf total_exp_plot.pdf
	pdflatex $<
	pdflatex $<
#	bibtex $<
#	pdflatex $<

out.txt:main
	./$< > $@

main:main.c ex.c

true_exp_plot.pdf: true_exp_plot.gpi out.txt
	pyxplot $<
simple_exp_plot.pdf: simple_exp_plot.gpi out.txt
	pyxplot $<
total_exp_plot.pdf: total_exp_plot.gpi out.txt
	pyxplot $<
