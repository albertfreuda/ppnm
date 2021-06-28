#Prepare compiler, and load gsl
CFLAGS = -Wall -std=gnu11 -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

#Create default output target and pdf-file for uploading to eksamen.au.dk
.PHONEY:default
default: out.txt exam_project.pdf

#Create pdf-file
exam_project.pdf: exam_project.tex
	pdflatex $<
	pdflatex $<

#Create output file
out.txt: main
	./$< > $@

#Link and compile c-code
main: main.c cholesky.c

#Define clean target
.PHONEY:clean
clean:
	$(RM) out* main 
