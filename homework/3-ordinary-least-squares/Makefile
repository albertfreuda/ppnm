CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: out.txt


out.txt: main
	./$< > $@

main: main.c functions.c 

.PHONEY:clean
clean: 
	$(RM) out* main *.png
