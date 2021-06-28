#Prepare compiler, and load gsl
CFLAGS = -Wall -std=gnu11 -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

#Create default output target
.PHONEY:default
default: out.txt

out.txt: main
	./$< > $@

main: main.c cholesky.c

#Define clean target
.PHONEY:clean
clean:
	$(RM) out* main
