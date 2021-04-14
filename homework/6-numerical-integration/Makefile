CFLAGS = -Wall -std=gnu11 -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY:default
default: out.txt

out.txt: main
	./$< > $@

main: main.c quad.c

.PHONEY:clean
clean:
	$(RM) out* main 
