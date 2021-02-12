CC = gcc
CFLAGS = -O -std=gnu11
CFLAGS += -Wall
LDLIBS = -lm

out.txt: main
	./$< > $@

clean:
	$(RM) main out.txt
