CC = gcc
CFLAGS = -O -std=gnu11
CFLAGS += -Wall
LDLIBS = -lm

out.txt: main
	./$< > $@

main: main.o equal.o
	$(CC) $^ -o $@ $(LDLIBS)

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

equal.o: equal.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) main out.txt *.o
