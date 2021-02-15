CC = gcc
CFLAGS = -std=gnu11 -O -Wall
LDLIBS = -lm

out.txt: main
	./$< > $@

main: main.o komplex.o
	$(CC) $^ -o $@ $(LDLIBS)

main.o: main.c komplex.h
	$(CC) $(CFLAGS) -c $< -o $@

komplex.o: komplex.c komplex.h
	$(CC) $(CFLAGS) -c $< -o $@

clean: 
	$(RM) *.o out* main
