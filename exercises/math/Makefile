CC = gcc
CFLAGS = -std=gnu99 -O -Wall
LDLIBS = -lm

out.txt: ex_math
	./ex_math > out.txt


ex_math: ex_math.o
	$(CC) $(LDFLAGS) $< -o $@ $(LDLIBS)


ex_math.o: ex_math.c
	$(CC) $(CFLAGS) -c ex_math.c -o ex_math.o


clean:
	$(RM) *.o out* ex_mat
