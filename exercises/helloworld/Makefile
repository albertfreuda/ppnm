CC = gcc
CFLAGS = -O -std=gnu11
LDLIBS = -lm


default: out.txt
	cat out.txt

out.txt: hello
	./hello > out.txt


hello: hello.o
	$(CC) -o hello hello.o $(LDLIBS)


hello.o: hello.c
	$(CC) $(CFLAGS) -c hello.c

.PHONEY: clean
clean: 
	$(RM) -f hello.o hello out.txt


test:
	echo $(LDLIBS)
	echo $(CC)
