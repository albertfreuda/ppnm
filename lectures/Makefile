CFLAGS = -Wall -O1 -std=gnu11 
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: 

.PHONEY:clean
clean:
	$(RM) out*
