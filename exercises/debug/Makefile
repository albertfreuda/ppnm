CFLAGS = -O -w -std=gnu99
CFLAGS += `gsl-config --cflags`
# Add more flags for debugging
CFLAGS += -Wall -Werror -DNDEBUG
LDLIBS += $(shell gsl-config --libs)

all: output ; cat $<
output: main ; ./$< > $@
clean: ; $(RM) main output

