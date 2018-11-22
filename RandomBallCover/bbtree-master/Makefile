CC=gcc
CCFLAGS=-Wall -O3
LINKFLAGS=-lm
SOURCES=search.c utils.c bregdiv.c bbtree.c test.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=testBBT
all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CCFLAGS) $(OBJECTS) -o $@ $(LINKFLAGS)

%.o:%.c
	$(CC) $(CCFLAGS) -c $+

clean:
	rm -rf *.o
	rm -rf testBBT
