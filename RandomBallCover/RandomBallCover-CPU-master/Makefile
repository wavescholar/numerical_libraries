CC=gcc
CCFLAGS= -fopenmp -funroll-loops -msse -Wall -O3
LINKFLAGS= -lgsl -lgslcblas -lm
EXACT_EXECUTABLE=exactRBC
EXACT_DRIVER=exactDriver.c
ONESHOT_EXECUTABLE=oneShotRBC
ONESHOT_DRIVER=oneShotDriver.c
SOURCES= utils.c brute.c rbc.c dists.c
OBJECTS=$(SOURCES:.c=.o) 
OS_OBJECTS=$(ONESHOT_DRIVER:.c=.o) 
E_OBJECTS=$(EXACT_DRIVER:.c=.o)

all: $(ONESHOT_EXECUTABLE) $(EXACT_EXECUTABLE)

$(ONESHOT_EXECUTABLE): $(OBJECTS) $(OS_OBJECTS)
	$(CC) $(CCFLAGS) $(OBJECTS) $(OS_OBJECTS) -o $@ $(LINKFLAGS)

$(EXACT_EXECUTABLE): $(OBJECTS) $(E_OBJECTS)
	$(CC) $(CCFLAGS) $(OBJECTS) $(E_OBJECTS) -o $@ $(LINKFLAGS)

%.o:%.c
	$(CC) $(CCFLAGS) -c $+ 

clean:
	rm -f *.o
	rm -f $(EXECUTABLE)
