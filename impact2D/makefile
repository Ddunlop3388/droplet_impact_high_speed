CC = qcc

#Define number of cores to use
NUM_CORES = 1

#Compiler Prefix flags
CPFLAGS = -Wall -O2 -g 

#Compiler suffix flags
CSFLAGS = -L$(BASILISK)/gl -lglutils -lfb_tiny -lm

#Define the source files
SRC = droplet_impact.c

#Define the object files
OBJS = $(SRCS:.c=.o)

#Define the executable name
EXEC = droplet_impact

#Default rule
all: $(EXEC) clean

#Link to the executable
$(EXEC): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(EXEC) $(CSFLAGS)
	#mpirun -np $(NUM_CORES) ./$(EXEC)	
	./$(EXEC)
#Clean rule to remove generated files
clean:
	rm -f $(OBJS) $(EXEC)
