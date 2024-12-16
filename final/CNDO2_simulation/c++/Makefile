# Define the C++ compiler
CXX = g++

# Compilation flags: -Wall for compiler warnings, -g for debugging information
CXXFLAGS = -Wall -g

# Include directories for header files
INCLUDES = -I/usr/include

# Library paths, if Armadillo is not in a standard location
LFLAGS = -L/usr/lib

# Libraries to link against, including Armadillo
LIBS = -larmadillo

# Source files
SRCS = md_simulation_cndo.cpp CNDO.cpp util.cpp AO.cpp

# Object files, derived from SRCS
OBJS = $(SRCS:.cpp=.o)

# Executable name
MAIN = md_simulation_cndo

# Targets
.PHONY: depend clean

all:    $(MAIN) run_simulation
	@echo  Simple compiler named md_simulation_cndoz has been compiled

$(MAIN): $(OBJS) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# Rule to convert .cpp to .o
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

run_simulation: $(MAIN)
	./$(MAIN) H2O.txt 1000

# DO NOT DELETE
