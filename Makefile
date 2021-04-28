# CC := g++
CC := g++ # typical homebrew compiler name 

# Check for environment definitions of compiler 
# e.g., on CC = g++-6 on OSX
ifdef PHYSICELL_CPP 
	CC := $(PHYSICELL_CPP)
endif

ifndef FILE
	FILE = src/main.cpp
endif

# To run Debugging mode, `make "DEBUG=-DDEBUG"`
 
CFLAGS     :=  -O3 -std=c++11 -I./src -Wall -g
# FFTW library
FFTWFLAG   := -I/usr/local/include  -L/usr/local/lib -lfftw3
# opencv library: for image
OPENCVFLAG := `pkg-config --cflags --libs opencv4`
# for parallelization
# OPENMPFLAG := -Xclang -fopenmp -lomp
 
COMPILE_COMMAND := $(CC) $(CFLAGS) $(FFTWFLAG) $(OPENCVFLAG) $(OPENMPFLAG) $(DEBUG)
 
OUTPUT := sinkhorn
 
all: 
	$(COMPILE_COMMAND) -o $(OUTPUT) $(FILE)

clean:
	rm -f *.o $(OUTPUT).*
