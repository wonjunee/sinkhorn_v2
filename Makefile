# CC := g++
CC := g++ # typical homebrew compiler name 

# Check for environment definitions of compiler 
# e.g., on CC = g++-6 on OSX
ifdef PHYSICELL_CPP 
	CC := $(PHYSICELL_CPP)
endif

 
CFLAGS     :=  -O3 -std=c++11 -I./src -Wall
# FFTW library
FFTWFLAG   := -I/usr/local/include  -L/usr/local/lib -lfftw3
# opencv library: for image
OPENCVFLAG := `pkg-config --cflags --libs opencv4`
# for parallelization
# OPENMPFLAG := -Xclang -fopenmp -lomp
 
COMPILE_COMMAND := $(CC) $(CFLAGS) $(FFTWFLAG) $(OPENCVFLAG) $(OPENMPFLAG)
 
OUTPUT := bfm
 
all: src/main.cpp
	$(COMPILE_COMMAND) -o $(OUTPUT) src/main.cpp

clean:
	rm -f *.o $(OUTPUT).*
