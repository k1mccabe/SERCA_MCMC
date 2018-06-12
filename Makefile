# Makefile for Writing Make Files Example

# *****************************************************
# Variables to control Makefile operation

CXX = mpic++
CXXFLAGS = -Wall -g -I/usr/include/openmpi-x86_64 #for use on amarolab computers
# ****************************************************
# Targets needed to bring the executable up to date
all: main 

main: main.o get_Residual.o update_States.o lastRun.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm *.o $(objects) main

