# Makefile for Writing Make Files Example

# *****************************************************
# Variables to control Makefile operation

CXX = g++
CXXFLAGS = -Wall -g

# ****************************************************
# Targets needed to bring the executable up to date
all: main 

main: main.o get_Residual.o update_States.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm *.o $(objects) main

