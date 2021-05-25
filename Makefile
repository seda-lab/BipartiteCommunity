#!/bin/bash
CXX=g++
DIRSRC= ./
CXXFLAGS= -O5 -Wall -Wno-unused-variable
EXEC=main

all: $(EXEC)

main : $(OBJ1) $(DIRSRC)main.o
	$(CXX) -o $@ $^ $(CXXFLAGS)
	
##########################################
# Generic rules
##########################################

%.o: %.cpp %.h
	$(CXX) -o  $@ -c $< $(CXXFLAGS)

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS)

clean:
	rm -f $(DIRSRC)*.o main

mrproper: clean
	rm -f *~ $(EXEC)
