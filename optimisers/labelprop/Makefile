#!/bin/bash
CXX=g++
DIRSRC= ./
CXXFLAGS= -DDEBUG -O5 -Wall -Wno-unused-variable
EXEC=test

all: $(EXEC)

test : $(DIRSRC)test.o
	$(CXX) -o $@ $^ $(CXXFLAGS)
	
##########################################
# Generic rules
##########################################

%.o: %.cpp %.h
	$(CXX) -o  $@ -c $< $(CXXFLAGS)

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS)

clean:
	rm -f $(DIRSRC)*.o test

	
mrproper: clean
	rm -f *~ $(EXEC)
