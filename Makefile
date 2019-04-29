#!/bin/bash
CXX=g++
CXXFLAGS= -ansi -O5 -Wall -fpermissive #-Wno-unused-variable
DIRSRC= ./src/
EXEC=louvain 
OBJ1= $(DIRSRC)graph.o $(DIRSRC)graph_binary.o $(DIRSRC)louvain.o $(DIRSRC)quality.o $(DIRSRC)modularity.o $(DIRSRC)bimodularity.o $(DIRSRC)barbermodularity.o 

all: $(EXEC)

louvain : $(OBJ1) $(DIRSRC)main_louvain.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

##########################################
# Generic rules
##########################################

%.o: %.cpp %.h
	$(CXX) -o  $@ -c $< $(CXXFLAGS)

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS)

clean:
	rm -f $(DIRSRC)*.o louvain

mrproper: clean
	rm -f *~ $(EXEC)
