SHELL := /bin/bash


#default compiler settings
CC = g++
OPT = -O3 -g -std=c++0x
LDFLAGS = -lm

# set pattern conversion name
Gen_EXE = test
Gen_SRC = main.cc
Gen_OBJ = Test.o

# compilation for runs
all: track generate link

 track : CLCT.cc CLCT.h
	$(CC) $(OPT) -c CLCT.cc -o CLCT.o

 generate : main.cc CLCT.cc CLCT.h
	$(CC) $(OPT) -c main.cc -o main.o
#	$(CC)  $(OPT) $(Gen_SRC) -o $(Gen_EXE) $(LDFLAGS) 

 link : CLCT.o main.o
	$(CC) main.o CLCT.o -o $(Gen_EXE) -lm
	/bin/rm -rf *.o

clean:
	rm -rf *.o $(Gen_EXE)
