# Makefile for SSgff

INCLUDE= ./include
CDIR= ./src
OBJ = ./objects
BIN = ./bin

HEADERS = $(INCLUDE)/ssgff.h
PROGRAM= ssgff
PRODUCT= $(BIN)/$(PROGRAM)

#normal compile with optimization
#OPT1: gcc
CC= gcc
OPTS= -I$(INCLUDE)  -gdwarf-2 -g3 -Wall -Wextra  
#-Wshadow -Wstrict-overflow -fno-strict-aliasing

#OPTS=-I$(INCLUDE) -Wall -Wextra -O3

#OPT2: clang
#scan-build make
#CC= clang
#OPTS= -I$(INCLUDE) -std=c11  -gdwarf-2 -g3 -Wpedantic -Wall -Wextra  -Wshadow -Wstrict-overflow -fno-strict-aliasing

#OPT3: Sun Studio cc
#CC= cc
#OPTS= -I$(INCLUDE) -std=c11

#-Wmisleading-indentation -Werror -std=c11

#######

OBJECTS = $(OBJ)/input.o $(OBJ)/output.o $(OBJ)/scan.o $(OBJ)/dictionary.o $(OBJ)/readargv.o

#######

$(PRODUCT): $(BIN) $(OBJ) $(OBJ)/$(PROGRAM).o $(OBJECTS) $(HEADERS)
	$(CC) $(OPTS) -o $(PRODUCT) $(OBJ)/$(PROGRAM).o $(OBJECTS) -lm

$(BIN) :
	mkdir $(BIN);

$(OBJ) :
	mkdir $(OBJ);

$(OBJ)/$(PROGRAM).o : $(CDIR)/$(PROGRAM).c  $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/$(PROGRAM).c -o $(OBJ)/$(PROGRAM).o

$(OBJ)/input.o :  $(CDIR)/input.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/input.c -o $(OBJ)/input.o

$(OBJ)/output.o :  $(CDIR)/output.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/output.c -o $(OBJ)/output.o

$(OBJ)/scan.o : $(CDIR)/scan.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/scan.c -o $(OBJ)/scan.o

$(OBJ)/dictionary.o : $(CDIR)/dictionary.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/dictionary.c -o $(OBJ)/dictionary.o

$(OBJ)/readargv.o : $(CDIR)/readargv.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/readargv.c -o $(OBJ)/readargv.o

clean:
	rm -f $(OBJ)/*.o $(PRODUCT) *~ $(INCLUDE)/*~ $(CDIR)/*~ core;
#	rmdir $(BIN) $(OBJ);
