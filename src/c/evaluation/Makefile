# Makefile for Evaluation

INCLUDE= ./include
CDIR= ./src
OBJ = ./objects
BIN = ./bin

HEADERS = $(INCLUDE)/evaluation.h 
PROGRAM= evaluation
PRODUCT= bin/$(PROGRAM)
CC= gcc 
#OPTS=-I$(INCLUDE) -Wall 

#scan-build make 
#CC= clang
#OPTS=-I$(INCLUDE) -Wall -Wextra -gdwarf-2 -g3 -ggdb
OPTS=-I$(INCLUDE) -std=c11 -gdwarf-2 -g2 -Wpedantic -Wall -Wextra  -Wshadow -Wstrict-overflow -fno-strict-aliasing  
#-Wduplicated-cond 


#######
MKDIR_P = mkdir -p

.PHONY: directories

all: directories ${OBJECTS} ${PRODUCT}

directories: ${OBJ} ${BIN}

${OBJ}:
	${MKDIR_P} ${OBJ}

${BIN}:
	${MKDIR_P} ${BIN}

OBJECTS = $(OBJ)/readargv.o $(OBJ)/Output.o\
	$(OBJ)/ReadExonsGFF.o $(OBJ)/ComputeTPe.o\
	$(OBJ)/ComputeTP.o $(OBJ)/ComputeTPg.o\
#######

$(PRODUCT): $(OBJ)/$(PROGRAM).o $(OBJECTS) $(HEADERS)
	$(CC) $(OPTS) -o $(PRODUCT) $(OBJ)/$(PROGRAM).o $(OBJECTS) -lm

$(OBJ)/$(PROGRAM).o : $(CDIR)/$(PROGRAM).c  $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/$(PROGRAM).c -o $(OBJ)/$(PROGRAM).o 

$(OBJ)/readargv.o :  $(CDIR)/readargv.c $(HEADERS) 
	$(CC) -c $(OPTS) $(CDIR)/readargv.c -o $(OBJ)/readargv.o

$(OBJ)/Output.o :  $(CDIR)/Output.c $(HEADERS) 
	$(CC) -c $(OPTS) $(CDIR)/Output.c -o $(OBJ)/Output.o

$(OBJ)/ReadExonsGFF.o : $(CDIR)/ReadExonsGFF.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/ReadExonsGFF.c -o $(OBJ)/ReadExonsGFF.o

$(OBJ)/ComputeTP.o : $(CDIR)/ComputeTP.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/ComputeTP.c -o $(OBJ)/ComputeTP.o

$(OBJ)/ComputeTPe.o : $(CDIR)/ComputeTPe.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/ComputeTPe.c -o $(OBJ)/ComputeTPe.o

$(OBJ)/ComputeTPg.o : $(CDIR)/ComputeTPg.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/ComputeTPg.c -o $(OBJ)/ComputeTPg.o

clean:
	rm -f $(OBJ)/*.o $(PRODUCT) *~ $(INCLUDE)/*~ $(CDIR)/*~ core;


