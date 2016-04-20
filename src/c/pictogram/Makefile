#Makefile

#CC= gcc
#CFLAGS= -O
#scan-build make 
CC= clang
CFLAGS= -O -std=c11 -gdwarf-2 -g2 -march=native -Wpedantic -Wall -Wextra  -Wshadow -Wstrict-overflow -fno-strict-aliasing

IFLAGS= -I./include/
LFLAGS= -lm

BIN= .

.KEEP_STATE:

PROGRAMS= 	pictogram

.c.o:
	$(CC) $(CFLAGS) $(IFLAGS) -c $<

all: $(PROGRAMS)

pictogram :	pictogram.c
		$(CC) $(CFLAGS) $(IFLAGS) pictogram.c $(LFLAGS) -o $(BIN)/pictogram



