# Makefile

OBJS	= bison.o lex.o main.o
OBJS2	= bison2.o lex2.o main2.o

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    FFLAGS = -lfl
else
    FFLAGS = -ll
endif

CC	= g++
CFLAGS	= -g -Wall -ansi -pedantic

TARGETS: calc calc2

calc:		$(OBJS)
		$(CC) $(CFLAGS) $(OBJS) -o calc $(FFLAGS)

calc2:		$(OBJS2)
		$(CC) $(CFLAGS) $(OBJS2) -o calc2 $(FFLAGS)

lex.o:		lex.c
		$(CC) $(CFLAGS) -c lex.c -o lex.o

lex.c:		calc.lex 
		flex calc.lex
		cp lex.yy.c lex.c

lex2.o:		lex.c
		$(CC) $(CFLAGS) -c lex.c -o lex2.o

bison.o:	bison.c
		$(CC) $(CFLAGS) -c bison.c -o bison.o

bison.c:	calc.y
		bison -d -v calc.y
		cp calc.tab.c bison.c
		cmp -s calc.tab.h tok.h || cp calc.tab.h tok.h

bison2.o:	bison2.c
		$(CC) $(CFLAGS) -c bison2.c -o bison2.o

bison2.c:	calc2.y
		bison -d -v calc2.y
		cp calc2.tab.c bison2.c
		cmp -s calc2.tab.h tok2.h || cp calc2.tab.h tok2.h

main.o:		main.cc
		$(CC) $(CFLAGS) -c main.cc -o main.o

main2.o:		main.cc
		$(CC) $(CFLAGS) -c main.cc -o main2.o

lex.o yac.o main.o	: heading.h
lex.o main.o		: tok.h
lex2.o main2.o : tok2.h

clean:
	rm -f *.o *~ lex.c lex.yy.c bison.c bison2.c tok.h tok2.h calc.tab.c calc2.tab.c calc.tab.h calc2.tab.h calc.output calc calc2 calc2.output
