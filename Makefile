CC=gcc
CXX=g++
CXXFLAGS=-g -O3
#CXX=g++
#CXXFLAGS=-O3 -g
INCLUDE= -I${HOME}/daten/workdir/head/lime-1.3.2/include/ -I. -I.. 
#LIBS=-L${BETA_LIBRARY_32_PATH} -L${HOME}/daten/workdir/lime-1.2.3-32/lib/ -lboost_program_options-gcc -llime -lpthread -lm
LIBS=-lboost_program_options -L${HOME}/daten/workdir/head/lime-1.3.2/lib/ -llime -lm

all: vector_ff ppcor

io.o: io.cc io.hh Makefile
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o $@ $<

local_vec.o: %.o: %.cc %.hh geometry.hh Makefile
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o $@ $<

gamma.o: gamma.cc gamma.hh Makefile
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o $@ $<

vector_ff.o: vector_ff.cc  Makefile local_vec.hh io.hh gamma.hh geometry.hh
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o $@ $<

ppcor.o: ppcor.cc Makefile io.hh
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o $@ $<

vector_ff: vector_ff.o local_vec.o io.o gamma.o dml.o DML_crc32.o 
	${CXX} ${CXXFLAGS} -o $@ vector_ff.o io.o local_vec.o gamma.o dml.o DML_crc32.o ${LIBS}

ppcor: ppcor.o dml.o DML_crc32.o
	${CXX} ${CXXFLAGS} -o $@ ppcor.o io.o dml.o DML_crc32.o ${LIBS}

DML_crc32.o: DML_crc32.cc
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o $@ $<

dml.o: dml.cc
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o $@ $<

clean:
	rm -f *.o *.d *~ vector_ff
