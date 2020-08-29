SHELL=/bin/sh
BENCHMARK=sp
BENCHMARKU=SP

include ../config/make.def

OBJS =	sp.o \
	${COMMON}/c_print_results.o ${COMMON}/c_timers.o ${COMMON}/c_wtime.o

include ../sys/make.common

${PROGRAM}: config ${OBJS}
	${CLINK} ${CLINKFLAGS} -o ${PROGRAM} ${OBJS} ${C_LIB}

sp.o:		sp.cpp npbparams.hpp
	${CCOMPILE} sp.cpp

clean:
	- rm -f *.o *~ mputil*
	- rm -f npbparams.hpp core
