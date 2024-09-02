SHELL=/bin/sh
BENCHMARK=is
BENCHMARKU=IS

include ../config/make.def

OBJS = is.o ${COMMON}/c_print_results.o ${COMMON}/c_${RAND}.o \
       ${COMMON}/c_timers.o ${COMMON}/c_wtime.o

include ../sys/make.common

${PROGRAM}: config ${OBJS}
	${CLINK} ${CLINKFLAGS} -o ${PROGRAM} ${OBJS} ${C_LIB}


is.o:		is.cpp npbparams.hpp
	${CCOMPILE} is.cpp

clean:
	- rm -f *.o *~ 
	- rm -f npbparams.hpp core
