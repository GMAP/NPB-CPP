SHELL=/bin/sh
BENCHMARK=ep
BENCHMARKU=EP

include ../config/make.def

OBJS = ${BENCHMARK}${VERSION}.o ${COMMON}/c_print_results.o ${COMMON}/c_${RAND}.o \
       ${COMMON}/c_timers.o ${COMMON}/c_wtime.o

include ../sys/make.common

${PROGRAM}: config ${OBJS}
	${CLINK} ${CLINKFLAGS} -o ${PROGRAM} ${OBJS} ${C_LIB}

${BENCHMARK}${VERSION}.o:		${BENCHMARK}${VERSION}.cpp npbparams.hpp
	${CCOMPILE} ${BENCHMARK}${VERSION}.cpp

clean:
	- rm -f *.o *~ 
	- rm -f npbparams.hpp core
