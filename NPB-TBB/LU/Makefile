SHELL=/bin/sh
BENCHMARK=lu
BENCHMARKU=LU

include ../config/make.def

OBJS = lu.o ${COMMON}/c_print_results.o \
       ${COMMON}/c_timers.o ${COMMON}/c_wtime.o

include ../sys/make.common

${PROGRAM}: config ${OBJS}
	${CLINK} ${CLINKFLAGS} -o ${PROGRAM} ${OBJS} ${C_LIB}

lu.o:		lu.cpp npbparams.hpp
	${CCOMPILE} lu.cpp

clean:
	- /bin/rm -f npbparams.hpp
	- /bin/rm -f *.o *~
