pkginclude_HEADERS = Rabit.h gsl_util.h lin.h lm.h util.h
instdir = $(prefix)/lib

inst_LIBRARIES = librabit.a
librabit_a_SOURCES=Rabit.cpp Rabit.h gsl_util.c gsl_util.h lin.c lin.h lm.c lm.h util.cpp util.h

bin_PROGRAMS=Rabit
Rabit_SOURCES=main.cpp Matrixmap.cpp Matrixmap.h
Rabit_LDADD = librabit.a -lgsl -lgslcblas -llapack -lblas -lgfortran
