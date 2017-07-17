TARGETS=test_csc lookup_test csc_learn_test
SRC1=$(wildcard src/*.c*)
SRC2=$(SRC1:.c=.o)
OBJ=$(SRC2:.cpp=.o)


EXTLIBS=-lgsl -lgslcblas
LDFLAGS=$(EXTLIBS)
SOURCES=$(wildcard src/*.h)
#CFLAGS=-O2 -pg -Wall -DPARALLEL for profiling

ifdef DEBUG
CFLAGS=-DDEBUG -ggdb -O0 -Wall -mtune=nocona -fno-inline -fno-omit-frame-pointer # -ggdb: maximum debug info -pg: profiling for gprof
export PARALLEL=
else
CFLAGS=-O2 -Wall -mtune=nocona
endif

ifdef USE_CBLAS
CFLAGS+= -DUSE_CBLAS
endif

ifdef PARALLEL
#error EEE
CFLAGS+= -fopenmp
EXTLIBS+= -lgomp
endif

ifdef PROFILE
  LDFLAGS+= -pg
  CFLAGS+= -pg
else
  ifdef DEBUG
  CFLAGS+= -fno-omit-frame-pointer
  endif
endif

all: $(TARGETS)

%.o:%.c
	gcc $(CFLAGS) -c $< -o $@

%.o:%.cpp
	g++ $(CFLAGS) -c $< -o $@


csc_learn_test: csc_learn_test.cpp $(OBJ) $(SOURCES)
	g++ $(CFLAGS) $(LDFLAGS) -o $@  $< $(OBJ)  $(EXTLIBS)

lookup_test: lookup_test.cpp $(OBJ) $(SOURCES)
	g++ $(CFLAGS) $(LDFLAGS) -o $@  $< $(OBJ)  $(EXTLIBS)

test_csc: test_csc.cpp $(OBJ) $(SOURCES)
	g++ $(LDFLAGS) -o $@  $< $(OBJ)  $(EXTLIBS)

clean:
	-rm $(OBJ) $(TARGETS)  2> /dev/null
