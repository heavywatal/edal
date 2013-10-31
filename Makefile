## Directories and Files
PACKAGE := $(notdir ${CURDIR})
SRCDIR := .
OBJDIR := build
INCLUDEDIR := ${HOME}/local/include
PROGRAM := anolis.out


## Directories and Files
SRCS := $(wildcard ${SRCDIR}/*.cpp)
OBJS := $(notdir $(SRCS:.cpp=.o))
OBJS := $(addprefix ${OBJDIR}/,${OBJS})
vpath %.cpp ${SRCDIR}


## Options
GXX := g++-4.8
UNAME := $(shell uname)
ifeq (${UNAME}, Darwin)
  CXX := clang++
#/usr/local/openmpi/bin/mpic++
else
  CXX := ${GXX}
endif

CC := $(CXX)
CPPFLAGS := -Wall -Wextra -Wno-unused-parameter -fno-strict-aliasing -iquote ${INCLUDEDIR} ${CPPDBG}
CXXFLAGS := -std=c++11 -O3 ${CXXDBG}
LDFLAGS := -L${HOME}/local/lib -L/usr/local/lib
LDLIBS := -lsfmt -lboost_program_options -lboost_serialization -lboost_system -lboost_filesystem -lboost_iostreams -lboost_regex -lz
TARGET_ARCH := -march=core2 -m64 -msse -msse2 -msse3

ifneq (,$(filter $(CXX), ${GXX}))
  CXXFLAGS += -mfpmath=sse
  CPPFLAGS += -pthread -isystem /usr/local/boost-gcc/include -D_GLIBCXX_USE_NANOSLEEP
  LDFLAGS += -L/usr/local/boost-gcc/lib -static-libstdc++
else
  CPPFLAGS += -isystem /usr/local/boost-clang/include -isystem /usr/local/include/c++/v1 -ftemplate-depth=512
  CXXFLAGS += -stdlib=libc++
  LDFLAGS += -L/usr/local/boost-clang/lib 
  ifeq (${UNAME}, Linux)
    CPPFLAGS += -ftemplate-depth=512
    LDLIBS += -lpthread -lsupc++
  endif
endif

export CXX CC TARGET_ARCH

## Targets
.PHONY: all clean run test help
.DEFAULT_GOAL := all

all:
	${MAKE} -j3 ${PROGRAM}
	cp -f ${PROGRAM} a.out

${PROGRAM}: ${OBJS}
	${LINK.cpp} ${OUTPUT_OPTION} $^ ${LOADLIBES} ${LDLIBS}

clean:
	${RM} ${OBJS} ${PROGRAM} ${PROGRAM}.sh ${PROGRAM}.bc ${PROGRAM}.s

run:
	./${PROGRAM}

test:
	./${PROGRAM} --test

help:
	./${PROGRAM} --help


.PHONY: debug release parallel full clang instruments
debug:
	${MAKE} CXXDBG="-g" all

release:
	${MAKE} CXX=${GXX} CPPDBG="-DNDEBUG" all
#	${MAKE} CXXDBG="-flto" CPPDBG="-DNDEBUG" all

parallel:
	${MAKE} CXX=${GXX} CXXDBG="-fopenmp" all

full:
	${MAKE} CXX=${GXX} CXXDBG="-fopenmp" CPPDBG="-DNDEBUG" all
#	${MAKE} CXX=${GXX} CXXDBG="-fopenmp -flto" CPPDBG="-DNDEBUG" all

llvm:
	${MAKE} CXX=clang++ CPPDBG="-DNDEBUG" all

instruments: release
	@echo $(shell gdate +%F_%T)
	instruments -t "/Applications/Xcode.app/Contents/Applications/Instruments.app/Contents/Resources/templates/Time Profiler.tracetemplate" -D ~/tmp/profile$(shell gdate +%F_%T) ${PROGRAM}

${OBJDIR}/%.o: | ${OBJDIR}
	$(COMPILE.cpp) ${OUTPUT_OPTION} $<

${OBJDIR}:
	mkdir $@

## Dependencies
-include Dependfile

.PHONY: Depend
Depend:
	${CXX} -MM ${CPPFLAGS} ${CXXFLAGS} ${TARGET_ARCH} ${SRCS} | sed 's|\(\w*\.o:\)|${OBJDIR}/\1|' > Dependfile

## misc.
.PHONY: open

mainsrcs := $(addprefix ${SRCDIR}/,\
optimum.h \
stimulus.h \
environment.h \
gene.h \
individual.h \
population.h \
main.h \
optimum.cpp \
stimulus.cpp \
environment.cpp \
gene.cpp \
individual.cpp \
population.cpp \
main.cpp \
global.hpp \
)

libsrcs := $(addprefix ${INCLUDEDIR}/cxxwtils/,\
test.cpp \
debug.hpp \
iostr.hpp \
gz.hpp \
os.hpp \
algorithm.hpp \
prandom.hpp \
genetic.hpp \
grn.hpp \
)

open:
	open -a mi ${mainsrcs} ${libsrcs} python/run.py
