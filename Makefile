## Directories and Files
PACKAGE := $(notdir ${CURDIR})
SRCDIR := .
OBJDIR := build
INCLUDEDIR := -I/usr/local/include -I${HOME}/local/include
PROGRAM := a.out


## Directories and Files
SRCS := $(wildcard ${SRCDIR}/*.cpp)
OBJS := $(notdir $(SRCS:.cpp=.o))
OBJS := $(addprefix ${OBJDIR}/,${OBJS})
vpath %.cpp ${SRCDIR}


## Options
GXX := $(notdir $(firstword $(foreach x,g++-6 g++-5 g++,$(shell which $x 2>/dev/null))))
CXX_ARRAY := clang++ ${GXX}
CXX := $(firstword $(foreach x,${CXX_ARRAY},$(shell which $x)))
CC := $(CXX)
CPPFLAGS := -Wall -Wextra -Wno-unused-parameter -fno-strict-aliasing ${INCLUDEDIR} ${CPPDBG} -ftemplate-depth=512
CXXFLAGS := -std=c++14 -O3 ${CXXDBG}
LDFLAGS = -L${HOME}/local/lib -L${BOOST}/lib -Wl,-rpath,${BOOST}/lib
LDLIBS := -lsfmt -lboost_program_options -lboost_filesystem -lboost_system -lboost_iostreams -lboost_zlib
TARGET_ARCH := -march=core2 -m64 -msse -msse2 -msse3

ifneq (,$(filter $(CXX), ${GXX}))
  CXXFLAGS += -mfpmath=sse
		BOOST := ${HOME}/local/boost-gcc
  CPPFLAGS += -pthread -D_GLIBCXX_USE_NANOSLEEP
else
  CXXFLAGS += -stdlib=libc++
		BOOST := ${HOME}/local/boost-clang
  ifeq ($(shell uname), Linux)
    LDLIBS += -lpthread -lsupc++
  endif
endif

export CXX CC TARGET_ARCH

## Targets
.PHONY: all clean run test help
.DEFAULT_GOAL := all

all:
	${MAKE} -j3 ${PROGRAM}

${PROGRAM}: ${OBJS}
	${LINK.cpp} ${OUTPUT_OPTION} $^ ${LOADLIBES} ${LDLIBS}

clean:
	${RM} ${OBJS} ${PROGRAM}

run:
	@./${PROGRAM}

test:
	./${PROGRAM} --test

help:
	./${PROGRAM} --help


.PHONY: debug release instruments
debug:
	${MAKE} CXXDBG="-g" all

release:
	${MAKE} CPPDBG="-DNDEBUG" all

instruments: release
	@echo $(shell gdate +%F_%T)
	instruments -t "/Applications/Xcode.app/Contents/Applications/Instruments.app/Contents/Resources/templates/Time Profiler.tracetemplate" -D ~/tmp/profile$(shell gdate +%F_%T) ${PROGRAM} -T200


.PHONY: html pdf pandoc
html:
	$(RM) -r docs/*
	doxygen

pdf:
	pdflatex --output-directory=tex tex/anolis.tex && open tex/anolis.pdf

pandoc:
	pandoc tex/anolis.tex -s --mathjax -o docs/model.html

docs: html pandoc
	@:

${OBJDIR}/%.o: | ${OBJDIR}
	$(COMPILE.cpp) ${OUTPUT_OPTION} $<

${OBJDIR}:
	mkdir $@

## Dependencies
-include Dependfile

.PHONY: Depend
Depend:
	${CXX} -MM ${CPPFLAGS} ${CXXFLAGS} ${TARGET_ARCH} ${SRCS} | sed 's|\(\w*\.o:\)|${OBJDIR}/\1|' > Dependfile
