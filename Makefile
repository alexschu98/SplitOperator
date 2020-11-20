UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
	CC	= g++
	INCLUDE   = /usr/include
	LIBCC     = /usr/lib
	LIBFLAG   = -lfftw3
	CCFLAGS   = -std=c++11 -march=native -O2 -Wall
endif
ifeq ($(UNAME), Darwin)
	CC	= clang++
	INCLUDE   = /usr/local/include
	LIBCC     = /usr/local/lib
	LIBFLAG   = -lfftw3
	CCFLAGS   = -std=c++11 -m64 -march=native -O2 -Wall
endif

all: splitop1d splitop2d splitop3d

splitop1d:
	$(CC) -o bin/splitop1d src/SplitOperator1D.cpp -I$(INCLUDE) -L$(LIBCC) $(CCFLAGS) $(LIBFLAG)

splitop2d:
	$(CC) -o bin/splitop2d src/SplitOperator2D.cpp -I$(INCLUDE) -L$(LIBCC) $(CCFLAGS) $(LIBFLAG)

splitop3d:
	$(CC) -o bin/splitop3d src/SplitOperator3D.cpp -I$(INCLUDE) -L$(LIBCC) $(CCFLAGS) $(LIBFLAG)
