CXX = g++
CXXFLAGS = -Wall -std=c++11 -O2
SRCDIR = src
OUTDIR = bin/Debug

_FILES = main.cpp \
      safio.cpp vec_math.cpp space_math.cpp\
      lattice.cpp ion.cpp \
      potentials.cpp hameq.cpp \
      scat.cpp
FILES = $(patsubst %,$(SRCDIR)/%,$(_FILES))

OUTNAME = Sea-Safari
OUTPUT = $(patsubst %,$(OUTDIR)/%,$(OUTNAME))

Sea-Safari: $(FILES)
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) $(FILES)