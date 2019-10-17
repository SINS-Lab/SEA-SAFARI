CXX = g++

#-Ofast cuts runtime to approximately 1/3.
CXXFLAGS = -Wall -std=c++11 -Ofast

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