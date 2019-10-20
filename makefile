CXX = g++

#-Ofast cuts runtime to approximately 1/3.
CXXFLAGS_D = -Wall -std=c++11 -O3 -pg
CXXFLAGS_R = -Wall -std=c++11 -O3

SRCDIR = src
OUTDIR_D = bin/Debug
OUTDIR_R = bin/Release

_FILES = main.cpp \
      safio.cpp vec_math.cpp space_math.cpp\
      lattice.cpp ion.cpp \
      potentials.cpp hameq.cpp \
      scat.cpp traj.cpp
FILES = $(patsubst %,$(SRCDIR)/%,$(_FILES))

OUTNAME = Sea-Safari
OUTPUT_D = $(patsubst %,$(OUTDIR_D)/%,$(OUTNAME))
OUTPUT_R = $(patsubst %,$(OUTDIR_R)/%,$(OUTNAME))

all: Sea-Safari-Release Sea-Safari-Debug

Sea-Safari-Release: $(FILES)
	$(CXX) $(CXXFLAGS_R) -o $(OUTPUT_R) $(FILES)

Sea-Safari-Debug: $(FILES)
	$(CXX) $(CXXFLAGS_D) -o $(OUTPUT_D) $(FILES)