CXX = g++

#-Ofast cuts runtime to approximately 1/3.
CXXFLAGS_D = -Wall -std=c++11 -O3 -pg
CXXFLAGS_R = -Wall -std=c++11 -O3

SRCDIR = src
OUTDIR_D = bin/Debug
OUTDIR_R = bin/Release

_FILES_SAFARI = main.cpp \
                safio.cpp string_utils.cpp \
                vec_math.cpp space_math.cpp \
                lattice.cpp ion.cpp \
                potentials.cpp hameq.cpp \
                scat.cpp traj.cpp
FILES_SAFARI = $(patsubst %,$(SRCDIR)/%,$(_FILES_SAFARI))

_FILES_XYZ = xyz_process.cpp xyz.cpp \
             string_utils.cpp
FILES_XYZ = $(patsubst %,$(SRCDIR)/%,$(_FILES_XYZ))

OUTNAME = Sea-Safari
OUTPUT_D = $(patsubst %,$(OUTDIR_D)/%,$(OUTNAME))
OUTPUT_R = $(patsubst %,$(OUTDIR_R)/%,$(OUTNAME))

#xyz: XYZ_PROCESSOR
all: Sea-Safari-Debug Sea-Safari-Release XYZ_PROCESSOR

Sea-Safari-Debug: $(FILES_SAFARI)
	$(CXX) $(CXXFLAGS_D) -o $(OUTPUT_D) $(FILES_SAFARI)

Sea-Safari-Release: $(FILES_SAFARI)
	$(CXX) $(CXXFLAGS_R) -o $(OUTPUT_R) $(FILES_SAFARI)

XYZ_PROCESSOR: $(FILES_XYZ)
	$(CXX) $(CXXFLAGS_R) -o $(OUTDIR_R) $(FILES_XYZ)