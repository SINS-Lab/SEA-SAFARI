CXX = g++

CXXFLAGS_D = -Wall -std=c++11 -O4 -fopenmp -march=native -m64 -msse4.2 -mavx2 -ffast-math -mfpmath=sse -freciprocal-math -ffinite-math-only -mrecip=all -mstackrealign -mpc64 -malign-double -pg
CXXFLAGS_X = -Wall -std=c++11 -O4 -fopenmp -march=native -m64 -msse4.2 -mavx2 -ffast-math -mfpmath=sse -freciprocal-math -ffinite-math-only -mrecip=all -mstackrealign -mpc64 -malign-double
CXXFLAGS_R = -Wall -std=c++11 -O4 -fopenmp -march=native -m64 -msse4.2 -mavx2 -ffast-math -mfpmath=sse -freciprocal-math -ffinite-math-only -mrecip=all -mstackrealign -mpc64 -malign-double

SRCDIR = src
OUTDIR_D = bin/Debug
OUTDIR_A = analysis
OUTDIR_R = bin/Release

_FILES_SAFARI_R = main.cpp \
                safio.cpp string_utils.cpp \
                vec_math.cpp space_math.cpp \
                lattice.cpp ion.cpp \
                potentials.cpp hameq.cpp hameq_cascade.cpp \
                scat.cpp traj.cpp tests.cpp temps.cpp traj_cascade.cpp
FILES_SAFARI_R = $(patsubst %,$(SRCDIR)/%,$(_FILES_SAFARI_R))

_FILES_SAFARI_D = main.cpp \
                safio.cpp string_utils.cpp \
                vec_math.cpp space_math.cpp \
                lattice.cpp ion.cpp \
                potentials.cpp hameq.cpp hameq_cascade.cpp \
                scat.cpp traj.cpp tests.cpp temps.cpp traj_cascade.cpp
FILES_SAFARI_D = $(patsubst %,$(SRCDIR)/%,$(_FILES_SAFARI_D))

_FILES_XYZ = xyz_process.cpp xyz.cpp \
             string_utils.cpp vec_math.cpp
FILES_XYZ = $(patsubst %,$(SRCDIR)/%,$(_FILES_XYZ))

OUTNAME = Sea-Safari
OUTPUT_D = $(patsubst %,$(OUTDIR_D)/%,$(OUTNAME))
OUTPUT_R = $(patsubst %,$(OUTDIR_R)/%,$(OUTNAME))
OUTNAME_Z = XYZ
OUTPUT_X = $(patsubst %,$(OUTDIR_A)/%,$(OUTNAME_Z))

#xyz: XYZ_PROCESSOR
#debug: Sea-Safari-Debug
all: XYZ_PROCESSOR Sea-Safari-Debug Sea-Safari-Release

Sea-Safari-Debug:
	$(CXX) $(CXXFLAGS_D) -o $(OUTPUT_D) $(FILES_SAFARI_D)

Sea-Safari-Release:
	$(CXX) $(CXXFLAGS_R) -o $(OUTPUT_R) $(FILES_SAFARI_R)

XYZ_PROCESSOR:
	$(CXX) $(CXXFLAGS_X) -o $(OUTPUT_X) $(FILES_XYZ)
