#
#  Makefile
#

CXX = g++
# CXX = clang++
# NVCXX = nvc++

# CXXFLAGS = -fPIC -I.. -std=c++17
CXXOPT = -O2
# CXXOPT = -O3
CXXFLAGS = -fPIC -I.. -I$(ROOTSYS)/include
CXXWARN = -Wall
SHARED = -shared
# NFLAGS = -fPIC -shared -acc -Minfo=accel

SRCS = Geom.cpp Lens.cpp Source.cpp Fractal.cpp Motion.cpp Lcurve.cpp
OBJS = Geom.o Lens.o Source.o Fractal.o Motion.o Lcurve.o
TARGET1 = ../glib/libMLFR.so
TARGET2 = ../glib/libRTGR.so

all: MLFR RTGR

MLFR: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(SHARED) -o $(TARGET1)

RTGR: $(OBJS) $(TARGET1) RT-Graph.o
	$(CXX) $(CXXFLAGS) $(SHARED) RT-Graph.o -o $(TARGET2)

Geom.o: Geom.cpp Geom.h
	$(CXX) -c $(CXXWARN) $(CXXOPT) $(CXXFLAGS) Geom.cpp

Motion.o:  Motion.cpp Motion.h
	$(CXX) -c $(CXXWARN) $(CXXOPT) $(CXXFLAGS) Motion.cpp

Lens.o: Lens.cpp Lens.h
	$(CXX) -c $(CXXWARN) $(CXXOPT) $(CXXFLAGS) Lens.cpp

Source.o: Source.cpp Source.h
	$(CXX) -c $(CXXWARN) $(CXXOPT) $(CXXFLAGS) Source.cpp

Fractal.o:  Fractal.cpp Fractal.h
	$(CXX) -c $(CXXWARN) $(CXXOPT) $(CXXFLAGS) Fractal.cpp

Lcurve.o:  Lcurve.cpp Lcurve.h
	$(CXX) -c $(CXXWARN) $(CXXOPT) $(CXXFLAGS) Lcurve.cpp

RT-Graph.o: RT-Graph.h RT-Graph.cpp
	$(CXX) -c $(CXXWARN) $(CXXOPT) $(CXXFLAGS) RT-Graph.cpp

clean:
	rm -f $(OBJS) RT-Graph.o $(TARGET1) $(TARGET2)


