#
#
#

# CXX = g++
CXX = g++
# CXX = clang++
CXXOPT = -O2
# CXXOPT = -O3
CXXFLAGS = $(CXXOPT) -Wall -I../src -I$(ROOTSYS)/include -std=c++17
LNLIBS = $(shell root-config --libs --glibs) -L../glib -lMLFR -lRTGR

all: ccca41 lc4101 
	
ccca41: ccca41.o
		$(CXX) $(CXXFLAGS) ccca41.o -o ccca41 $(LNLIBS)

lc4101: lc4101.o
		$(CXX) $(CXXFLAGS) lc4101.o -o lc4101 $(LNLIBS)

ccca41.o : ccca41.cpp
		$(CXX) -c $(CXXFLAGS) ccca41.cpp -o ccca41.o  $(LNLIBS)

lc4101.o : lc4101.cpp
		$(CXX) -c $(CXXFLAGS) lc4101.cpp -o lc4101.o  $(LNLIBS)

clean: 
	rm -f ccca41.o  lc4101.o \
	ccca41 lc4101 \
	*.png

