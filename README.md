Fractal algorithm in gravitational microlensing
===============================================

**Fumio Abe**

Here,  we present a new algorithm in modeling of gravitational microlensing.
This method is particulary usefull for multiple planet analysis.

# Requirement
C++ compiler  : g++, clang++
ROOT, pyroot (https://root.cern.ch) : interpreter and library for C++ and interface to python
If you don't use ROOT and pyroot, omit compiling RT-Graph.cpp

# Directories

src         : Source directory of libraries
glib        : Library directory
exam-root   : Example directory for ROOT
test-python : test directory for python

# To make glib
$ cd src
$ make
Then, libMLFR.so and libRTGR.so are generated in glib directory

# To try c++ example
$ cd exam-c++
$ make
Then,ccca41 and lc4101 are generated. 
To run those program, 
$ source setup.sh
Then, 
$ ./ccca41
etc.

# To run root example
$ cd exam-root
$ root
root> .L setup_o.cpp         or  .L setup_s.cpp   (if no glib)
root> .L Exam.cpp
root> lc41_01.process_while(1000)
root> RTLcurve rtlc(lc41_01)
root> auto mag = rtlc.mag()
root> mag.Draw()








