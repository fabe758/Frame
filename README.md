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
$ cd src <br>
$ make <br>
Then, libMLFR.so and libRTGR.so are generated in glib directory

# To try c++ example
$ cd exam-c++ <br>
$ make <br>
Then,ccca41 and lc4101 are generated. <br>
To run those program, <br>
$ source setup.sh <br>
Then, <br>
$ ./ccca41 <br>
etc.

# To run root example
$ cd exam-root <br>
$ root <br>
root> .L setup_o.cpp         or  .L setup_s.cpp   (if no glib) <br>
root> .L Exam.cpp <br>
root> lc41_01.process_while(1000) <br>
root> RTLcurve rtlc(lc41_01) <br>
root> auto mag = rtlc.mag() <br>
root> mag.Draw() <br>








