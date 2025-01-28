#
#
#
import ROOT as RT
RT.gROOT.ProcessLine('#include "../src/Geom.h"')
RT.gROOT.ProcessLine('#include "../src/Lens.h"')
RT.gROOT.ProcessLine('#include "../src/Source.h"')
RT.gROOT.ProcessLine('#include "../src/Motion.h"')
RT.gROOT.ProcessLine('#include "../src/Fractal.h"')
RT.gROOT.ProcessLine('#include "../src/Lcurve.h"')
RT.gROOT.ProcessLine('#include "../src/RT-Graph.h"')
RT.gSystem.Load('../glib/libMLFR.so')
RT.gSystem.Load('../glib/libRTGR.so')

