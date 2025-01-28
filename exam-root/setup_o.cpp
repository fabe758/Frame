{
  // gROOT->ProcessLine(".L set_lim.cpp");
  gROOT->ProcessLine("#include \"../src/Geom.h\"");
  gROOT->ProcessLine("#include \"../src/Lens.h\"");
  gROOT->ProcessLine("#include \"../src/Source.h\"");
  gROOT->ProcessLine("#include \"../src/Motion.h\"");
  gROOT->ProcessLine("#include \"../src/Fractal.h\"");
  gROOT->ProcessLine("#include \"../src/Lcurve.h\"");
  gROOT->ProcessLine("#include \"../src/RT-Graph.h\"");
  gROOT->ProcessLine("#include <Math/SpecFuncMathMore.h>");
  //  gROOT->ProcessLine("#include \"Test.h\"");
  // gSystem->Load("/usr/lib/x86_64-linux-gnu/libstdc++.so.6");
  gSystem->Load("../glib/libMLFR.so");
  gSystem->Load("../glib/libRTGR.so");
  // gSystem->Load("/usr/lib/root/libMathMore.so");
  //  gSystem->Load("libTest.so");
}
