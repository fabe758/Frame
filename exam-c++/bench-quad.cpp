//
//
//
#include "Lcurve.h"
#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TMarker.h>
#include <TStyle.h>
#include <string>

int main(int argc, char *argv[]) {
  std::string file11 = "figs/cont-quad.png";
  std::string file12 = "figs/cont-quad.pdf";
  std::string file21 = "figs/cont-quad2.png";
  std::string file22 = "figs/cont-quad2.pdf";
  if (argc > 1)
    std::string file = argv[1];

  // Qyad lens
  Mlens<double> ml = {Lens<double>(0.9, V2d<double>(0.0, 0.0)),
                      Lens<double>(0.04, V2d<double>(1.13, 0.11)),
                      Lens<double>(0.04, V2d<double>(0.98, -0.21)),
                      Lens<double>(0.02, V2d<double>(1.22, -0.22))};
  // Circular uniform source at (0.1, 0.01),  radius is 0.01
  Source<double> src(Circle<double>(V2d<double>(0.1, 0.01), 0.01));
  Source<double> src001(Circle<double>(V2d<double>(0.1, 0.01), 0.001));

  Time<double> t0(0.0);
  Time<float> ft0(0.0);
  Time<long double> lt0(0.0);
  auto t = t0.range(1000, -100.0, 100.0);
  auto ft = ft0.range(1000, -100.0, 100.0);
  auto lt = lt0.range(1000, -100.0, 100.0);
  dLinear ln(-0.3, t0, M_PI / 2.0, 20.0);
  fLinear fln(-0.3, ft0, M_PI / 2.0, 20.0);
  lLinear dln(-0.3, lt0, M_PI / 2.0, 20.0);
  FrParam<double> par10;
  par10.rd_init = 10;

  dSrcMotion sm(src, ln);
  // dSrcMotion sm(src001, ln);
  fSrcMotion fsm(src.to_f(), ln.to_f());
  lSrcMotion lsm(src.to_l(), ln.to_l());

  dLcurve lc({sm}, ml, t, par10);
  lc.set_bench(true);
  fLcurve flc({fsm}, ml.to_f(), ft, par10.to_f());
  flc.set_bench(true);
  lLcurve llc({lsm}, ml.to_l(), lt, par10.to_l());
  llc.set_bench(true);

  // lc.process(15);
  // flc.process(15);
  // llc.process(15);

  lc.process(22);
  flc.process(22);
  llc.process(22);

  // lc.process(25);
  // flc.process(25);
  // llc.process(25);

  auto dct = lc.ctime();
  auto fct = flc.ctime();
  auto lct = llc.ctime();

  std::cout << " CPU time in double = " << dct << std::endl;
  std::cout << " CPU time in float = " << fct << std::endl;
  std::cout << " CPU time in long double = " << lct << std::endl;
  std::cout << "long double : start = " << llc.lcp().start
            << ", stop = " << llc.lcp().stop << std::endl;
  std::cout << "CLOCKS_PER_SEC = " << CLOCKS_PER_SEC << std::endl;

  t.clear();
  ft.clear();
  lt.clear();
  lc.clear();
  flc.clear();
  llc.clear();
  ml.clear();

  return 0;
}
