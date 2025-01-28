#include "RT-Graph.h"

int main() {
  dMlens ml41 = {dLens(0.9, dV2d(0.0, 0.0)), dLens(0.04, dV2d(1.13, 0.11)),
                 dLens(0.04, dV2d(0.98, -0.21)),
                 dLens(0.02, dV2d(1.22, -0.22))};

  dSource src01(dCircle(dV2d(0.1, 0.00), 0.01));

  Time<double> t0(0.0);

  dLinear ln1(-0.3, t0, M_PI / 2.0, 20.0);

  dSrcMotion sm01(src01, ln1);

  dLcurve lc41_01({sm01}, dMlMotion(ml41), t0.range(200, -20.0, 20.0));

  lc41_01.process_while(1000);

  RTLcurve rtlc(lc41_01);

  auto mag = rtlc.mag();

  TCanvas c("Light curve example", "Light curve example", 800, 400);

  mag.Draw();

  c.SaveAs("lc41_01.png");

  c.Close();
}
