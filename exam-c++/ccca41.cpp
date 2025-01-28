#include "RT-Graph.h"

int main() {
  RTCanvas rtc;
  dMlens ml41 = {dLens(0.9, dV2d(0.0, 0.0)), dLens(0.04, dV2d(1.13, 0.11)),
                 dLens(0.04, dV2d(0.98, -0.21)),
                 dLens(0.02, dV2d(1.22, -0.22))};

  std::vector<RTCC_Sel> sel41 = {
      {Square<double>(V2d<double>(0.1, 0.0), 1.5), {{21, 10.0}, {5, 500.0}}},
      {Square<double>(V2d<double>(1.0, 0.0), 0.5), {{21, 2.0}, {5, 100.0}}}};
  RTMlens rtml(ml41);
  auto grccca = rtml.grccca(sel41);
  rtc.frame = Square<double>(V2d<double>(0.1, 0.0), 1.5);
  RTMkCanvas(rtc);
  rtc.file = "ccca41.png";
  grccca.cc.SetMarkerColor(1);
  grccca.ca.SetMarkerColor(2);
  grccca.lens.SetMarkerColor(4);
  grccca.cc.Draw("p");
  grccca.ca.Draw("p");
  grccca.lens.Draw("*");
  RTSvCanvas(rtc);
}
