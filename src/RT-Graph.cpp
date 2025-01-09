/**
 * @file RT-Graph.cpp
 +
 * @brief Class library to draw graphs with CERN ROOT
 * @author Fumio Abe
 * @date 08 July,  2024
 *
 * @details class library for making plots of data
 * @note
 *
 */

/**
 * include headers  for this project
 * all c++ headers are loaded in Geom.h
 * Geom.h,  Lens.h,  and Source.h, etc are loaded in RT-Graph.h
 *
 */
#include "RT-Graph.h"
#include <TGraph2D.h>

TGraph make1dplt(Data1d data) {
  return TGraph(data.x.size(), data.x.data(), data.y.data());
}

TGraph2D make2dplt(Data2d data) {
  return TGraph2D(data.x.size(), data.x.data(), data.y.data(), data.z.data());
}

TGraph2D *p_make2dplt(Data2d data) {
  return new TGraph2D(data.x.size(), data.x.data(), data.y.data(),
                      data.z.data());
}

void SetFrame(TVirtualPad *Pad, dV2d min, dV2d max) {
  Pad->DrawFrame(min.x(), min.y(), max.x(), max.y());
}

void SetFrame(TVirtualPad *Pad, dRectangle rec) {
  gPad->DrawFrame(rec.min().x(), rec.min().y(), rec.max().x(), rec.max().y());
}

TMarker RTV2d::tmarker() const { return TMarker(c_v2d.x(), c_v2d.y(), 2); }

TLine RTSeg2d::tline() const {
  return TLine(c_seg.ends()[0].x(), c_seg.ends()[0].y(), c_seg.ends()[1].x(),
               c_seg.ends()[1].y());
}

std::vector<TLine> RTTriangle::tline() const {
  auto en0 = RTSeg2d(c_tri.sides()[0]);
  auto en1 = RTSeg2d(c_tri.sides()[1]);
  auto en2 = RTSeg2d(c_tri.sides()[2]);
  return {en0.tline(), en1.tline(), en2.tline()};
}

TPolyLine RTTriangle::tpolyline() const {
  auto ap = c_tri.apex();
  double x[4] = {ap[0].x(), ap[1].x(), ap[2].x(), ap[0].x()};
  double y[4] = {ap[0].y(), ap[1].y(), ap[2].y(), ap[0].y()};
  return TPolyLine(4, x, y);
}

TBox RTRectangle::tbox() const {
  TBox tb(c_rect.min().x(), c_rect.min().y(), c_rect.max().x(),
          c_rect.max().y());
  tb.SetLineWidth(2);
  tb.SetLineColor(2);
  return tb;
}

TEllipse RTCircle::tellipse() const {
  return TEllipse(c_circ.cent().x(), c_circ.cent().y(), c_circ.r(), c_circ.r());
}

// RTEllipse::RTEllipse(Ellipse<double> ell) { c_ell = ell; }

TEllipse RTEllipse::tellipse() const {
  return TEllipse(c_ell.cent().x(), c_ell.cent().y(), c_ell.ab().x(),
                  c_ell.ab().y(), 0.0, 360.0, 180.0 / M_PI * c_ell.ph());
}

TMarker RTLens::marker() const {
  RTV2d rtv(c_lens.p());
  TMarker mk = rtv.tmarker();
  mk.SetMarkerSize(c_param.marker.size);
  mk.SetMarkerStyle(c_param.marker.style);
  mk.SetMarkerColorAlpha(c_param.marker.color, c_param.marker.trans);
  return mk;
}

TEllipse RTLens::er() const {
  RTCircle rtc(c_lens.circ());
  TEllipse te = rtc.tellipse();
  te.SetLineColorAlpha(c_param.er_line.color, c_param.er_line.trans);
  te.SetLineWidth(c_param.er_line.width);
  te.SetLineStyle(c_param.er_line.style);
  te.SetFillColorAlpha(c_param.er_fill.color, c_param.er_fill.trans);
  return te;
}

TGraph RTLens::ger() const {
  std::vector<double> vx, vy;
  for (double ph = 0.0; ph < 2.0 + M_PI; ph += 2.0 * M_PI / c_param.er_ndiv) {
    vx.push_back(c_lens.re() * cos(ph) + c_lens.p().x());
    vy.push_back(c_lens.re() * sin(ph) + c_lens.p().y());
  }
  TGraph gr(vx.size(), vx.data(), vy.data());
  gr.SetLineWidth(c_param.er_line.width);
  gr.SetLineColorAlpha(c_param.er_line.color, c_param.er_line.trans);
  gr.SetLineWidth(c_param.er_line.width);
  return gr;
}

RTCC_Sel2 RTLens::sel() const {
  RTCC_Sel2 sel2;
  sel2.sq = c_lens.sq(2.0);
  sel2.max = 30000;
  sel2.sel = 10000;
  return sel2;
};

std::vector<TMarker> RTMlens::markers() const {
  std::vector<TMarker> vm;
  for (size_t i = 0; i < c_ml.size(); i++) {
    RTLens rtl = RTLens(c_ml[i]);
    vm.push_back(rtl.marker());
  }
  return vm;
};

TGraph RTMlens::centers() const {
  std::vector<double> xl;
  std::vector<double> yl;
  for (size_t i = 0; i < c_ml.size(); i++) {
    xl.push_back(c_ml[i].p().x());
    yl.push_back(c_ml[i].p().y());
  }
  return TGraph(xl.size(), xl.data(), yl.data());
}

TGraph RTMlens::er(int il, int ndiv) const {
  std::vector<double> vx, vy;
  for (double ph = 0.0; ph < 2.0 * M_PI; ph += 2.0 * M_PI / ndiv) {
    vx.push_back(c_ml[il].re() * cos(ph) + c_ml[il].p().x());
    vy.push_back(c_ml[il].re() * sin(ph) + c_ml[il].p().y());
  }
  TGraph gr(vx.size(), vx.data(), vy.data());
  gr.SetLineColor(3);
  return gr;
}

TGraph RTMlens::er(int il) const { return this->er(il, 500); }

std::vector<TGraph> RTMlens::aller(int ndiv) const {
  std::vector<TGraph> vg;
  for (size_t i = 0; i < c_ml.size(); i++)
    vg.push_back(this->er(i, ndiv));
  return vg;
}

std::vector<TGraph> RTMlens::aller() const {
  std::vector<TGraph> vg;
  for (size_t i = 0; i < c_ml.size(); i++)
    vg.push_back(this->er(i));
  return vg;
}

std::vector<TGraph> RTMlens::allger() const {
  std::vector<TGraph> vg;
  for (size_t i = 0; i < c_ml.size(); i++) {
    RTLens rtl(c_ml[i]);
    vg.push_back(rtl.ger());
  }
  return vg;
}

std::vector<TEllipse> RTMlens::ellipse() const {
  std::vector<TEllipse> ve;
  for (size_t i = 0; i < c_ml.size(); i++)
    ve.push_back(
        TEllipse(c_ml[i].p().x(), c_ml[i].p().y(), c_ml[i].re(), c_ml[i].re()));
  return ve;
}

std::vector<TGraph> RTMlens::allger(RTLensParam param) const {
  std::vector<TGraph> vg;
  for (size_t i = 0; i < c_ml.size(); i++) {
    RTLens rtl(c_ml[i], param);
    vg.push_back(rtl.ger());
  }
  return vg;
}

Data2d RTMlens::dbt(Rectangle<double> rec, V2d<double> bt, int nx, int ny,
                    double zlim) const {
  Data2d d2d;
  for (double x = rec.min().x(); x < rec.max().x(); x += rec.d() / nx) {
    for (double y = rec.min().y(); y < rec.max().y(); y += rec.d() / ny) {
      double z = c_ml.db(V2d<double>(x, y), bt);
      d2d.x.push_back(x);
      d2d.y.push_back(y);
      if (z < zlim) {
        d2d.z.push_back(z);
      } else {
        d2d.z.push_back(zlim);
      }
    }
  }
  return d2d;
}

Data2d RTMlens::dbt(Rectangle<double> rec, V2d<double> bt, double zlim) const {
  return this->dbt(rec, bt, 1000, 1000, zlim);
}

Data2d RTMlens::dbt(V2d<double> bt, double zlim) const {
  return this->dbt(
      Rectangle<double>(V2d<double>(-1.5, -1.5), V2d<double>(1.5, 1.5)), bt,
      zlim);
}

Data2d RTMlens::dbt(V2d<double> bt) const {
  return this->dbt(
      Rectangle<double>(V2d<double>(-1.5, -1.5), V2d<double>(1.5, 1.5)), bt,
      0.1);
}

TGraph2D RTMlens::cont_dbt(Rectangle<double> rec, V2d<double> bt, int nx,
                           int ny, double zlim) const {
  return make2dplt(this->dbt(rec, bt, nx, ny, zlim));
}

TGraph2D *RTMlens::p_cont_dbt(Rectangle<double> rec, V2d<double> bt,
                              double zlim) const {
  return p_make2dplt(this->dbt(rec, bt, zlim));
}

TGraph2D RTMlens::cont_dbt(V2d<double> bt) const {
  return make2dplt(this->dbt(bt));
}

TGraph2D *RTMlens::p_cont_dbt(V2d<double> bt) const {
  return p_make2dplt(this->dbt(bt));
}

TGraph2D *RTMlens::p_cont_dbt(V2d<double> bt, double zlim) const {
  return p_make2dplt(this->dbt(bt, zlim));
}

TGraph RTMlens::scat_dbt(Rectangle<double> rec, V2d<double> bt, int nx, int ny,
                         double rho) const {
  Data1d data1;
  double dx = (rec.max().x() - rec.min().x()) / nx;
  double dy = (rec.max().y() - rec.min().y()) / ny;
  for (double x = rec.min().x(); x < rec.max().x(); x += dx) {
    for (double y = rec.min().y(); y < rec.max().y(); y += dy) {
      double z = c_ml.db(V2d<double>(x, y), bt);
      if (z < rho) {
        data1.x.push_back(x);
        data1.y.push_back(y);
      }
    }
  }
  TGraph gr = make1dplt(data1);
  return gr;
}

TGraph RTMlens::scat_dbt(Rectangle<double> rec, V2d<double> bt,
                         double rho) const {
  return this->scat_dbt(rec, bt, 2000, 2000, rho);
}

TGraph RTMlens::scat_dbt(V2d<double> bt, double rho) const {
  Rectangle<double> rec(V2d<double>(-1.5, -1.5), V2d<double>(1.5, 1.5));
  return this->scat_dbt(rec, bt, rho);
}

std::string RTMlens::gp_dbt(Rectangle<double> rec, V2d<double> bt, int nx,
                            int ny) const {
  std::string gp;
  gp.append("unset surface\n");
  gp.append("set contour\n");
  gp.append("set size square\n");
  gp.append("splot '-' with line\n");
  for (double x = rec.min().x(); x < rec.max().x(); x += rec.d() / nx) {
    for (double y = rec.min().y(); y < rec.max().y(); y += rec.d() / ny) {
      double z = c_ml.db(V2d<double>(x, y), bt);
      if (z > 0.1)
        z = 0.1;
      gp.append(std::to_string(x));
      gp.append(" ");
      gp.append(std::to_string(y));
      gp.append(" ");
      gp.append(std::to_string(z));
      gp.append("\n");
    }
    gp.append("\n");
  }
  gp.append("\n");
  return gp;
}

std::string RTMlens::gp_dbt(Rectangle<double> rec, V2d<double> bt) const {
  return this->gp_dbt(rec, bt, 1000, 1000);
}

std::string RTMlens::gp_dbt(V2d<double> bt) const {
  return this->gp_dbt(
      Rectangle<double>(V2d<double>(-1.5, -1.5), V2d<double>(1.5, 1.5)), bt);
}

std::string RTMlens::gp_dbt(Lens<double> lens, V2d<double> bt, int nx,
                            int ny) const {
  return this->gp_dbt(lens.rec(1.5), bt, nx, ny);
}

std::string RTMlens::gp_dbt(Lens<double> lens, V2d<double> bt) const {
  return this->gp_dbt(lens.rec(1.5), bt);
}

void RTMlens::save_gp(std::string file, std::string &data) const {
  std::ofstream out(file);
  out << data;
  out.close();
}

TMQue<double> RTMlens::ccca(std::vector<RTCC_Sel> sel) const {
  TMQue<double> ccca;
  auto ml = c_ml;
  for (size_t ia = 0; ia < sel.size(); ia++) {
    TMQue<double> tmq(sel[ia].sq, c_ml);
    for (size_t is = 0; is < sel[ia].sel.size(); is++) {
      tmq.split(sel[ia].sel[is].round, ml);
      auto tmq2 = tmq.ccca(sel[ia].sel[is].mu);
      tmq = tmq2;
    }
    ccca += tmq;
  }
  return ccca;
}

RTCCCA RTMlens::grccca(std::vector<RTCC_Sel> sel) const {
  RTCCCA gr;
  gr.lens = this->centers();
  auto tmq = this->ccca(sel);
  RTTMQue rttmq(tmq);
  gr.cc = rttmq.gr_lens();
  gr.ca = rttmq.gr_src();

  return gr;
}

std::vector<RTCC_Sel2> RTMlens::sel2() const {
  std::vector<RTCC_Sel2> sel2;
  for (size_t i = 0; i < c_ml.size(); i++) {
    RTLens l(c_ml[i]);
    sel2.push_back(l.sel());
  }
  return sel2;
}

std::vector<TMQue<double>> RTMlens::ccca2_1(int n,
                                            std::vector<RTCC_Sel2> sel2) const {
  auto ml = c_ml;
  std::vector<TMQue<double>> tmq;
  for (size_t il = 0; il < sel2.size(); il++)
    tmq.push_back(TMQue<double>(sel2[il].sq, c_ml));
  for (size_t il = 0; il < sel2.size(); il++) {
    for (int is = 0; is < n; is++) {
      while (tmq[il].size() < sel2[il].max)
        tmq[il].split(ml);
      std::sort(tmq[il].begin(), tmq[il].end(),
                [](auto a, auto b) { return a.mu() > b.mu(); });
      tmq[il].erase(tmq[il].begin() + sel2[il].sel, tmq[il].end());
    }
  }
  return tmq;
}
TMQue<double> RTMlens::ccca2(int n, std::vector<RTCC_Sel2> sel2) const {
  auto tmq = this->ccca2_1(n, sel2);
  TMQue<double> tmq2;
  for (size_t i = 0; i < tmq.size(); i++)
    tmq2 += tmq[i];
  return tmq2;
}

RTCCCA RTMlens::grccca2(int n, std::vector<RTCC_Sel2> sel2) const {
  RTCCCA gr;
  gr.lens = this->centers();
  auto tmq = this->ccca2(n, sel2);
  RTTMQue rttmq(tmq);
  gr.cc = rttmq.gr_lens();
  gr.ca = rttmq.gr_src();
  return gr;
}

RTCCCA RTMlens::grccca2(std::vector<RTCC_Sel2> sel2) const {
  return this->grccca2(5, sel2);
}

RTCCCA RTMlens::grccca2() const { return this->grccca2(5, this->sel2()); }

TMarker RTSource::mark() const {
  TMarker mk(c_src.cent().x(), c_src.cent().y(), 3);
  mk.SetMarkerStyle(c_param.marker.style);
  mk.SetMarkerColor(c_param.marker.color);
  mk.SetMarkerSize(c_param.marker.size);
  return mk;
}

TEllipse RTSource::ell() const {
  TEllipse el;
  if (c_src.shape() == Shape::Elliptical) {
    el = TEllipse(c_src.cent().x(), c_src.cent().y(), c_src.ell().ab().x(),
                  c_src.ell().ab().y(), 0.0, 360,
                  180.0 / M_PI * c_src.ell().ph());
  } else {
    el = TEllipse(c_src.cent().x(), c_src.cent().y(), c_src.circ().r());
  }
  el.SetLineStyle(c_param.ell_line.style);
  el.SetLineWidth(c_param.ell_line.width);
  el.SetLineColor(c_param.ell_line.color);
  el.SetFillStyle(c_param.ell_fill.style);
  el.SetFillColorAlpha(c_param.ell_fill.color, c_param.ell_fill.trans);
  return el;
}

TGraph2D RTSource::sb(int n, dV2d min, dV2d max) {
  double dx = (max.x() - min.x()) / n;
  double dy = (max.y() - min.y()) / n;
  Data2d data;
  for (double x = min.x(); x < max.x(); x += dx) {
    for (double y = min.y(); y < max.y(); y += dy) {
      data.x.push_back(x);
      data.y.push_back(y);
      data.z.push_back(c_src.sb(dV2d(x, y)));
    }
  }
  return make2dplt(data);
}

TGraph2D RTSource::sb(dV2d min, dV2d max) { return this->sb(300, min, max); }

TGraph2D RTSource::sb(double fact) {
  dV2d dmc(c_src.circ().r(), c_src.circ().r());
  if (c_src.shape() == Shape::Elliptical)
    dmc = dV2d(c_src.ell().ab().x(), c_src.ell().ab().x());

  dmc = dmc * fact;
  dV2d min = c_src.cent() - dmc;
  dV2d max = c_src.cent() + dmc;
  return this->sb(min, max);
}

TMarker RTSrcMotion::marker(Time<double> t) const {
  return TMarker(c_smot.p(t).x(), c_smot.p(t).y(), 29);
}

TGraph RTSrcMotion::graph(std::vector<Time<double>> t) const {
  std::vector<double> xs;
  std::vector<double> ys;
  for (size_t i = 0; i < t.size(); i++) {
    xs.push_back(c_smot.p(t[i]).x());
    ys.push_back(c_smot.p(t[i]).y());
  }
  return TGraph(xs.size(), xs.data(), ys.data());
}

std::vector<TMarker> RTSrcMotion::markers(std::vector<Time<double>> t) const {
  std::vector<TMarker> mks;
  for (size_t i = 0; i < t.size(); i++)
    mks.push_back(this->marker(t[i]));
  return mks;
}

TMarker RTLensMotion::marker(Time<double> t) const {
  return TMarker(c_lm.p(t).x(), c_lm.p(t).y(), 34);
};

TGraph RTLensMotion::graph(std::vector<Time<double>> t) const {
  std::vector<double> xl;
  std::vector<double> yl;
  for (size_t i = 0; i < t.size(); i++) {
    xl.push_back(c_lm.p(t[i]).x());
    yl.push_back(c_lm.p(t[i]).y());
  }
  return TGraph(xl.size(), xl.data(), yl.data());
};

std::vector<TMarker> RTLensMotion::markers(std::vector<Time<double>> t) const {
  std::vector<TMarker> mks;
  for (size_t i = 0; i < t.size(); i++)
    mks.push_back(this->marker(t[i]));
  return mks;
};

std::vector<TMarker> RTMlMotion::markers(Time<double> t) const {
  std::vector<TMarker> mks;
  for (size_t i = 0; c_mlm.size(); i++) {
    RTLensMotion rtlm(c_mlm[i]);
    mks.push_back(rtlm.marker(t));
  }
  return mks;
};

std::vector<TGraph> RTMlMotion::graphs(std::vector<Time<double>> t) const {
  std::vector<TGraph> grps;
  for (size_t i = 0; i < c_mlm.size(); i++) {
    RTLensMotion rtlm(c_mlm[i]);
    grps.push_back(rtlm.graph(t));
  }
  return grps;
};

std::vector<std::vector<TMarker>>
RTMlMotion::markers(std::vector<Time<double>> t) const {
  std::vector<std::vector<TMarker>> mks;
  for (size_t i = 0; i < c_mlm.size(); i++) {
    RTLensMotion rtlm(c_mlm[i]);
    mks.push_back(rtlm.markers(t));
  }
  return mks;
};

TGraph2D RTSource::sb() { return this->sb(1.3); }

TPolyLine RTTriL::tpolyline() const {
  auto p = c_tril.apex();
  double x[4] = {p[0].x(), p[1].x(), p[2].x(), p[0].x()};
  double y[4] = {p[0].y(), p[1].y(), p[2].y(), p[0].y()};
  return TPolyLine(4, x, y);
}

TPolyLine RTTriS::tpolyline() const {
  auto p = c_tris.apex();
  double x[4] = {p[0].x(), p[1].x(), p[2].x(), p[0].x()};
  double y[4] = {p[0].y(), p[1].y(), p[2].y(), p[0].y()};
  return TPolyLine(4, x, y);
}

TPolyLine RTTriMap::lens() const {
  auto l = RTTriL(c_trimap.tril());
  return l.tpolyline();
}

TPolyLine RTTriMap::src() const {
  auto l = RTTriS(c_trimap.tris());
  return l.tpolyline();
}

std::vector<TPolyLine> RTTMQue::lens() const {
  std::vector<TPolyLine> vl;
  for (size_t i = 0; i < c_tmq.size(); i++) {
    auto l = this->c_tmq[i].tril();
    auto tl = RTTriL(l);
    vl.push_back(tl.tpolyline());
  }
  return vl;
}

TGraph RTTMQue::gr_lens() const {
  std::vector<double> x;
  std::vector<double> y;
  for (size_t i = 0; i < c_tmq.size(); i++) {
    x.push_back(c_tmq[i].tril().cent().x());
    y.push_back(c_tmq[i].tril().cent().y());
  }
  return TGraph(x.size(), x.data(), y.data());
}

TGraph RTTMQue::gr_src() const {
  std::vector<double> x;
  std::vector<double> y;
  for (size_t i = 0; i < c_tmq.size(); i++) {
    x.push_back(c_tmq[i].tris().cent().x());
    y.push_back(c_tmq[i].tris().cent().y());
  }
  return TGraph(x.size(), x.data(), y.data());
}

std::vector<TPolyLine> RTTMQue::src() const {
  std::vector<TPolyLine> vs;
  for (size_t i = 0; i < c_tmq.size(); i++) {
    auto s = this->c_tmq[i].tris();
    auto ts = RTTriS(s);
    vs.push_back(ts.tpolyline());
  }
  return vs;
}

void RTMkCanvas(RTCanvas &rtc) {
  rtc.c = new TCanvas(rtc.name.data(), rtc.title.data(), rtc.width, rtc.height);
  if (rtc.frame.d() != 0.0)
    rtc.c->DrawFrame(rtc.frame.min().x(), rtc.frame.min().y(),
                     rtc.frame.max().x(), rtc.frame.max().y());
}

void RTSvCanvas(RTCanvas &rtc) { rtc.c->SaveAs(rtc.file.data()); }

void RTClCanvas(RTCanvas &rtc) { rtc.c->Close(); }

std::vector<TPolyLine> RTFrTri::in_lens() const {
  std::vector<TPolyLine> vtri = this->in().lens();
  for (auto tri : vtri) {
    tri.SetLineStyle(c_param.in_lens_line.style);
    tri.SetLineWidth(c_param.in_lens_line.width);
    tri.SetLineColor(c_param.in_lens_line.color);
    tri.SetFillStyle(c_param.in_lens_fill.style);
    tri.SetFillColorAlpha(c_param.in_lens_fill.color,
                          c_param.in_lens_fill.trans);
  }
  return vtri;
}

std::vector<TPolyLine> RTFrTri::out_lens() const {
  std::vector<TPolyLine> vtri = this->out().lens();
  for (auto tri : vtri) {
    tri.SetLineStyle(c_param.out_lens_line.style);
    tri.SetLineWidth(c_param.out_lens_line.width);
    tri.SetLineColor(c_param.out_lens_line.color);
    tri.SetFillStyle(c_param.out_lens_fill.style);
    tri.SetFillColorAlpha(c_param.out_lens_fill.color,
                          c_param.out_lens_fill.trans);
  }
  return vtri;
}

std::vector<TPolyLine> RTFrTri::unc_lens() const {
  std::vector<TPolyLine> vtri = this->uncertain().lens();
  for (auto tri : vtri) {
    tri.SetLineStyle(c_param.unc_lens_line.style);
    tri.SetLineWidth(c_param.unc_lens_line.width);
    tri.SetLineColor(c_param.unc_lens_line.color);
    tri.SetFillStyle(c_param.unc_lens_fill.style);
    tri.SetFillColorAlpha(c_param.unc_lens_fill.color,
                          c_param.unc_lens_fill.trans);
  }
  return vtri;
}

std::vector<TPolyLine> RTFrTri::in_src() const {
  std::vector<TPolyLine> vtri = this->in().src();
  for (auto tri : vtri) {
    tri.SetLineStyle(c_param.in_src_line.style);
    tri.SetLineWidth(c_param.in_src_line.width);
    tri.SetLineColor(c_param.in_src_line.color);
    tri.SetFillStyle(c_param.in_src_fill.style);
    tri.SetFillColorAlpha(c_param.in_src_fill.color, c_param.in_src_fill.trans);
  }
  return vtri;
}

std::vector<TPolyLine> RTFrTri::out_src() const {
  std::vector<TPolyLine> vtri = this->out().src();
  for (auto tri : vtri) {
    tri.SetLineStyle(c_param.out_src_line.style);
    tri.SetLineWidth(c_param.out_src_line.width);
    tri.SetLineColor(c_param.out_src_line.color);
    tri.SetFillStyle(c_param.out_src_fill.style);
    tri.SetFillColorAlpha(c_param.out_src_fill.color,
                          c_param.out_src_fill.trans);
  }
  return vtri;
}

std::vector<TPolyLine> RTFrTri::unc_src() const {
  std::vector<TPolyLine> vtri = this->uncertain().src();
  for (auto tri : vtri) {
    tri.SetLineStyle(c_param.unc_src_line.style);
    tri.SetLineWidth(c_param.unc_src_line.width);
    tri.SetLineColor(c_param.unc_src_line.color);
    tri.SetFillStyle(c_param.unc_src_fill.style);
    tri.SetFillColorAlpha(c_param.unc_src_fill.color,
                          c_param.unc_src_fill.trans);
  }
  return vtri;
}

RTFTLensObj RTFrTri::mk_lensobj() const {
  RTFTLensObj lo;
  RTMlens rtml(c_frtri.ml());
  lo.lens_mk = rtml.markers();
  lo.lens_er = rtml.allger();
  RTSource rts(c_frtri.src());
  lo.src_mk = rts.mark();
  lo.in_lens = this->in_lens();
  lo.out_lens = this->out_lens();
  lo.unc_lens = this->unc_lens();
  return lo;
}

std::vector<RTFTLensObj> RTFrTri::mk_lensobj(int round) {
  std::vector<RTFTLensObj> vlo;
  for (int i = 0; i < round; i++) {
    c_frtri.process();
    auto lo = this->mk_lensobj();
    this->set_draw_lens(lo);
    vlo.push_back(lo);
  }
  return vlo;
}

std::vector<RTFTSrcObj> RTFrTri::mk_srcobj(int round) {
  std::vector<RTFTSrcObj> vso;
  for (int i = 0; i < round; i++) {
    c_frtri.process();
    auto so = this->mk_srcobj();
    this->set_draw_src(so);
    vso.push_back(so);
  }
  return vso;
}

RTFTSrcObj RTFrTri::mk_srcobj() const {
  RTFTSrcObj so;
  RTSource rts(c_frtri.src());
  so.src_mk = rts.mark();
  so.src_sb = rts.sb();
  so.in_src = this->in_src();
  so.out_src = this->out_src();
  so.unc_src = this->unc_src();
  return so;
}

void RTFrTri::set_draw_lens(RTFTLensObj &lo) const {
  for (size_t i = 0; i < lo.lens_mk.size(); i++) {
    lo.lens_mk[i].SetMarkerStyle(c_param.lens.marker.style);
    lo.lens_mk[i].SetMarkerSize(c_param.lens.marker.size);
    lo.lens_mk[i].SetMarkerColorAlpha(c_param.lens.marker.color,
                                      c_param.lens.marker.trans);
  }
  for (size_t i = 0; i < lo.lens_er.size(); i++) {
    lo.lens_er[i].SetLineStyle(c_param.lens.er_line.style);
    lo.lens_er[i].SetLineWidth(c_param.lens.er_line.width);
    lo.lens_er[i].SetLineColorAlpha(c_param.lens.er_line.color,
                                    c_param.lens.er_line.trans);
  };
  lo.src_mk.SetMarkerStyle(c_param.src.marker.style);
  lo.src_mk.SetMarkerSize(c_param.src.marker.size);
  lo.src_mk.SetMarkerColorAlpha(c_param.src.marker.color,
                                c_param.src.marker.trans);

  for (size_t i = 0; i < lo.in_lens.size(); i++) {
    lo.in_lens[i].SetLineStyle(c_param.in_lens_line.style);
    lo.in_lens[i].SetLineWidth(c_param.in_lens_line.width);
    lo.in_lens[i].SetLineColorAlpha(c_param.in_lens_line.color,
                                    c_param.in_lens_line.trans);
    lo.in_lens[i].SetFillStyle(c_param.in_lens_fill.style);
    lo.in_lens[i].SetFillColorAlpha(c_param.in_lens_fill.color,
                                    c_param.in_lens_fill.trans);
  };
  for (size_t i = 0; i < lo.out_lens.size(); i++) {
    lo.out_lens[i].SetLineStyle(c_param.out_lens_line.style);
    lo.out_lens[i].SetLineWidth(c_param.out_lens_line.width);
    lo.out_lens[i].SetLineColorAlpha(c_param.out_lens_line.color,
                                     c_param.out_lens_line.trans);
    lo.out_lens[i].SetFillStyle(c_param.out_lens_fill.style);
    lo.out_lens[i].SetFillColorAlpha(c_param.out_lens_fill.color,
                                     c_param.out_lens_fill.trans);
  }
  for (size_t i = 0; i < lo.unc_lens.size(); i++) {
    lo.unc_lens[i].SetLineStyle(c_param.unc_lens_line.style);
    lo.unc_lens[i].SetLineWidth(c_param.unc_lens_line.width);
    lo.unc_lens[i].SetLineColorAlpha(c_param.unc_lens_line.color,
                                     c_param.unc_lens_line.trans);
    lo.unc_lens[i].SetFillStyle(c_param.unc_lens_fill.style);
    lo.unc_lens[i].SetFillColorAlpha(c_param.unc_lens_fill.color,
                                     c_param.unc_lens_fill.trans);
  }
}

void RTFrTri::set_draw_src(RTFTSrcObj &so) const {
  so.src_mk.SetMarkerStyle(c_param.src.marker.style);
  so.src_mk.SetMarkerSize(c_param.src.marker.size);
  so.src_mk.SetMarkerColorAlpha(c_param.src.marker.color,
                                c_param.src.marker.trans);
  so.src_ell.SetLineStyle(c_param.src.ell_line.style);
  so.src_ell.SetLineWidth(c_param.src.ell_line.width);
  so.src_ell.SetLineColorAlpha(c_param.src.ell_line.color,
                               c_param.src.ell_line.trans);
  so.src_ell.SetFillStyle(c_param.src.ell_fill.style);
  so.src_ell.SetFillColorAlpha(c_param.src.ell_fill.color,
                               c_param.src.ell_fill.trans);
  so.src_sb.SetLineStyle(c_param.src.sb_line.style);
  so.src_sb.SetLineWidth(c_param.src.sb_line.width);
  so.src_sb.SetLineColorAlpha(c_param.src.sb_line.color,
                              c_param.src.sb_line.trans);
  so.src_sb.SetFillStyle(c_param.src.sb_fill.style);
  so.src_sb.SetFillColorAlpha(c_param.src.sb_fill.color,
                              c_param.src.sb_fill.trans);
  for (size_t i = 0; i < so.in_src.size(); i++) {
    so.in_src[i].SetLineStyle(c_param.in_src_line.style);
    so.in_src[i].SetLineWidth(c_param.in_src_line.width);
    so.in_src[i].SetLineColorAlpha(c_param.in_src_line.color,
                                   c_param.in_src_line.trans);
    so.in_src[i].SetFillStyle(c_param.in_src_fill.style);
    so.in_src[i].SetFillColorAlpha(c_param.in_src_fill.color,
                                   c_param.in_src_fill.trans);
  };
  for (size_t i = 0; i < so.out_src.size(); i++) {
    so.out_src[i].SetLineStyle(c_param.out_src_line.style);
    so.out_src[i].SetLineWidth(c_param.out_src_line.width);
    so.out_src[i].SetLineColorAlpha(c_param.out_src_line.color,
                                    c_param.out_lens_line.trans);
    so.out_src[i].SetFillStyle(c_param.out_src_fill.style);
    so.out_src[i].SetFillColorAlpha(c_param.out_src_fill.color,
                                    c_param.out_src_fill.trans);
  }
  for (size_t i = 0; i < so.unc_src.size(); i++) {
    so.unc_src[i].SetLineStyle(c_param.unc_src_line.style);
    so.unc_src[i].SetLineWidth(c_param.unc_src_line.width);
    so.unc_src[i].SetLineColorAlpha(c_param.unc_src_line.color,
                                    c_param.unc_src_line.trans);
    so.unc_src[i].SetFillStyle(c_param.unc_src_fill.style);
    so.unc_src[i].SetFillColorAlpha(c_param.unc_src_fill.color,
                                    c_param.unc_src_fill.trans);
  }
}

void RTFrTri::draw_lens(RTCanvas &c, RTFTLensObj &lo) const {
  for (size_t i = 0; i < lo.lens_mk.size(); i++) {
    lo.lens_mk[i].Draw(c_param.lens.marker.draw_opt.data());
  }
  for (size_t i = 0; i < lo.lens_er.size(); i++) {
    lo.lens_er[i].Draw(c_param.lens.er_fill.draw_opt.data());
  };
  lo.src_mk.Draw(c_param.src.marker.draw_opt.data());

  for (size_t i = 0; i < lo.in_lens.size(); i++) {
    lo.in_lens[i].Draw(c_param.in_lens_fill.draw_opt.data());
  };
  for (size_t i = 0; i < lo.out_lens.size(); i++) {
    lo.out_lens[i].Draw(c_param.out_lens_fill.draw_opt.data());
  }
  for (size_t i = 0; i < lo.unc_lens.size(); i++) {
    lo.unc_lens[i].Draw(c_param.unc_lens_fill.draw_opt.data());
  }
}

void RTFrTri::draw_lens(RTCanvas &c, RTFTLensObj &lo, RTCCCA &ccca) const {
  for (size_t i = 0; i < lo.lens_mk.size(); i++) {
    lo.lens_mk[i].Draw(c_param.lens.marker.draw_opt.data());
  }
  ccca.cc.Draw("p");
  lo.src_mk.Draw(c_param.src.marker.draw_opt.data());

  for (size_t i = 0; i < lo.in_lens.size(); i++) {
    lo.in_lens[i].Draw(c_param.in_lens_fill.draw_opt.data());
  };
  for (size_t i = 0; i < lo.out_lens.size(); i++) {
    lo.out_lens[i].Draw(c_param.out_lens_fill.draw_opt.data());
  }
  for (size_t i = 0; i < lo.unc_lens.size(); i++) {
    lo.unc_lens[i].Draw(c_param.unc_lens_fill.draw_opt.data());
  }
}

void RTFrTri::mk_anim_lens(RTCanvas &c, std::vector<RTFTLensObj> &vlo,
                           std::string fprefix) const {
  auto sti = [](int i) {
    std::string sti = "   ";
    std::sprintf(sti.data(), "%03d", i);
    return sti;
  };
  for (size_t i = 0; i < vlo.size(); i++) {
    RTMkCanvas(c);
    c.file = fprefix + sti(i) + ".png";
    this->draw_lens(c, vlo[i]);
    RTSvCanvas(c);
    RTClCanvas(c);
  }
}

void RTFrTri::mk_anim_lens(RTCanvas &c, std::vector<RTFTLensObj> &vlo,
                           RTCCCA &ccca, std::string fprefix) const {
  auto sti = [](int i) {
    std::string sti = "   ";
    std::sprintf(sti.data(), "%03d", i);
    return sti;
  };
  for (size_t i = 0; i < vlo.size(); i++) {
    RTMkCanvas(c);
    c.file = fprefix + sti(i) + ".png";
    this->draw_lens(c, vlo[i], ccca);
    RTSvCanvas(c);
    RTClCanvas(c);
  }
}

void RTFrTri::draw_src(RTCanvas &c, RTFTSrcObj &so) const {
  so.src_ell.Draw(c_param.src.ell_fill.draw_opt.data());
}

RTTMQue RTFrTri::in() const { return RTTMQue(c_frtri.in()); }

RTTMQue RTFrTri::uncertain() const { return RTTMQue(c_frtri.uncertain()); }

RTTMQue RTFrTri::out() const { return RTTMQue(c_frtri.out()); }

void RTFrTri::process() { c_frtri.process(); }

void RTFrTri::process(int n) { c_frtri.process(n); }

RTFTRoff RTFrTri::rtftroff(int rdmin, int rdmax) const {
  RTFTRoff roff;
  int min = rdmin - c_frtri.param().rd_init;
  int max = rdmax - c_frtri.param().rd_init;
  auto ftf = c_frtri.to_f();
  auto ftd = c_frtri;
  auto ftl = c_frtri.to_l();
  for (int i = 0; i < max; i++) {
    ftf.process();
    ftd.process();
    ftl.process();
    if (ftd.round() > rdmin) {
      roff.area_f.push_back(ftf.area_in() + ftf.area_uncertain() / 2.0);
      roff.area_d.push_back(ftd.area_in() + ftd.area_uncertain() / 2.0);
      roff.area_l.push_back(ftl.area_in() + ftl.area_uncertain() / 2.0);
      roff.area_sn_f.push_back(ftf.area_sn());
      roff.area_sn_d.push_back(ftd.area_sn());
      roff.area_sn_l.push_back(ftl.area_sn());
      roff.round.push_back(ftd.round());
    }
  }
  std::vector<double> round;
  std::vector<double> arfl;
  std::vector<double> ardl;
  std::vector<double> snfl;
  std::vector<double> sndl;
  for (int i = 0; i < roff.round.size(); i++) {
    round.push_back(static_cast<double>(roff.round[i]));
    arfl.push_back(static_cast<double>(roff.area_f[i] / roff.area_l[i] - 1.0));
    sndl.push_back(static_cast<double>(roff.area_d[i] / roff.area_l[i] - 1.0));
    snfl.push_back(
        static_cast<double>(roff.area_sn_f[i] / roff.area_sn_l[i] - 1.0));
    sndl.push_back(
        static_cast<double>(roff.area_sn_d[i] / roff.area_sn_l[i] - 1.0));
  }
  roff.ar_fl = TGraph(round.size(), round.data(), arfl.data());
  roff.ar_dl = TGraph(round.size(), round.data(), ardl.data());
  roff.ar_sn_fl = TGraph(round.size(), round.data(), snfl.data());
  roff.ar_sn_dl = TGraph(round.size(), round.data(), sndl.data());
  return roff;
}

void RTDraw_images(RTCanvas &rtc, FrTri<double> &frt,
                   std::vector<Square<double>> &grp, std::string pre,
                   std::string type) {
  auto sti = [](int i) {
    std::string sti = "   ";
    std::sprintf(sti.data(), "%03d", i);
    return sti;
  };
  RTFrTri rtf(frt);
  for (size_t i = 0; i < grp.size(); i++) {
    rtc.frame = grp[i];
    RTMkCanvas(rtc);
    rtc.name = "round " + sti(frt.round()) + "; group " + sti(i);
    rtc.title = "round " + sti(frt.round()) + "; group " + sti(i);
    rtc.file = pre + "-" + sti(i) + "-" + sti(frt.round()) + "." + type;
    RTFTLensObj lo = rtf.mk_lensobj();
    rtf.set_draw_lens(lo);
    rtf.draw_lens(rtc, lo);
    RTSvCanvas(rtc);
    RTClCanvas(rtc);
  }
}

void RTDraw_images(RTCanvas &rtc, FrTri<double> frt,
                   std::vector<Square<double>> grp) {
  std::string pre = "figs/qlens";
  std::string type = "png";
  RTDraw_images(rtc, frt, grp, pre, type);
}

void RTDraw_src_tri(RTCanvas &rtc, FrTri<double> &frt,
                    std::vector<Square<double>> &grp, Square<double> sqr,
                    std::string pre, std::string type) {
  auto sti = [](int i) {
    std::string sti = "   ";
    std::sprintf(sti.data(), "%03d", i);
    return sti;
  };
  rtc.frame = sqr;
  RTFrTri rtf(frt);
  for (size_t i = 0; i < grp.size(); i++) {
    RTMkCanvas(rtc);
    rtc.name = "round " + sti(frt.round()) + "; group " + sti(i);
    rtc.title = "round " + sti(frt.round()) + "; group " + sti(i);
    rtc.file = pre + "-" + sti(i) + "-" + sti(frt.round()) + "." + type;
    RTFTSrcObj so = rtf.mk_srcobj();
    rtf.set_draw_src(so);
    rtf.draw_src(rtc, so);
    RTSvCanvas(rtc);
    RTClCanvas(rtc);
  }
}

TGraph RTLcParam::in_sz_vs_rd() const {
  std::vector<double> in_sz;
  for (auto insz : c_lcp.in_size[0])
    in_sz.push_back(static_cast<double>(insz));
  std::vector<double> rd0;
  for (auto rd : c_lcp.round[0])
    rd0.push_back(static_cast<double>(rd));
  TGraph gr(c_lcp.round[0].size(), rd0.data(), in_sz.data());
  gr.SetTitle("in_size vs round");
  return gr;
};

void RTFrPrec::add(FrTri<double> frtri) {
  c_round.push_back(static_cast<double>(frtri.round()));
  c_d.push_back(frtri.d());
  c_area_in.push_back(frtri.area_in());
  c_area_unc.push_back(frtri.area_uncertain());
  c_area_sn.push_back(frtri.area_sn());
  c_size_in.push_back(static_cast<double>(frtri.in().size()));
  c_size_unc.push_back(static_cast<double>(frtri.uncertain().size()));
};

std::vector<double> RTFrPrec::round() const { return c_round; };

std::vector<double> RTFrPrec::d() const { return c_d; };

std::vector<double> RTFrPrec::area_in() const { return c_area_in; };

std::vector<double> RTFrPrec::area_unc() const { return c_area_unc; };

std::vector<double> RTFrPrec::area_cent() const {
  std::vector<double> ac;
  for (size_t i = 0; i < c_area_in.size(); i++)
    ac.push_back(c_area_in[i] + c_area_unc[i] / 2.0);
  return ac;
};

std::vector<double> RTFrPrec::area_sn() const { return c_area_sn; };

std::vector<double> RTFrPrec::size_in() const { return c_size_in; };

std::vector<double> RTFrPrec::size_unc() const { return c_size_unc; };

std::vector<double> RTFrPrec::dar() const {
  std::vector<double> dar0;
  for (size_t i = 0; i < c_round.size(); i++)
    dar0.push_back(std::abs(this->area_cent()[i] - this->area_sn()[i]));
  return dar0;
};

std::vector<double> RTFrPrec::rer() const {
  std::vector<double> rer0;
  for (size_t i = 0; i < c_round.size(); i++)
    rer0.push_back(std::abs(
        (this->area_cent()[i] - this->area_sn()[this->area_sn().size() - 1]) /
        this->area_sn()[this->area_sn().size() - 1]));
  return rer0;
};

std::vector<double> RTFrPrec::rdar() const {
  std::vector<double> rdar0;
  for (size_t i = 0; i < c_round.size(); i++)
    rdar0.push_back(std::abs((this->area_cent()[i] - this->area_sn()[i]) /
                             this->area_sn()[i]));
  return rdar0;
};

std::vector<double> RTFrPrec::rersn() const {
  std::vector<double> rersn0;
  for (size_t i = 0; i < c_round.size(); i++)
    rersn0.push_back(std::abs(
        (this->area_sn()[i] - this->area_sn()[this->area_sn().size() - 1]) /
        this->area_sn()[this->area_sn().size() - 1]));
  return rersn0;
};

TGraph RTFrPrec::d_vs_rd() const {
  TGraph gr(this->round().size(), this->round().data(), c_d.data());
  gr.SetTitle("d vs round");
  return gr;
};

TGraph RTFrPrec::arc_vs_rd() const {
  TGraph gr(this->round().size(), this->round().data(),
            this->area_cent().data());
  gr.SetTitle("Area vs round");
  return gr;
};

TGraph RTFrPrec::arsn_vs_rd() const {
  TGraph gr(this->round().size(), this->round().data(), this->area_sn().data());
  gr.SetTitle("Estimated area vs round");
  return gr;
};

TGraph RTFrPrec::szi_vs_rd() const {
  TGraph gr(this->round().size(), this->round().data(), this->size_in().data());
  gr.SetTitle("Number of triangle pairs vs round");
  return gr;
};

TGraph RTFrPrec::szu_vs_rd() const {
  TGraph gr(this->round().size(), this->round().data(),
            this->size_unc().data());
  gr.SetTitle("Number of triangle pairs in uncertain queue vs round");
  return gr;
};

TGraph RTFrPrec::dar_vs_rd() const {
  TGraph gr(this->round().size(), this->round().data(), this->dar().data());
  gr.SetTitle("Difference of mean and estimated areas");
  return gr;
};

TGraph RTFrPrec::rdar_vs_rd() const {
  TGraph gr(this->round().size(), this->round().data(), this->rdar().data());
  gr.SetTitle("Relative difference of areas vs round");
  return gr;
};

TGraph RTFrPrec::rer_vs_rd() const {
  TGraph gr(this->round().size(), this->round().data(), this->rer().data());
  gr.SetTitle("Convergence of relative difference of estimated area");
  return gr;
};

TGraph RTFrPrec::rersn_vs_rd() const {
  return TGraph(this->round().size(), this->round().data(),
                this->rersn().data());
};

std::vector<SrcMotion<double>> RTLcurve::srcms() const { return c_lc.srcm(); };

SrcMotion<double> RTLcurve::srcm() const { return c_lc.srcm()[0]; };

Source<double> RTLcurve::src(int it) const {
  return this->srcm().src((c_lc.time()[it]));
}

Source<double> RTLcurve::src(Time<double> t) const {
  return this->srcm().src(t);
}

FrTri<double> RTLcurve::frtri(int it) const { return c_lc.frtri(it); }

FrTri<double> RTLcurve::frtri(Time<double> t) const { return c_lc.frtri(t); }

RTFrTri RTLcurve::rtfrtri(int it) const { return RTFrTri(this->frtri(it)); }

RTFrTri RTLcurve::rtfrtri(Time<double> t) const {
  return RTFrTri(this->frtri(t));
}

RTSource RTLcurve::rtsrc(int it) const {
  return RTSource(this->frtri(it).src());
}

RTSource RTLcurve::rtsrc(Time<double> t) const {
  return RTSource(this->frtri(t).src());
}

RTMlens RTLcurve::rtml(int it) const { return RTMlens(this->frtri(it).ml()); }

RTMlens RTLcurve::rtml(Time<double> t) const {
  return RTMlens(this->frtri(t).ml());
}

RTCCCA RTLcurve::rtccca2(int it) const {
  RTMlens rtml(this->ml(it));
  return rtml.grccca2();
}

RTCCCA RTLcurve::rtccca2(Time<double> t) const {
  RTMlens rtml(this->ml(t));
  return rtml.grccca2();
}

TGraph2D *RTLcurve::p_cont_dbt(int it) const {
  RTMlens rtml(this->ml(it));
  return rtml.p_cont_dbt(this->src(it).cent());
}

TGraph2D *RTLcurve::p_cont_dbt(Time<double> t) const {
  RTMlens rtml(this->ml(t));
  return rtml.p_cont_dbt(this->src(t).cent());
}

RTFrPrec RTLcurve::prec(int it) const {
  auto frtri = this->frtri(it);
  RTFrPrec prec;
  for (int i = 0; i < 30; i++) {
    frtri.process();
    prec.add(frtri);
  }
  return prec;
}

RTFrPrec RTLcurve::prec(Time<double> t) const {
  auto frtri = this->frtri(t);
  RTFrPrec prec;
  for (int i = 0; i < 30; i++) {
    frtri.process();
    prec.add(frtri);
  }
  return prec;
}

RTFTLensObj RTLcurve::lensobj(int it, int sel) const {
  auto frtri = this->frtri(it);
  for (int i = 0; i < sel; i++)
    frtri.process();
  RTFrTri rtfrtri(frtri);
  auto lensobj = rtfrtri.mk_lensobj();
  rtfrtri.set_draw_lens(lensobj);
  return lensobj;
}

RTFTLensObj RTLcurve::lensobj_w(int it, int in_size) const {
  auto frtri = this->frtri(it);
  frtri.process_while(in_size);
  RTFrTri rtfrtri(frtri);
  auto lensobj = rtfrtri.mk_lensobj();
  rtfrtri.set_draw_lens(lensobj);
  return lensobj;
}

RTFTLensObj RTLcurve::lensobj(Time<double> t, int sel) const {
  auto frtri = this->frtri(t);
  for (int i = 0; i < sel; i++)
    frtri.process();
  RTFrTri rtfrtri(frtri);
  auto lensobj = rtfrtri.mk_lensobj();
  rtfrtri.set_draw_lens(lensobj);
  return lensobj;
}

RTFTLensObj RTLcurve::lensobj_w(Time<double> t, int in_size) const {
  auto frtri = this->frtri(t);
  frtri.process_while(in_size);
  RTFrTri rtfrtri(frtri);
  auto lensobj = rtfrtri.mk_lensobj();
  rtfrtri.set_draw_lens(lensobj);
  return lensobj;
}

std::vector<RTFTLensObj> RTLcurve::mk_lensobj_s(int round) const {
  std::vector<RTFTLensObj> vlo;
  for (size_t it = 0; it < c_lc.time().size(); it++) {
    vlo.push_back(this->lensobj(it, round));
  }
  return vlo;
}

std::vector<RTFTLensObj> RTLcurve::mk_lensobj_w(int in_size) const {
  std::vector<RTFTLensObj> vlo;
  for (size_t it = 0; it < c_lc.time().size(); it++) {
    vlo.push_back(this->lensobj_w(it, in_size));
  }
  return vlo;
}

RTFTSrcObj RTLcurve::srcobj(int it, int sel) const {
  auto frtri = this->frtri(it);
  for (int i = 0; i < sel; i++)
    frtri.process();
  RTFrTri rtfrtri(frtri);
  auto srcobj = rtfrtri.mk_srcobj();
  rtfrtri.set_draw_src(srcobj);
  return srcobj;
}

RTFTSrcObj RTLcurve::srcobj(Time<double> t, int sel) const {
  auto frtri = this->frtri(t);
  for (int i = 0; i < sel; i++)
    frtri.process();
  RTFrTri rtfrtri(frtri);
  auto srcobj = rtfrtri.mk_srcobj();
  rtfrtri.set_draw_src(srcobj);
  return srcobj;
}

void RTLcurve::set_rtc(RTCanvas &rtc, int it, int il, double fact) const {
  auto ml = this->ml(it);
  rtc.frame = ml[il].sq(fact);
}

TGraph RTLcurve::srctrj() const {
  std::vector<double> x;
  std::vector<double> y;
  for (size_t i = 0; i < c_lc.time().size(); i++) {
    x.push_back(this->srcm().p(c_lc.time()[i]).x());
    y.push_back(this->srcm().p(c_lc.time()[i]).y());
  }
  TGraph gr(x.size(), x.data(), y.data());
  gr.SetTitle("Source trajectory");
  return gr;
};

void RTLcurve::draw_lens(RTCanvas &c, RTFTLensObj &lo) const {
  auto rtfrtri = this->rtfrtri(0);
  rtfrtri.set_draw_lens(lo);
  rtfrtri.draw_lens(c, lo);
}

void RTLcurve::draw_src(RTCanvas &c, RTFTSrcObj &so) const {
  auto rtfrtri = this->rtfrtri(0);
  rtfrtri.set_draw_src(so);
  rtfrtri.draw_src(c, so);
}

std::vector<RTFTLensObj> RTLcurve::mk_lensobj(int it, int round) const {
  auto rtfrtri = this->rtfrtri(it);
  return rtfrtri.mk_lensobj(round);
}

std::vector<RTFTLensObj> RTLcurve::mk_lensobj(Time<double> t, int round) const {
  auto rtfrtri = this->rtfrtri(t);
  return rtfrtri.mk_lensobj(round);
}

std::vector<RTFTSrcObj> RTLcurve::mk_srcobj(int it, int round) const {
  auto rtfrtri = this->rtfrtri(it);
  return rtfrtri.mk_srcobj(round);
}

std::vector<RTFTSrcObj> RTLcurve::mk_srcobj(Time<double> t, int round) const {
  auto rtfrtri = this->rtfrtri(t);
  return rtfrtri.mk_srcobj(round);
}

void RTLcurve::mk_anim_lens(RTCanvas &c, std::vector<RTFTLensObj> &lo,
                            std::string fprefix) const {
  auto rtfrtri = this->rtfrtri(0);
  rtfrtri.mk_anim_lens(c, lo, fprefix);
}

void RTLcurve::mk_anim_lens(RTCanvas &c, std::vector<RTFTLensObj> &lo,
                            RTCCCA &ccca, std::string fprefix) const {
  auto rtfrtri = this->rtfrtri(0);
  rtfrtri.mk_anim_lens(c, lo, ccca, fprefix);
}

void RTLcurve::show_anim_lens(RTCanvas &c, std::vector<RTFTLensObj> &lo,
                              int time) const {
  for (size_t i = 0; i < lo.size(); i++) {
    RTMkCanvas(c);
    c.c->Clear();
    // RTClCanvas(c);
    this->draw_lens(c, lo[i]);
    sleep(time);
    c.c->Clear();
  }
}

void RTLcurve::show_anim_lens(RTCanvas &c, std::vector<RTFTLensObj> lo) const {
  auto rtfrtri = this->rtfrtri(0);
  for (size_t i = 0; i < lo.size(); i++) {
    RTMkCanvas(c);
    c.c->Clear();
    rtfrtri.draw_lens(c, lo[i]);
    std::cout << "q for quit, other key to continue : ";
    std::string key;
    std::cin >> key;
    if (key == "q")
      break;
    c.c->Clear();
  }
}

MlMotion<double> RTLcurve::mlm() const { return c_lc.mlm(); };

Mlens<double> RTLcurve::ml(Time<double> t) const { return c_lc.mlm().ml(t); };

Mlens<double> RTLcurve::ml(int it) const {
  return c_lc.mlm().ml(c_lc.time()[it]);
};

Mlens<double> RTLcurve::ml() const { return c_lc.mlm().ml(c_lc.time()[0]); };

TGraph RTLcurve::gml(Time<double> t) const {
  RTMlens rtml(this->ml(t));
  return rtml.centers();
};

TGraph RTLcurve::gml() const {
  RTMlens rtml(this->ml());
  return rtml.centers();
};

RTCCCA RTLcurve::grccca(Time<double> t, std::vector<RTCC_Sel> sel) const {
  RTMlens rtml(this->ml(t));
  return rtml.grccca(sel);
};

RTCCCA RTLcurve::grccca(std::vector<RTCC_Sel> sel) const {
  RTMlens rtml(this->ml());
  return rtml.grccca(sel);
};

RTCCCA RTLcurve::grccca2() const {
  RTMlens rtml(this->ml());
  return rtml.grccca2();
}

std::vector<TGraph> RTLcurve::lenstrj() const {
  std::vector<TGraph> glens;
  std::vector<double> x;
  std::vector<double> y;
  for (size_t il = 0; il < this->ml().size(); il++) {
    for (size_t it = 0; it < c_lc.time().size(); it++) {
      x.push_back(this->ml(c_lc.time()[it])[il].p().x());
      y.push_back(this->ml(c_lc.time()[it])[il].p().y());
    }
    glens.push_back(TGraph(x.size(), x.data(), y.data()));
    x.clear();
    y.clear();
  }
  return glens;
};

std::vector<TGraph> RTLcurve::mags() const {
  std::vector<TGraph> gmag;
  std::vector<double> t;
  for (size_t i = 0; i < c_lc.time().size(); i++)
    t.push_back(c_lc.time()[i].day());
  for (size_t is = 0; is < c_lc.mag().size(); is++)
    gmag.push_back(TGraph(t.size(), t.data(), c_lc.mags()[is].data()));
  return gmag;
};

TGraph RTLcurve::mag() const {
  std::vector<double> t;
  for (size_t i = 0; i < c_lc.time().size(); i++)
    t.push_back(c_lc.time()[i].day());
  TGraph gr(t.size(), t.data(), c_lc.mag().data());
  gr.SetTitle("Mganification vs. time (day)");
  return gr;
};

std::vector<TGraph> RTLcurve::mags2() const {
  std::vector<TGraph> gmag(c_lc.mag2().size());
  std::vector<double> t;
  for (size_t i = 0; i < c_lc.time().size(); i++)
    t.push_back(c_lc.time()[i].day());
  for (size_t is = 0; is < c_lc.mag().size(); is++)
    gmag[is] = TGraph(t.size(), t.data(), c_lc.mags2()[is].data());
  return gmag;
}

TGraph RTLcurve::mag2() const {
  std::vector<double> t;
  for (size_t i = 0; i < c_lc.time().size(); i++)
    t.push_back(c_lc.time()[i].day());
  TGraph gr(t.size(), t.data(), c_lc.mag2().data());
  gr.SetTitle("Mganification vs. time (day)");
  return gr;
};

TGraph RTLcurve::as() const {
  std::vector<double> x, y;
  for (auto as : c_lc.as()) {
    x.push_back(as.x());
    y.push_back(as.y());
  }
  TGraph gr(x.size(), x.data(), y.data());
  gr.SetTitle("Astrometric shift");
  return gr;
};

TGraph RTLcurve::as2() const {
  std::vector<double> x, y;
  for (auto as : c_lc.as2()) {
    x.push_back(as.x());
    y.push_back(as.y());
  }
  TGraph gr(x.size(), x.data(), y.data());
  gr.SetTitle("Astrometric shift");
  return gr;
};

// double RTLcYoo_un(double u, double rho) {
//   double z = u / rho;
//   double k = std::min(1.0/z, 1.0);
//   auto b0 = [k](double zz){return 4.0 / M_PI * zz *
//   ROOT::Math::ellint_2(zz, k);}; auto q = [k](double zz)
//     {return 1.0 / 3.0 * (7.0 - 8.0 * zz*zz - 4.0 * (1.0 - zz*zz) *
//         ROOT::Math::ellint_1(zz, k) / ROOT::Math::ellint_2(zz, k));};
//   auto a = [](double u){return (u*u + 2.0) / (u * sqrt(u*u + 4.0));};
//   return a(u) * b0(z) * (1.0 + rho*rho/8.0 * q(z));
// };
//
//
// std::vector<double> RTLcYoo_un(std::vector<double> u, Source<double> src) {
//   std::vector<double> mag;
//   for(size_t i=0; i<u.size(); i++) mag.push_back(RTLcYoo_un(u[i],
//   src.circ().r())); return mag;
// };

double RTMagMax(double rho) { return sqrt(rho * rho + 4.0) / rho; };

double RTMagMax(Source<double> src) { return RTMagMax(src.circ().r()); };

double RTLoop::area_tzx() const {
  double area = 0.0;
  double dx, yy;
  for (size_t i = 0; i < c_loop.size(); i++) {
    if (i == c_loop.size() - 1) {
      dx = c_loop[0].x() - c_loop[i].x();
      yy = c_loop[0].y() + c_loop[i].y();
    } else {
      dx = c_loop[i + 1].x() - c_loop[i].x();
      yy = c_loop[i + 1].y() + c_loop[i].y();
    }
    area += dx * yy / 2.0;
  }
  return area;
};

std::array<std::vector<V2d<double>>, 2> RT_img(int nd, Mlens<double> ml,
                                               Source<double> src) {
  double dp = 2.0 * M_PI / nd;
  std::array<std::vector<V2d<double>>, 2> img;
  for (int i = 0; i < nd; i++) {
    V2d<double> edge(src.cent().x() + src.circ().r() * cos(dp * i),
                     src.cent().y() + src.circ().r() * sin(dp * i));
    Source<double> sedg(Circle<double>(edge, 0.01));
    auto si = single_image(ml, sedg);
    img[0].push_back(si[0]);
    img[1].push_back(si[1]);
  }
  return img;
};

// double RTLcWitt(double u, double rho) {
//   double n = 4.0 * u * rho;
//   double k =  sqrt(4.0 * n / ((4.0 + u * rho) * (4.0 + u * rho)));
//   double mu1 = 2.0 / M_PI * (M_PI + )
//   double z = u / rho;
//   double k = std::min(1.0/z, 1.0);
//   auto b0 = [k](double zz){return 4.0 / M_PI * zz *
//   ROOT::Math::ellint_2(k, zz);}; auto q = [k](double zz)
//     {return 1.0 / 3.0 * (7.0 - 8.0 * zz*zz - 4.0 * (1.0 - zz*zz) *
//         ROOT::Math::ellint_1(k, zz) / ROOT::Math::ellint_2(k, zz));};
//   auto a = [](double u){return (u*u + 2.0) / (u * sqrt(u*u + 4.0));};
//   return a(u) * b0(z) * (1.0 + rho*rho/8.0 * q(z));
// };
