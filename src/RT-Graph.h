/**
 * @file RT-Graph.h
 *
 * @brief header file for drawing graphs with CERN root
 * @author Fumio Abe
 *
 */

// define _INC_RT_GRAPH_h_ to load only once
#ifndef _INC_RT_GRAPH_h_
#define _INC_RT_GRAPH_h_

#include "Fractal.h"
#include "Geom.h"
#include "Lcurve.h"
#include <TArrow.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMarker.h>
#include <TPolyLine.h>
#include <TVirtualPad.h>
// #include <Math/SpecFuncMathMore.h>
#include <fstream>

/// @brief struct to express 1d data array
struct Data1d {
  /// @brief x array
  std::vector<double> x;
  /// @brief y array
  std::vector<double> y;
};

/// @brief struct to express 2d data array
struct Data2d {
  /// @brief x array
  std::vector<double> x;
  /// @brief y array
  std::vector<double> y;
  /// @brief z array
  std::vector<double> z;
};

/// @brief Make TGraph from Data1d
/// @param[in] data : 1d data
/// @return TGraph make1dplt(Data1d data) : TGraph to make plot
TGraph make1dplt(Data1d data);
/// @brief Make TGraph2D from Data2d
/// @param[in] data : 2d data
/// @return TGraph2D make2dplt(Data2d data) : TGraph2D to make plot
TGraph2D make2dplt(Data2d data);
/// @brief Make TGraph2D pointer from Data2d
/// @param[in] data : 2d data
/// @return TGraph2D *p_make2dplt(Data2d data) : TGraph2D pointer to make plot
TGraph2D *p_make2dplt(Data2d data);
/// @brief Set draw frame with min and max
/// @param[in] *Pad : pointer to gPad
/// @param[in] min : lower left corner of the draw frame
/// @param[in] max : upper right corner of the draw frame
void SetFrame(TVirtualPad *Pad, dV2d min, dV2d max);
/// @brief Set draw frame with min and max
/// @param[in] *Pad : pointer to gPad
/// @param[in] rec : rectangle to express draw frame
void SetFrame(TVirtualPad *Pad, dRectangle rec);

/**
 * class RTV2d
 * @brief class to draw V2d
 *
 */
class RTV2d {
private:
  /// @brief V2d
  V2d<double> c_v2d = V2d<double>(0.0, 0.0);

public:
  /// @brief Constructor
  /// @param[in] v : 2d vector
  RTV2d(V2d<double> v) : c_v2d(v) {};
  /// @brief Return RMarker
  /// @return TMarker tmarker() : TMaker to plot
  TMarker tmarker() const;
};

/**
 * class RTSeg2d
 * @brief class to draw Seg2d
 *
 */
class RTSeg2d {
private:
  /// @brief 2d segment
  Seg2d<double> c_seg;

public:
  /// @brief Constructor
  /// @param[in] seg : 2d segment
  RTSeg2d(Seg2d<double> seg) : c_seg(seg) {};
  /// @brief Return TLine to draw
  /// @return TLine tline() : TLine
  TLine tline() const;
};

/**
 * class RTTriangle
 * @brief class to draw Triangle
 *
 */
class RTTriangle {
private:
  /// @brief Triangle parameters
  Triangle<double> c_tri;

public:
  /// @brief Constructor
  /// @param[in] tri : Triangle
  RTTriangle(Triangle<double> tri) : c_tri(tri) {};
  /// @brief Return Tline vector to draw triangle
  /// @return std::vector<TLine> tline() : TLine vector to draw triangle
  std::vector<TLine> tline() const;
  /// @brief Return TPolyLine to draw triangle
  /// @return TPolyLine tpolyline() : TPolyLine to draw triangle
  TPolyLine tpolyline() const;
};

/**
 * class RTRectangle
 * @brief class to draw Rectangle
 *
 */
class RTRectangle {
private:
  /// @brief Rectangle parameters
  Rectangle<double> c_rect;

public:
  /// @brief Constructor
  /// @param[in] rec : rectangle parameters
  RTRectangle(Rectangle<double> rec) : c_rect(rec) {};
  /// @brief Return TBox to draw rectangle
  TBox tbox() const;
};

/**
 * class RTCircle
 * @brief class to draw Circle
 *
 */
class RTCircle {
private:
  /// @brief Circle parameters
  Circle<double> c_circ;

public:
  /// @brief Constructor
  /// @param[in] circ : circle parameters
  RTCircle(Circle<double> circ) : c_circ(circ) {};
  /// @brief Return TEllipse to draw circle
  TEllipse tellipse() const;
};

/**
 * class RTEllipse
 * @brief class to draw Ellipse
 *
 */
class RTEllipse {
private:
  /// @brief Ellipse parameters
  Ellipse<double> c_ell;

public:
  /// @brief Constructor
  /// @param[in] ell : ellipse parameters
  RTEllipse(Ellipse<double> ell) : c_ell(ell) {};
  /// @brief Return TEllipse to draw ellipse
  TEllipse tellipse() const;
};

/// @brief struct of parameters to draw line
struct RTLineParam {
  /// @brief Line style
  int style = 1;
  /// @brief Line width
  int width = 1;
  /// @brief Line color
  int color = 1;
  /// @brief Line color transparency
  float trans = 0.7;
  /// @brief Draw option
  std::string draw_opt = "same";
};

/// @brief struct of parameters to fill area
struct RTFillParam {
  /// @brief Fill style
  int style = 1001;
  /// @brief Fill color
  int color = 1;
  /// @brief Fill color transparency
  float trans = 0.3;
  /// @brief Draw option
  std::string draw_opt = "same";
};

struct RTMarkerParam {
  int style = 1;
  float size = 1.0;
  int color = 3;
  float trans = 1.8;
  std::string draw_opt = "same";
};

/// @brief struct of parameters to draw lens
struct RTLensParam {
  /// @brief Number of dots of Einstein ring
  int er_ndiv = 500;
  /// @brief Einstein ring line parameters
  RTLineParam er_line = {1, 1, 3, 0.3, "same"};
  /// @brief Einstein ring fill parameters
  RTFillParam er_fill = {1001, 3, 0.1, "same"};
  /// @brief Lens center marker parameters
  RTMarkerParam marker = {1, 1.0, 3, 0.8, "same"};
};

struct RTCC_SelStep {
  int round = 18;
  double mu = 100.0;
};

struct RTCC_Sel {
  Square<double> sq = Square<double>(V2d<double>(0.0, 0.0), 1.5);
  std::vector<RTCC_SelStep> sel;
};

struct RTCC_Sel2 {
  Square<double> sq;
  int max;
  int sel;
};

/**
 * class RTLens
 * @brief class to draw Lens
 *
 */
class RTLens {
private:
  /// @brief Lens parameters
  Lens<double> c_lens;
  /// @brief Draw parameters;
  RTLensParam c_param;

public:
  /// @brief Constructor
  /// @param[in] lens : lens parameters
  RTLens(Lens<double> lens) : c_lens(lens) {};
  /// @brief Constructor
  /// @param[in] lens : lens parameters
  /// @param[in] param : draw parameters
  RTLens(Lens<double> lens, RTLensParam param)
      : c_lens(lens), c_param(param) {};
  /// @brief Set draw parameters
  /// @param[in] param ; draw parameters
  void SetParam(RTLensParam param) { c_param = param; };
  /// @brief Return draw parameters
  /// @return Draw parameters
  RTLensParam param() const { return c_param; };
  /// @brief Return TMarker to draw center of the lens
  /// @return TMarker marker() : marker of the center
  TMarker marker() const;
  /// @brief Return TEllipse to draw the Einstein ring
  /// @return TEllipse er() : Einstein ring
  TEllipse er() const;
  /// @brief Return TGraph to draw the Einstein ring
  /// @return TEllipse ger() : Einstein ring
  TGraph ger() const;
  RTCC_Sel2 sel() const;
};

struct RTCCCA {
  TGraph cc;
  TGraph ca;
  TGraph lens;
};

/**
 * class RTMlens
 * @brief class to draw Mlens
 *
 */
class RTMlens {
private:
  /// @brief Multiple lens parameters
  Mlens<double> c_ml;

public:
  /// @brief Constructor
  /// @param[in] ml : Multiple lens parameters
  RTMlens(Mlens<double> ml) : c_ml(ml) {};
  /// @brief Return markers of the center of lenses
  /// @return std::vector<TMarker> markers() : markers
  std::vector<TMarker> markers() const;
  TGraph centers() const;
  /// @brief Return points on the Einstein ring of the il'th lens to draw dotted
  /// Einstein ring
  /// @param[in] il : index of the lens
  /// @param[in] ndiv : number of division
  /// @return TGraph er(int il, int ndiv) : Einstein ring of il'th lens
  TGraph er(int il, int ndiv) const;
  /// @brief Return points on the Einstein ring of the il'th lens to draw dotted
  /// Einstein ring
  /// @param[in] il : index of the lens
  /// @return TGraph er(int il) : Einstein ring of il'th lens with default
  /// division (ndiv = 500)
  TGraph er(int il) const;
  /// @brief Return all Einstein rings with ndiv division
  /// @param[in] ndiv : number of division
  /// @return std::vector<TGraph> aller(int ndiv) : all Einstein rings
  std::vector<TGraph>
  aller(int ndiv) const; // Einstein rings (ndiv : number of dots)
  /// @brief Return all Einstein rings
  /// @return std::vector<TGraph> aller() : all Einstein rings
  std::vector<TGraph> aller() const;
  /// @brief Return all Einstein rings
  /// @return std::vector<TGraph> aller() : all Einstein rings
  std::vector<TGraph> allger() const;
  /// @brief Return all Einstein rings
  /// @param[in] param : Draw parameters
  /// @return std::vector<TGraph> aller() : all Einstein rings
  std::vector<TGraph> allger(RTLensParam param) const;
  /// @brief Return all Einstein rings with solid line
  /// @return std::vector<TEllipse> ellipse() : all Einstein rings with solid
  /// line
  std::vector<TEllipse> ellipse() const;
  /// @brief Return 2d data to make contour of the images
  /// @param[in] rec : the area to calculate
  /// @param[in] bt : the source position
  /// @param[in] nx : number of divisions toward x direction
  /// @param[in] ny : number of divisions toward y direction
  /// @param[in] zlim : lower limit of z
  /// @return Data2d dbt() : 2d data of dbt
  Data2d dbt(Rectangle<double> rec, V2d<double> bt, int nx, int ny,
             double zlim) const;
  /// @brief Return 2d data to make contour of the images
  /// @param[in] rec : the area to calculate
  /// @param[in] bt : the source position
  /// @param[in] zlim : lower limit of z
  /// @return Data2d dbt() : 2d data of dbt
  Data2d dbt(Rectangle<double> rec, V2d<double> bt, double zlim) const;
  /// @brief Return 2d data to make contour of the images
  /// @param[in] bt : the source position
  /// @param[in] zlim : lower limit of z
  /// @return Data2d dbt() : 2d data of dbt
  Data2d dbt(V2d<double> bt, double zlim) const;
  /// @brief Return 2d data to make contour of the images
  /// @param[in] bt : the source position
  /// @return Data2d dbt() : 2d data of dbt
  Data2d dbt(V2d<double> bt) const;
  /// @brief Return 2d graph to make contour of the images
  /// @param[in] rec : the area to calculate
  /// @param[in] bt : the source position
  /// @param[in] nx : number of divisions toward x direction
  /// @param[in] ny : number of divisions toward y direction
  /// @param[in] zlim : lower limit of z
  /// @return TGraph2D cont_dbt() : 2d contour data to draw
  TGraph2D cont_dbt(Rectangle<double> rec, V2d<double> bt, int nx, int ny,
                    double zlim) const;
  /// @brief Return 2d graph to make contour of the images
  /// @param[in] rec : the area to calculate
  /// @param[in] bt : the source position
  /// @param[in] zlim : lower limit of z
  /// @return TGraph2D cont_dbt() : 2d contour data to draw
  TGraph2D *p_cont_dbt(Rectangle<double> rec, V2d<double> bt,
                       double zlim) const;
  /// @brief Return 2d graph to make contour of the images
  /// @param[in] bt : the source position
  /// @return TGraph2D cont_dbt() : 2d contour data to draw
  TGraph2D cont_dbt(V2d<double> bt) const;
  /// @brief Return 2d graph pointer (for dynamic memory) to make contour of the
  /// images
  /// @param[in] bt : the source position
  /// @param[in] zlim : lower limit of z
  /// @return TGraph2D *p_cont_dbt() : 2d contour data pointer to draw
  TGraph2D *p_cont_dbt(V2d<double> bt, double zlim) const;
  /// @brief Return 2d graph pointer (for dynamic memory) to make contour of the
  /// images
  /// @param[in] bt : the source position
  /// @return TGraph2D *p_cont_dbt() : 2d contour data pointer to draw
  TGraph2D *p_cont_dbt(V2d<double> bt) const;
  /// @brief Return graph to make scatter plot of the images
  /// @param[in] rec : the area to calculate
  /// @param[in] bt : the source position
  /// @param[in] nx : number of divisions toward x direction
  /// @param[in] ny : number of divisions toward y direction
  /// @param[in] rho : source radious divided by the Einstein radious
  /// @return TGraph scat_dbt() : scatter plot data to draw
  TGraph scat_dbt(Rectangle<double> rec, V2d<double> bt, int nx, int ny,
                  double rho) const;
  /// @brief Return graph to make scatter plot of the images
  /// @param[in] rec : the area to calculate
  /// @param[in] bt : the source position
  /// @param[in] rho : source radious divided by the Einstein radious
  /// @return TGraph scat_dbt() : scatter plot data to draw
  TGraph scat_dbt(Rectangle<double> rec, V2d<double> bt, double rho) const;
  /// @brief Return graph to make scatter plot of the images
  /// @param[in] bt : the source position
  /// @param[in] rho : source radious divided by the Einstein radious
  /// @return TGraph scat_dbt() : scatter plot data to draw
  TGraph scat_dbt(V2d<double> bt, double rho) const;
  /// @brief Return GNUPLOT data to make scatter plot of the images
  /// @param[in] rec : the area to calculate
  /// @param[in] bt : the source position
  /// @param[in] nx : number of divisions toward x direction
  /// @param[in] ny : number of divisions toward y direction
  /// @return std::string gp_dbt() : scatter plot data for GNUPLOT
  std::string gp_dbt(Rectangle<double> rec, V2d<double> bt, int nx,
                     int ny) const;
  /// @brief Return GNUPLOT data to make scatter plot of the images
  /// @param[in] rec : the area to calculate
  /// @param[in] bt : the source position
  /// @return std::string gp_dbt() : scatter plot data for GNUPLOT
  std::string gp_dbt(Rectangle<double> rec, V2d<double> bt) const;
  /// @brief Return GNUPLOT data to make scatter plot of the images around a
  /// lens
  /// @param[in] lens : the lens to draw
  /// @param[in] bt : the source position
  /// @param[in] nx : number of divisions toward x direction
  /// @param[in] ny : number of divisions toward y direction
  /// @return std::string gp_dbt() : scatter plot data for GNUPLOT
  std::string gp_dbt(Lens<double> lens, V2d<double> bt, int nx, int ny) const;
  /// @brief Return GNUPLOT data to make scatter plot of the images around a
  /// lens
  /// @param[in] lens : the lens to draw
  /// @param[in] bt : the source position
  /// @return std::string gp_dbt() : scatter plot data for GNUPLOT
  std::string gp_dbt(Lens<double> lens, V2d<double> bt) const;
  /// @brief Return GNUPLOT data to make scatter plot of the images around a
  /// lens
  /// @param[in] bt : the source position
  /// @return std::string gp_dbt() : scatter plot data for GNUPLOT
  std::string gp_dbt(V2d<double> bt) const;
  /// @brief Save GNUPLOT data to a file
  /// @param[in] file : file name to save
  /// @param[in] data : data to output
  void save_gp(std::string file, std::string &data) const;
  /// @brief Return TMQue<double> corresponding to critical curve
  /// @param[in] sel : selection parameters
  /// @return TMQue<double> ccca(std::vector<RTCC_Sel> sel) : TMQue<double>
  /// corresponding to critical curves
  TMQue<double> ccca(std::vector<RTCC_Sel> sel) const;
  /// @brief Return Tgraphs of critical curves and caustics
  /// @param[in] sel : selection parameters
  /// @return RTCCCA ccca(std::vector<RTCC_Sel> sel) : TMQue<double>
  /// corresponding to critical curves
  RTCCCA grccca(std::vector<RTCC_Sel> sel) const;
  std::vector<RTCC_Sel2> sel2() const;
  std::vector<TMQue<double>> ccca2_1(int n, std::vector<RTCC_Sel2> sel2) const;
  TMQue<double> ccca2(int n, std::vector<RTCC_Sel2> sel2) const;
  RTCCCA grccca2(int n, std::vector<RTCC_Sel2> sel2) const;
  RTCCCA grccca2(std::vector<RTCC_Sel2> sel2) const;
  RTCCCA grccca2() const;
};

/// @brief struct of parameters to draw source
struct RTSourceParam {
  /// @brief Center marker parameters
  RTMarkerParam marker = {3, 1, 2, 0.8, "same"};
  /// @brief Shape line parameters
  RTLineParam ell_line = {1, 1, 2, 0.8, "same"};
  /// @brief Shape fill parameters
  RTFillParam ell_fill = {1001, 2, 0.3, "same"};
  /// @brief Surface brightness line parameters
  RTLineParam sb_line = {1, 1, 2, 0.8, "same"};
  /// @brief Surface brightness fill parameters
  RTFillParam sb_fill = {1001, 2, 0.3, "same"};
};

/**
 * class RTSource
 * @brief class to draw Source
 *
 */
class RTSource {
private:
  /// @brief Source parameters
  Source<double> c_src;
  /// @brief RTSource parameters
  RTSourceParam c_param;

public:
  /// @brief Constructor
  /// @param[in] src : source parameters
  RTSource(Source<double> src) : c_src(src) {};
  /// @brief Constructor
  /// @param[in] src : source parameters
  /// @param[in] param : draw parameters
  RTSource(Source<double> src, RTSourceParam param)
      : c_src(src), c_param(param) {};
  /// @brief Return marker of the center
  /// @return TMarker mark() : marker
  TMarker mark() const;
  /// @brief Return the shape
  /// @return TEllipse ell() : shape
  TEllipse ell() const;
  /// @brief Return the surface brightness
  /// @param[in] n : number of divisions
  /// @param[in] min : lower left corner of the graph
  /// @param[in] max : upper right corner of the graph
  /// @return TGraph2D sb(int n, dV2d min, dV2d max) : surface brightness
  TGraph2D sb(int n, dV2d min, dV2d max);
  /// @brief Return the surface brightness
  /// @param[in] min : lower left corner of the graph
  /// @param[in] max : upper right corner of the graph
  /// @return TGraph2D sb(dV2d min, dV2d max) : surface brightness
  TGraph2D sb(dV2d min, dV2d max);
  /// @brief Return the surface brightness
  /// @param[in] fact : factor to enlarge area
  /// @return TGraph2D sb(double fact) : surface brightness
  TGraph2D sb(double fact);
  /// @brief Return the surface brightness
  /// @return TGraph2D sb() : surface brightness
  TGraph2D sb();
};

class RTSrcMotion {
private:
  /// @brief SrcMotion
  SrcMotion<double> c_smot;

public:
  /// @brief Constructor for RTSrcMotion
  /// @param[in] smot : SrcMotion
  RTSrcMotion(SrcMotion<double> smot) : c_smot(smot) {};
  /// @brief  Returns Tmaker of source
  /// @param[i] t : time
  TMarker marker(Time<double> t) const;
  /// @brief Returns Tgraph for source trajectory
  /// @param[in] t : a series of time
  TGraph graph(std::vector<Time<double>> t) const;
  /// @brief Returns a series of markers
  /// @param[in] t : a series of time
  std::vector<TMarker> markers(std::vector<Time<double>> t) const;
};

class RTLensMotion {
private:
  LensMotion<double> c_lm;

public:
  RTLensMotion(LensMotion<double> lm) : c_lm(lm) {};
  TMarker marker(Time<double> t) const;
  TGraph graph(std::vector<Time<double>> t) const;
  std::vector<TMarker> markers(std::vector<Time<double>> t) const;
};

class RTMlMotion {
private:
  MlMotion<double> c_mlm;

public:
  RTMlMotion(MlMotion<double> mlm) : c_mlm(mlm) {};
  std::vector<TMarker> markers(Time<double> t) const;
  std::vector<TGraph> graphs(std::vector<Time<double>> t) const;
  std::vector<std::vector<TMarker>> markers(std::vector<Time<double>> t) const;
};

/**
 * class RTTriL
 * @brief class to draw lens triangle
 *
 */
class RTTriL {
private:
  /// @brief Triangle parameters
  TriL<double> c_tril;

public:
  /// @brief Constructor
  /// @param[in] tril : lens triangle parameters
  RTTriL(TriL<double> tril) : c_tril(tril) {};
  /// @brief Return shape of the triangle
  /// @return TPolyLine tpolyline() : shape of the triangle
  TPolyLine tpolyline() const;
};

/**
 * class RTTriS
 * @brief class to draw source triangle
 *
 */
class RTTriS {
private:
  /// @brief Triangle parameters
  TriS<double> c_tris;

public:
  /// @brief Constructor
  /// @param[in] tris : source triangle parameters
  RTTriS(TriS<double> tris) : c_tris(tris) {};
  /// @brief Return shape of the triangle
  /// @return TPolyLine tpolyline() : shape of the triangle
  TPolyLine tpolyline() const; // Sape of TriS
};

/**
 * class RTTriMap
 * @brief class to draw TriMap (TriL and TriS)
 *
 */
class RTTriMap {
private:
  /// @brief TriMap parameters
  TriMap<double> c_trimap;

public:
  /// @brief Constructor
  /// @param[in] trimap : TriMap parameters (TriL and TriS)
  RTTriMap(TriMap<double> trimap) : c_trimap(trimap) {};
  /// @brief Return shape of the lens triangle
  /// @return  : lens triangle to draw
  TPolyLine lens() const;
  /// @brief Return shape of the source triangle
  /// @return  : source triangle to draw
  TPolyLine src() const;
};

/**
 * class RTTMQue
 * @brief class to draw TriMap queue
 *
 */
class RTTMQue {
private:
  /// @brief TriMap queue
  TMQue<double> c_tmq;

public:
  /// @brief Constructor
  /// @param[in] tmq : TriMap queue
  RTTMQue(TMQue<double> tmq) : c_tmq(tmq) {};
  /// @brief Return shapes of TriL's
  /// @return std::vector<TPolyLine> lens() : shapes of the TriL
  std::vector<TPolyLine> lens() const;
  /// @brief Return shapes of TriS's
  /// @return std::vector<TPolyLine> src() : shapes of the TriS
  std::vector<TPolyLine> src() const;
  /// @brief Return center points on lens plane
  /// @return TGraph gr_lens() : center points on lens plane
  TGraph gr_lens() const;
  /// @brief Return center points on source plane
  /// @return TGraph gr_src() : center points on source plane
  TGraph gr_src() const;
};

/// @brief struct of parameters to make TCanvas
struct RTCanvas {
  /// @brief Pointer to TCanvas
  TCanvas *c;
  /// @brief Name of TCanvas
  std::string name = "Plot";
  /// @brieff Title of TCanvas
  std::string title = "Plot";
  /// @brief Width of TCanvas
  int width = 800;
  /// @brief Height of TCanvas
  int height = 800;
  /// @brief Rectangle corredponding to DrawFrame
  // dRectangle frame = dRectangle(dV2d(-1.4, -1.5), dV2d(1.6, 1.5));
  dRectangle frame;
  /// @brief File name to save
  std::string file = "plot.png";
};

/// @brief Make TCanvas
/// @param[in] &rtc : TCanvas parameters
void RTMkCanvas(RTCanvas &rtc);

/// @brief Save TCanvas to the file
/// @param[in] &rtc : TCanvas parameters
void RTSvCanvas(RTCanvas &rtc);

/// @brief Close TCanvas
/// @param[in] &rtc : TCanvas parameters
void RTClCanvas(RTCanvas &rtc);

/// @brief struct of parameters to draw FrTri
struct RTFrTriParam {
  /// @brief Name of the plot
  std::string name = "FrTri plot";
  /// @brief Title of the plot
  std::string title = "FrTri plot";
  /// @brief Fine name to save
  std::string file = "frtri-plot.png";
  /// @brief Parameters to draw lenses
  RTLensParam lens;
  /// @brief Parameters to draw source
  RTSourceParam src;
  /// @brief Parameters to draw lines of "in" triangle on the lens plane
  RTLineParam in_lens_line = {1, 1, 2, 0.7, "same"};
  /// @brief Parameters to fill "in" triangle on the lens plane
  RTFillParam in_lens_fill = {1001, 2, 0.3, "same"};
  /// @brief Parameters to draw lines of "out" triangle on the lens plane
  RTLineParam out_lens_line = {1, 1, 5, 0.7, "same"};
  /// @brief Parameters to fill "out" triangle on the lens plane
  RTFillParam out_lens_fill = {1001, 5, 0.3, "same"};
  /// @brief Parameters to draw lines of "unc" triangle on the lens plane
  RTLineParam unc_lens_line = {1, 1, 7, 0.7, "same"};
  /// @brief Parameters to fill "unc" triangle on the lens plane
  RTFillParam unc_lens_fill = {1001, 7, 0.3, "same"};
  /// @brief Parameters to draw lines of "in" triangle on the source plane
  RTLineParam in_src_line = {1, 1, 2, 0.7, "same"};
  /// @brief Parameters to fill "in" triangle on the source plane
  RTFillParam in_src_fill = {1001, 2, 0.3, "same"};
  /// @brief Parameters to draw lines of "out" triangle on the source plane
  RTLineParam out_src_line = {1, 1, 5, 0.7, "same"};
  /// @brief Parameters to fill "out" triangle on the source plane
  RTFillParam out_src_fill = {1001, 5, 0.3, "same"};
  /// @brief Parameters to draw lines of "unc" triangle on the source plane
  RTLineParam unc_src_line = {1, 1, 7, 0.7, "same"};
  /// @brief Parameters to fill "unc" triangle on the source plane
  RTFillParam unc_src_fill = {1001, 7, 0.3, "same"};
};

/// @brief struct of object to draw on the lens plane
struct RTFTLensObj {
  /// @brief Lens markers
  std::vector<TMarker> lens_mk;
  /// @brief Lens Einstein rings
  std::vector<TGraph> lens_er;
  /// @brief Source marker
  TMarker src_mk;
  /// @brief "in" Triangles
  std::vector<TPolyLine> in_lens;
  /// @brief "out" Triangles
  std::vector<TPolyLine> out_lens;
  /// @brief "unc" Triangles
  std::vector<TPolyLine> unc_lens;
};

/// @brief struct of object to draw on the source plane
struct RTFTSrcObj {
  /// @brief Source marker
  TMarker src_mk;
  /// @brief Source shape
  TEllipse src_ell;
  /// @brief Source surface brightness
  TGraph2D src_sb;
  /// @brief "in" Triangles
  std::vector<TPolyLine> in_src;
  /// @brief "out" Triangles
  std::vector<TPolyLine> out_src;
  /// @brief "unc" Triangles
  std::vector<TPolyLine> unc_src;
};

struct RTFTRoff {
  std::vector<int> round;
  std::vector<float> area_f;
  std::vector<double> area_d;
  std::vector<long double> area_l;
  std::vector<float> area_sn_f;
  std::vector<double> area_sn_d;
  std::vector<long double> area_sn_l;
  TGraph ar_fl;
  TGraph ar_dl;
  TGraph ar_sn_fl;
  TGraph ar_sn_dl;
};

/**
 * class RTFrTri
 * @brief class to draw FrTri (Fractal algorithm)
 *
 */
class RTFrTri {
private:
  /// @brief Fractal parameters
  FrTri<double> c_frtri;
  RTFrTriParam c_param;

public:
  /// @brief Constructor
  /// @param[in] frtri : fractal parameters
  RTFrTri(FrTri<double> frtri) : c_frtri(frtri) {};
  /// @brief Constructor
  /// @param[in] frtri : fractal parameters
  /// @param[in] param : Draw parameters
  RTFrTri(FrTri<double> frtri, RTFrTriParam param)
      : c_frtri(frtri), c_param(param) {};
  /// @brief Return FrTriParam
  /// @return FrTriParam
  RTFrTriParam param() const { return c_param; };
  /// @brief Return in queue
  /// @return RTTMQue in() : "in" queue
  RTTMQue in() const;
  /// @brief Return uncertain queue
  /// @return RTTMQue uncertain() : "uncertain" queue
  RTTMQue uncertain() const;
  /// @brief Return out queue
  /// @return RTTMQue out() : "out" queue
  RTTMQue out() const;
  void process();
  void process(int n);
  /// @brief Return in triangles on the lens plane
  /// @return std::vector<TPolyLine> in_lens() : in triangles on the lens plane
  std::vector<TPolyLine> in_lens() const;
  /// @brief Return out triangles on the lens plane
  /// @return std::vector<TPolyLine> out_lens() : out triangles on the lens
  /// plane
  std::vector<TPolyLine> out_lens() const;
  /// @brief Return uncertain triangles on the lens plane
  /// @return std::vector<TPolyLine> unc_lens() : uncertain triangles on the
  /// lens plane
  std::vector<TPolyLine> unc_lens() const;
  /// @brief Return in triangles on the source plane
  /// @return std::vector<TPolyLine> in_lens() : in triangles on the source
  /// plane
  std::vector<TPolyLine> in_src() const;
  /// @brief Return out triangles on the source plane
  /// @return std::vector<TPolyLine> out_lens() : out triangles on the source
  /// plane
  std::vector<TPolyLine> out_src() const;
  /// @brief Return uncertain triangles on the source plane
  /// @return std::vector<TPolyLine> unc_lens() : uncertain triangles on the
  /// source plane
  std::vector<TPolyLine> unc_src() const;
  /// @brief Return FTLensObj
  /// @return FTLensObj
  RTFTLensObj mk_lensobj() const;
  std::vector<RTFTLensObj> mk_lensobj(int round);
  /// @brief Return FTSrcObj
  /// @return FTSrcObj
  RTFTSrcObj mk_srcobj() const;
  std::vector<RTFTSrcObj> mk_srcobj(int round);
  /// @brief Set drawing parameters of FTLensObj
  /// @param[in] lo : FTLensObj
  void set_draw_lens(RTFTLensObj &lo) const;

  /// @brief Set drawing parameters of FTSrcObj
  /// @param[in] so : FTSrcObj
  void set_draw_src(RTFTSrcObj &so) const;
  /// @brief Draw FTLensObj
  /// @param[in] c : Canvas
  /// @param[in] lo : FTLensObj
  void draw_lens(RTCanvas &c, RTFTLensObj &lo) const;
  void draw_lens(RTCanvas &c, RTFTLensObj &lo, RTCCCA &ccca) const;
  void mk_anim_lens(RTCanvas &c, std::vector<RTFTLensObj> &vlo,
                    std::string fprefix) const;
  void mk_anim_lens(RTCanvas &c, std::vector<RTFTLensObj> &vlo, RTCCCA &ccca,
                    std::string fprefix) const;
  /// @brief Draw FTSrcObj
  /// @param[in] c : Canvas
  /// @param[in] so : FTSrcObj
  void draw_src(RTCanvas &c, RTFTSrcObj &so) const;
  RTFTRoff rtftroff(int rdmin, int rdmax) const;
};

void RTDraw_images(RTCanvas &rtc, FrTri<double> &frt,
                   std::vector<Square<double>> &grp, std::string pre,
                   std::string type);
void RTDraw_images(RTCanvas &rtc, FrTri<double> frt,
                   std::vector<Square<double>> grp);

void RTDraw_src_tri(RTCanvas &rtc, FrTri<double> &frt,
                    std::vector<Square<double>> &grp, Square<double> sqr,
                    std::string pre, std::string type);

class RTFrPrec {
private:
  std::vector<double> c_round;
  std::vector<double> c_d;
  std::vector<double> c_area_in;
  std::vector<double> c_area_unc;
  std::vector<double> c_area_sn;
  std::vector<double> c_size_in;
  std::vector<double> c_size_unc;

public:
  RTFrPrec() {};
  void add(FrTri<double> frtri);
  std::vector<double> round() const;
  std::vector<double> d() const;
  std::vector<double> area_in() const;
  std::vector<double> area_unc() const;
  std::vector<double> area_cent() const;
  std::vector<double> area_sn() const;
  std::vector<double> size_in() const;
  std::vector<double> size_unc() const;
  std::vector<double> dar() const;
  std::vector<double> rdar() const;
  std::vector<double> rer() const;
  std::vector<double> rersn() const;
  TGraph d_vs_rd() const;
  TGraph arc_vs_rd() const;
  TGraph arsn_vs_rd() const;
  TGraph szi_vs_rd() const;
  TGraph szu_vs_rd() const;
  TGraph dar_vs_rd() const;
  TGraph rdar_vs_rd() const;
  TGraph rer_vs_rd() const;
  TGraph rersn_vs_rd() const;
};

class RTFTBench {
private:
  FrTri<double> c_dft;
  double c_ctime;
  int c_in_size;
  int c_unc_size;

public:
  RTFTBench(FrTri<double> dft) : c_dft(dft) {};
};

class RTLcParam {
private:
  LcParam<double> c_lcp;

public:
  RTLcParam(LcParam<double> lcp) : c_lcp(lcp) {};
  LcParam<double> lcp() const { return c_lcp; };
  TGraph in_sz_vs_rd() const;
  TGraph unc_sz_vs_rd() const;
  TGraph out_sz_vs_rd() const;
};

/**
 * class RTLcurve
 * @brief class to draw Lcurve
 *
 */
class RTLcurve {
private:
  /// @brief Lcurve parameters
  Lcurve<double> c_lc;

public:
  /// @brief Constructor
  /// @param[in] lc : Lcurve parameters
  RTLcurve(Lcurve<double> &lc) : c_lc(lc) {};
  /// @brief Return source motions
  /// @return std::vector<SrcMotion<double>> srcms() : source motions
  std::vector<SrcMotion<double>> srcms() const;
  /// @brief Return source motion of the first one
  /// @return SrcMotion<double> srcm() : motion of the first source
  SrcMotion<double> srcm() const;
  Source<double> src(int it) const;
  Source<double> src(Time<double> t) const;
  FrTri<double> frtri(int it) const;
  FrTri<double> frtri(Time<double> t) const;
  RTFrTri rtfrtri(int it) const;
  RTFrTri rtfrtri(Time<double> t) const;
  RTSource rtsrc(int it) const;
  RTSource rtsrc(Time<double> t) const;
  RTMlens rtml(int it) const;
  RTMlens rtml(Time<double> t) const;
  RTCCCA rtccca2(int it) const;
  RTCCCA rtccca2(Time<double> t) const;
  TGraph2D *p_cont_dbt(int it) const;
  TGraph2D *p_cont_dbt(Time<double> t) const;
  RTFrPrec prec(int it) const;
  RTFrPrec prec(Time<double> t) const;
  RTFTLensObj lensobj(int it, int sel) const;
  RTFTLensObj lensobj(Time<double> t, int sel) const;
  RTFTLensObj lensobj_w(int it, int in_sise) const;
  RTFTLensObj lensobj_w(Time<double> t, int in_size) const;
  RTFTSrcObj srcobj(int it, int sel) const;
  RTFTSrcObj srcobj(Time<double> t, int sel) const;
  void set_rtc(RTCanvas &rtc, int it, int il, double fact) const;
  void draw_lens(RTCanvas &c, RTFTLensObj &lo) const;
  void draw_src(RTCanvas &c, RTFTSrcObj &so) const;
  std::vector<RTFTLensObj> mk_lensobj(int it, int round) const;
  std::vector<RTFTLensObj> mk_lensobj(Time<double> t, int round) const;
  std::vector<RTFTLensObj> mk_lensobj_s(int round) const;
  std::vector<RTFTLensObj> mk_lensobj_w(int in_size) const;
  std::vector<RTFTSrcObj> mk_srcobj(int it, int round) const;
  std::vector<RTFTSrcObj> mk_srcobj(Time<double> t, int round) const;
  void mk_anim_lens(RTCanvas &c, std::vector<RTFTLensObj> &lo,
                    std::string fprefix) const;
  void mk_anim_lens(RTCanvas &c, std::vector<RTFTLensObj> &lo, RTCCCA &ccca,
                    std::string fprefix) const;
  void show_anim_lens(RTCanvas &c, std::vector<RTFTLensObj> &lo,
                      int time) const;
  void show_anim_lens(RTCanvas &c, std::vector<RTFTLensObj> lo) const;

  /// @brief Return graph of the first source trajectory
  /// @return TGraph srctrj() : graph of the first source trajectory
  TGraph srctrj() const;
  /// @brief Return MlMotion<double>
  /// @return MlMotion<double> mlm() : MlMotion
  MlMotion<double> mlm() const;
  /// @brief Return Mlens<double> at time t
  /// @param[in] t : time
  /// @return Mlens<double> ml(Time<double> t) : Mlens<double> at time t
  Mlens<double> ml(Time<double> t) const;
  Mlens<double> ml(int it) const;
  /// @brief Return Mlens<double> at t[0]
  /// @return  : Mlens<double> at time t[0]
  Mlens<double> ml() const;
  /// @brief Return graph of centers of lenses at time t
  /// @param[in] t : time
  /// @return TGraph gml(Time<double> t) : center of lenses
  TGraph gml(Time<double> t) const;
  /// @brief Return graph of centers of lenses at time t[0]
  /// @return TGraph gml() : center of lenses at time t[0]
  TGraph gml() const;
  /// @brief Return critical curves and caustics at time t
  /// @param t : time
  /// @param[in] sel : selection parameters
  /// @return RTCCCA ccca() : graphs of critical curve
  /// curves and causitics
  RTCCCA ccca(Time<double> t, std::vector<RTCC_Sel> sel) const;
  /// @brief Return critical curves and caustics at time t[0]
  /// @param[in] sel : selection parameters
  /// @return RTCCCA ccca() : graphs of critical curve
  /// curves and causitics
  RTCCCA ccca(std::vector<RTCC_Sel> sel) const;
  RTCCCA ccca() const;
  /// @brief Return trajectories of lenese
  /// @return std::vector<TGraph> lenstrj() : trajectories
  std::vector<TGraph> lenstrj() const;
  /// @brief Return magnification graph array
  /// @return std::vector<TGraph> mag() : array of magnification graphs
  std::vector<TGraph> mags() const;
  // std::vector<TGraph> mags2() const;
  /// @brief Return magnification graph for first source
  /// @return TGraph mag() : magnification graph
  TGraph mag() const;
  // TGraph mag() const;
  /// @brief Retuen graph of astrometric shift
  /// @return TGraph as() : astrometric shift
  TGraph as() const;
  // TGraph as() const;
};

/// brief Returns maximum magnification of the single lens event
/// @param[in] rho : source radius
double RTMagMax(double rho);

/// brief Returns maximum magnification of the single lens event
/// @param[in] src : source
double RTMagMax(Source<double> src);

/**
 * class RTLoop
 * @brief class to calculate area of loop
 *
 */
class RTLoop {
private:
  std::vector<V2d<double>> c_loop;

public:
  /// @brief Constructor
  /// @param[in] loop : a loop
  RTLoop(std::vector<V2d<double>> loop) : c_loop(loop) {};
  /// @brief Returns area of the loop
  /// @return double : integrated area
  double area_tzx() const;
};

std::array<std::vector<V2d<double>>, 2> RT_img(int nd, Mlens<double> ml,
                                               Source<double> src);

#endif
