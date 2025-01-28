/**
 * @file Fractal.h
 *
 * @brief header file for functions, classes to handle fractal algorithm
 * @author Fumio Abe
 *
 */

// define _INC_Fractal_h_ to load only once
#ifndef _INC_Fractal_h_
#define _INC_Fractal_h_

// #include "Source.h"
// #include "Source.h"
#include "Motion.h"
#include <ctime>
#include <fstream>

/// @brief enum Zone : in, uncertain, out
enum class Zone {
  /// @brief This object is securely inside
  in,
  /// @brief This object is uncertain
  uncertain,
  /// @brief This object is securely outside
  out
};

std::string zmap(Zone zone);

/**
 * template class TriL
 * @brief Triangle on the lens plane
 *
 * Type : float,  double,  long double
 * Aliases : fTriL (TriL<float),
 *           dTriL (TriL<double),
 *           lTriL (TriL<long double)
 *
 */
template <class Type> class TriL : public Triangle<Type> {
public:
  /// @brief Constructor : 	isosceles right triangle
  /// @param[in] apex0 :  right angle apex
  /// @param[in] apex1 :  left apex
  TriL(V2d<Type> apex0, V2d<Type> apex1)
      : Triangle<Type>(
            {apex0, apex1, apex0 + ((apex1 - apex0) ^ (0.5 * M_PI))}) {};
  /// @brief Return right angle apex
  /// @return V2d<Type> TriL<Type>::right() : right angle apex
  V2d<Type> right() const;
  /// @brief Return three apexes
  /// @return std::array<V2d<Type>, 3> TriL<Type>::apex() : apexes
  std::array<V2d<Type>, 3> apex() const;
  /// @brief Return length of the short segment
  /// @return Type TriL<Type>::d() : alength of the short segment
  Type d() const;
  /// @brief Return if the triangle is too large compared to the circle or not
  /// @param[in]  cl : the circle to compare
  /// @param[in]  &drl :  the mergin factor
  /// @return bool TriL<Type>::too_large() : the triangle is too large compared
  /// to the circle or not
  bool too_large(Circle<Type> cl, Type &drl) const;
  /// @brief Return area of the triangle
  /// @return Type TriL<Type>::area() : area of the triangle
  Type area() const;
  /// @brief Return zone of the triangle
  /// @param[in] cl : the circle to compare
  /// @return Zone TriL<Type>::zone() :  zone in the circle
  Zone zone(Circle<Type> cl) const;
  /// @brief Return zone of the triangle  (mergin : mrg)
  /// @param[in] cl : the circle to compare
  /// @param[in] mrg : the mergin to judge
  /// @return Zone TriL<Type>::zone() :  zone in the circle (mergin : mrg)
  Zone zone(Circle<Type> cl, Type mrg) const;
  /// @brief Return zone of the triangle  (mergin : mrg)
  /// @param[in] el : the ellipse to compare
  /// @param[in] mrg : the mergin to judge
  /// @return Zone TriL<Type>::zone() :  zone in the circle (mergin : mrg)
  Zone zone(Ellipse<Type> el, Type mrg) const;
  /// @brief Return if the triangle on the border on the circle or not
  /// @param[in] cl : the circle to compare
  /// @return bool TriL<Type>::border() : on the border (true) or not (false)
  bool border(Circle<Type> cl) const;
  /// @brief Return if the triangle is adjoint to another one (tril) or not
  /// @param[in] tril : another triangle to compare
  /// @return bool TriL<Type>::adjoin() : adjoin with tril (true) or not
  /// (false)
  bool adjoin(TriL<Type> tril) const;
  /// @brief Split this triangle to self similar two triangles
  /// @return std::array<TriL<Type>, 2> TriL<Type>::split() : self similar
  /// triangles
  std::array<TriL<Type>, 2> split() const;
  /// @brief returns converted angle to TriL<float> (32-bit floating point)
  /// @return TriL<float> : converted angle to TriL<float>
  TriL<float> to_f() const;
  /// @brief returns converted angle to TriL<double> (64-bit floating point)
  /// @return TriL<double> : converted angle to TriL<double>
  TriL<double> to_d() const;
  /// @brief returns converted angle to TriL<long double> (128-bit floating
  /// point)
  /// @return TriL<long double> : converted angle to TriL<long double>
  TriL<long double> to_l() const;
};

/// @brief Global const : null TriL
template <class Type>
const TriL<Type> TriL_0(TriL<Type>(V2d_0<Type>, V2d_0<Type>));

/// @brief struct for brightness
template <typename Type> struct Brightness {
  /// @brief mean value
  Type mean = 0.0;
  /// @brief Root-Means-Squares
  Type rms = 0.0;
};

/**
 * template class TriS
 * @brief Triangle on the source plane
 * Type : float,  double,  long double
 *
 */
template <class Type> class TriS : public Triangle<Type> {
public:
  /// @brief constructor TriS<Type>::TriS(TriL<Type> &tril, Mlens<Type> &ml)
  /// @param[in] tril : triangle on the lensing plane
  /// @param[in] &ml : multiple lens system
  TriS(TriL<Type> tril, Mlens<Type> ml)
      : Triangle<Type>({ml.beta(tril.apex()[0]), ml.beta(tril.apex()[1]),
                        ml.beta(tril.apex()[2])}) {};
  /// @brief Constructor to make divided triangle from TriL
  /// @param[in] tris : "original" TriS
  /// @param[in]  &right : the the right angle apex at lens plane
  /// @param[in]  &ml : the lens system
  TriS(TriS<Type> tris, V2d<Type> &right, Mlens<Type> ml)
      : Triangle<Type>({ml.beta(right), tris.apex()[1], tris.apex()[2]}) {};
  /// @brief Constructor from all apexes
  /// @param[in] aps : all apexes
  TriS(std::array<V2d<Type>, 3> aps) : Triangle<Type>(aps) {};
  /// @brief Returns apexes
  /// @return std::array<V2d<Type>, 3> apex() : apex array
  std::array<V2d<Type>, 3> apex() const;
  /// @brief Retuens length of maximum side
  /// @return Type maxd() : length of maximum side
  Type maxd() const;
  /// @brief Returns Zone relative to the circle
  /// @param[in] cl : Circle
  /// @return Zone zone(Circle<Type> cl) : Zone
  Zone zone(Circle<Type> cl) const;
  /// @brief Returns Zone relative to the circle with a margin
  /// @param[in] cl : Circle
  /// @param[in] mrg : margin
  /// @return Zone zone() : Zone
  Zone zone(Circle<Type> cl, Type mrg) const;
  /// @brief Returns Zone relative to the ellipse
  /// @param[in] el : ellipse
  /// @return Zone zone(Ellipse<Type> el) : Zone
  Zone zone(Ellipse<Type> el) const;
  /// @brief Returns Zone relative to the ellipse with a margin
  /// @param[in] el : ellipse
  /// @param[in] mrg : margin
  /// @return Zone zone() : Zone
  Zone zone(Ellipse<Type> el, Type mrg) const;
  /// @brief Returns whether the triangle is on border of the circle or not
  /// @param cl : circle
  /// @return bool border() : true if on border
  bool border(Circle<Type> cl) const;
  /// @brief Returns whether the triangle is on border of the ellipse or not
  /// @param el : ellipse
  /// @return bool border() : true if on border
  bool border(Ellipse<Type> el) const;
  /// @brief Return the surface brightnesses on the each apexes
  /// @param[in] src : the source
  /// @return std::vector<Type> TriS<Type>::sb() : the surface brightnesses
  std::vector<Type> sb(Source<Type> src) const;
  /// @brief Returns the brightnesses of the triangle
  /// @param[in] src : the source
  /// @return Brightness<Type> TriS<Type>::brightness() : total brightness of
  /// the triangle
  Brightness<Type> brightness(Source<Type> src) const;
  /// @brief Returns the splitted triangles
  /// @param[in] right : the mapped right angle corner
  /// @return std::array<TriS<Type>, 2> TriS<Type>::split() : splitted triangles
  std::array<TriS<Type>, 2> split(V2d<Type> right) const;
  /// @brief returns converted angle to TriS<float> (32-bit floating point)
  /// @return TriS<float> : converted angle to TriS<float>
  TriS<float> to_f() const;
  /// @brief returns converted angle to TriS<double> (64-bit floating point)
  /// @return TriS<double> : converted angle to TriS<double>
  TriS<double> to_d() const;
  /// @brief returns converted angle to TriS<long double> (128-bit floating
  /// point)
  /// @return TriS<long double> : converted angle to TriS<long double>
  TriS<long double> to_l() const;
};

/// @brief Global null const TriS
template <class Type>
const TriS<Type> TriS_0 = TriS<Type>(TriL_0<Type>, Mlens_0<Type>);

/**
 * template class TriMap
 * @brief Triangle on the lens plane and the mapped triangle on the source plane
 * Type : float,  double,  long double
 *
 */
template <class Type> class TriMap {
private:
  /// @brief triangle on the lens plane
  TriL<Type> c_tril;
  /// @brief mapped triangle on the source plane
  TriS<Type> c_tris;

public:
  /// @brief Blank constructor
  TriMap() : c_tril(TriL_0<Type>), c_tris(TriS_0<Type>) {};
  /// @brief Constructor
  /// @param[in] &tril : triangle on the lens plane
  /// @param[in] &ml : multiple lens system
  TriMap(TriL<Type> &tril, Mlens<Type> &ml)
      : c_tril(tril), c_tris(TriS(c_tril, ml)) {};
  /// @brief Constructor
  /// @param[in] &tril : triangle on the lens plane
  /// @param[in] tris : triangle on the source plane
  TriMap(TriL<Type> tril, TriS<Type> tris) : c_tril(tril), c_tris(tris) {};
  /// @brief Return tril
  /// @return TriL<Type> tril() : triangle on the lens plane
  TriL<Type> tril() const;
  /// @brief Return tris
  /// @return TriS<Type> tris() : triangle on the source plane
  TriS<Type> tris() const;
  /// @brief Return local magnification (the ratio of area of area on the lens
  /// plane and the source plane)
  /// @return Type mu() : the local magnification
  Type mu() const;
  /// @brief Return zone
  /// @param[in] &src : the source
  /// @param[in] mrg : mergin
  /// @param[in] &ml : the lens system
  /// @return Zone zone() : zone
  Zone zone(Source<Type> &src, Type mrg, Mlens<Type> &ml) const;
  /// @brief Return zone
  /// @param[in] &src : the source
  /// @param[in] &ml : the lens system
  /// @return Zone zone() : zone
  Zone zone(Source<Type> &src, Mlens<Type> &ml) const;
  /// @brief Return sensible number of the area of images
  /// @param[in] src : the source
  /// @return Type : sensible number of the area
  Type area_sn(Source<Type> src) const;
  /// @brief Return splitted triangles
  /// @param[in] &ml : the lens system
  /// @return std::array<TriMap<Type>, 2> : splitted trimap's
  std::array<TriMap<Type>, 2> split(Mlens<Type> &ml) const;
  /// @brief Return splitted triangles by n-times splitting
  /// @param[in] n : number of splittings
  /// @param[in] &ml : the lens system
  /// @return std::vector<TriMap<Type>> : splitted trimap's
  std::vector<TriMap<Type>> split(int n, Mlens<Type> &ml) const;
  /// @brief Prints parameters
  void print() const;
  /// @brief Return converted trimap to TriMap<float> (32-bit floating point)
  /// @return TriMap<float> : converted trimap
  TriMap<float> to_f() const;
  /// @brief Return converted trimap to TriMap<double> (64-bit floating point)
  /// @return TriMap<double> : converted trimap
  TriMap<double> to_d() const;
  /// @brief Return converted trimap to TriMap<long double> (128-bit floating
  /// point)
  /// @return TriMap<long double> : converted trimap
  TriMap<long double> to_l() const;
};

/// @brief Global const null TriMap
template <typename Type>
TriMap<Type> const TriMap_0 = TriMap<Type>(TriL_0<Type>, Mlens_0<Type>);

/**
 * template class TMQue
 * @brief Queue of TriMap
 * Type : float,  double,  long double
 *
 */
template <class Type> class TMQue : public std::deque<TriMap<Type>> {
  using size_type = typename std::deque<TriMap<Type>>::size_type;
  using std::deque<TriMap<Type>>::deque;

private:
  /// @brief Queue for grouping
  std::queue<TriMap<Type>> c_tmqp;
  /// @brief Queues of grou@s
  std::vector<TMQue<Type>> c_grps;

public:
  /// @brief Return specific TriMap
  /// @param[in] i : index of the TriMap
  /// @return TriMap<Type> &operator[] : The TriMap
  TriMap<Type> &operator[](size_type i);
  /// @brief Return specific TriMap
  /// @param[in] i : index of the TriMap
  /// @return TriMap<Type> &operator[] : The TriMap
  const TriMap<Type> &operator[](size_type i) const;
  /// @brief Constructor to initialize by a TriMap
  /// @param[in] tm : TriMap to initialize
  TMQue(TriMap<Type> tm) : std::deque<TriMap<Type>>({tm}) {};
  /// @brief Constructor to make from a square and a ml
  /// @param[in] sq : the square on the lens plane
  /// @param[in] ml : multiple lens
  TMQue(Square<Type> sq, Mlens<Type> ml) {
    auto tl0 = TriL<Type>(V2d<Type>(sq.max().x(), sq.min().y()), sq.max());
    auto tl1 = TriL<Type>(V2d<Type>(sq.min().x(), sq.max().y()), sq.min());
    auto tm0 = TriMap<Type>(tl0, ml);
    auto tm1 = TriMap<Type>(tl1, ml);
    *this = {tm0, tm1};
  }
  /// @brief Merge TMQue's
  /// @param[in] tm : TMQue to merge
  void operator+=(TMQue<Type> tm);
  /// @brief Return minimum x and y
  /// @return V2d<Type> min() : minimum x and y
  V2d<Type> min() const;
  /// @brief Return maximum x and y
  /// @return V2d<Type> min() : maximun x and y
  V2d<Type> max() const;
  /// @brief Return a quare including all
  /// @return Square<Type> include() : a square including all
  Square<Type> include() const;
  /// @brief Rtuen total area of TriL's
  /// @return Type area() : total area of TriL's
  Type area() const;
  /// @brief Return sensible number of the total area
  /// @param[in] src : the source
  /// @brief Type area_sn(Source<Type> src) : sensible number of the total area
  Type area_sn(Source<Type> src) const;
  /// @brief Return border TriMap's
  /// @param[in] src : the source
  /// @return TMQue<Type> border(Source<Type> src) : border TriMap's
  TMQue<Type> border(Source<Type> src) const;
  /// @brief Make critical curves and caustics
  /// @param[in] min_mu : lower limit of mu (local magnification)
  /// @return TMQue<Type> ccca(Type min_mu) : critical curves and caustics
  TMQue<Type> ccca(Type min_mu) const;
  /// @brief Initialize c_tmqp to start grouping
  void tmqp_init();
  /// @brief Pop c_tmqp
  void tmqp_pop();
  /// @brief Rotate c_tmqp queue
  void tmqp_rotate();
  /// @brief Return size of c_tmqp
  /// @return size_t tmqp_size() : size of c_tmqp
  size_t tmqp_size() const;
  /// @brief Return a group
  /// @return TMQue<Type> group() : a group
  TMQue<Type> group();
  /// @brief Return all groups
  /// @return std::vector<TMQue<Type>> groups() : all groups
  std::vector<TMQue<Type>> groups();
  /// @brief Return squares corresponding to groups
  /// @return Square's corresponding to groups
  std::vector<Square<Type>> squares();
  /// @brief Return tmqp
  /// @return std::queue<TriMap<Type>> tmqp() :  c_tmqp
  std::queue<TriMap<Type>> tmqp();
  /// @brief Return converted TMQue
  /// @return TMQue<Type> tmqp_to_tmque() : converted TMQue
  TMQue<Type> tmqp_to_tmque();
  /// @brief Return the front of c_tmqp
  /// @return TriMap<Type> tmqp_front() :  front of c_tmqp
  TriMap<Type> tmqp_front() const;
  /// @brief Return the back of c_tmqp
  /// @return TriMap<Type> tmqp_back() : back of c_tmqp
  TriMap<Type> tmqp_back() const;
  /// @brief Split this TMQue into two
  /// @param[in] &ml : multiple lens
  void split(Mlens<Type> &ml);
  /// @brief Split this TMQue n times
  /// @param[in] n : number of splitting
  /// @param[in] &ml : multiple lens
  void split(int n, Mlens<Type> &ml);
  /// @brief Print parameters
  void print() const;
};

/// @brief Global const null TMQue
template <typename Type>
const TMQue<Type> TMQue_0 = TMQue<Type>(TriMap_0<Type>);

/// @brief struct to express fractal parameters
template <typename Type> struct FrParam {
  /// @brief Mergin to identify Zone in the lens plane
  Type mrg_l = 0.1;
  /// @brief Mergin to identify Zone in the source plane
  Type mrg_s = 0.5;
  /// @brief Parameter to identify too_large in the lens plane
  Type max_drl = 0.2;
  /// @brief Lower limit of beta for single lens approximation
  Type sla_minbt = 0.8; // Lower limit of beta for SLA
  /// @brief Factor to identify the triangle is large compared to the source
  Type sla_fact_out = 3.0;
  /// @brief Factor to identify the source is in the triangle or not
  Type in_fact = 1.2;
  /// @brief A parameter. Currently dummy
  // Type mxdr = M_PI / 2.0;   // Currently dummy
  Type dummy2 = M_PI / 2.0; // Currently dummy
  /// @brief The factor to enlarge search area
  Type fact = 2.0;
  /// @brief Accumulate out queue
  bool acc_out = false;
  // bool acc_out = true;
  /// @brief Initial round number
  int rd_init = 10;
  /// @brief true if it's a benchmark
  bool bench = false;
  /// @brief benchmark value : start time in clock
  long begin = 0;
  /// @brief benchmark value : end time in clock
  long end = 0;
  /// @brief Returns FrParam<float>
  /// @return FrParam<float> to_f() : float parameters
  FrParam<float> to_f() const {
    FrParam<float> fparam;
    fparam.mrg_l = static_cast<float>(this->mrg_l);
    fparam.mrg_s = static_cast<float>(this->mrg_s);
    fparam.max_drl = static_cast<float>(this->max_drl);
    fparam.sla_minbt = static_cast<float>(this->sla_minbt);
    fparam.sla_fact_out = static_cast<float>(this->sla_fact_out);
    fparam.in_fact = static_cast<float>(this->in_fact);
    fparam.dummy2 = static_cast<float>(this->dummy2);
    fparam.fact = static_cast<float>(this->fact);
    fparam.acc_out = this->acc_out;
    fparam.acc_out = static_cast<float>(this->acc_out);
    fparam.rd_init = this->rd_init;
    fparam.bench = this->bench;
    fparam.begin = this->begin;
    fparam.end = this->end;
    return fparam;
  }
  /// @brief Returns FrParam<double>
  /// @return FrParam<double> to_d() : double parameters
  FrParam<double> to_d() const {
    FrParam<double> dparam;
    dparam.mrg_l = static_cast<double>(this->mrg_l);
    dparam.mrg_s = static_cast<double>(this->mrg_s);
    dparam.max_drl = static_cast<double>(this->max_drl);
    dparam.sla_minbt = static_cast<double>(this->sla_minbt);
    dparam.sla_fact_out = static_cast<double>(this->sla_fact_out);
    dparam.in_fact = static_cast<double>(this->in_fact);
    dparam.dummy2 = static_cast<double>(this->dummy2);
    dparam.fact = static_cast<double>(this->fact);
    dparam.acc_out = this->acc_out;
    dparam.acc_out = static_cast<double>(this->acc_out);
    dparam.rd_init = this->rd_init;
    dparam.bench = this->bench;
    return dparam;
  }
  /// @brief Returns FrParam<long double>
  /// @return FrParam<long double> to_l() : long double parameters
  FrParam<long double> to_l() const {
    FrParam<long double> lparam;
    lparam.mrg_l = static_cast<long double>(this->mrg_l);
    lparam.mrg_s = static_cast<long double>(this->mrg_s);
    lparam.max_drl = static_cast<long double>(this->max_drl);
    lparam.sla_minbt = static_cast<long double>(this->sla_minbt);
    lparam.sla_fact_out = static_cast<long double>(this->sla_fact_out);
    lparam.in_fact = static_cast<long double>(this->in_fact);
    lparam.dummy2 = static_cast<long double>(this->dummy2);
    lparam.fact = static_cast<long double>(this->fact);
    lparam.acc_out = this->acc_out;
    lparam.acc_out = static_cast<long double>(this->acc_out);
    lparam.rd_init = this->rd_init;
    lparam.bench = this->bench;
    lparam.begin = this->begin;
    lparam.end = this->end;
    lparam.begin = this->begin;
    lparam.end = this->end;
    return lparam;
  }
};

/// @brief Process branch
enum class EnProcess {
  /// @brief Continue
  Continue,
  /// @brief Break
  Break,
};

// Blank template class Time
// template <class Type> class Time;

/**
 * template class Snap
 * @brief class to take snap shot
 * Type : float,  double,  long double
 *
 */
template <class Type> class Snap {
private:
  /// @brief Multiple lens
  Mlens<Type> c_ml;
  /// @brief Source
  Source<Type> c_src;
  /// @brief Parameters
  FrParam<Type> c_param;
  /// @brief Number of round
  int c_round;
  /// @brief Time
  Time<Type> c_t;

public:
  /// @brief Constructor with ml, src, param, round, and t
  Snap(Mlens<Type> ml, Source<Type> src, FrParam<Type> param, int round,
       Time<Type> t)
      : c_ml(ml), c_src(src), c_param(param), c_round(round), c_t(t) {};
  /// @brief Constructor with ml, src, param, and round
  Snap(Mlens<Type> ml, Source<Type> src, FrParam<Type> param, int round)
      : c_ml(ml), c_src(src), c_param(param), c_round(round), c_t(0.0) {};
  /// @brief Set time
  void Set(Time<Type> t) { c_t = t; };
  /// @brief Set round
  void Set(int round) { c_round = round; };
  /// @brief Return ml
  /// @return Mlens<Type> ml() : multiple lens parameters
  Mlens<Type> ml() const { return c_ml; };
  /// @brief Return src
  /// @returnSource<Type> src()  : Source parameters
  Source<Type> src() const { return c_src; };
  /// @brief Return param
  /// @return FrParam<Type> param() : Fractal parameters
  FrParam<Type> param() const { return c_param; };
  /// @brief Return round
  /// @return int round() : round number
  int round() const { return c_round; };
  /// @brief Return time
  /// @return Time<Type> t() : time of the snap shot
  Time<Type> t() const { return c_t; };
  /// @brief Save snap shot to the file
  /// @param[im] file : file name to save
  void Save(std::string file) {
    std::ofstream out(file);
    out << c_ml.size() << std::endl;
    for (size_t i = 0; i < c_ml.size(); i++)
      out << c_ml[i].q() << " " << c_ml[i].p().x() << " " << c_ml[i].p().y()
          << std::endl;
    out << c_src.shape() << " " << c_src.gradation() << std::endl;
    out << c_src.brightness() << " " << c_src.limbdark() << std::endl;
    if (c_src.shape() == Shape::Circular)
      out << c_src.circ().cent().x() << " " << c_src.cent().y() << " "
          << c_src.circ().r() << std::endl;
    else
      out << c_src.ell().cent().x() << " " << c_src.ell().cent().y() << " "
          << c_src.ell().ab().x() << " " << c_src.ell().ab().y() << " "
          << c_src.ell().ph() << std::endl;
    out << c_round << std::endl;
    out << c_t.day() << std::endl;
    out.close();
  }
  /// @brief Restore snap shot from file
  /// @param[in] file : file name to read
  void Read(std::string file) {
    std::ifstream in(file);
    int n;
    in >> n;
    Type q, x, y;
    Mlens<Type> ml;
    for (int i = 0; i < n; i++) {
      in >> q >> x >> y;
      ml.push_back(Lens<Type>(q, V2d<Type>(x, y)));
    };
    c_ml = ml;
    Shape shp;
    Gradation gr;
    int ish, igr;
    in >> ish >> igr;
    if (ish == 0)
      shp = Shape::Circular;
    else
      shp = Shape::Elliptical;
    if (igr == 0)
      gr = Gradation::Uniform;
    else
      gr = Gradation::LimbDark;
    Type br, ld;
    in >> br >> ld;
    if (shp == Shape::Circular) {
      Type x, y, r;
      in >> x >> y >> r;
      if (gr == Gradation::Uniform)
        c_src = Source<Type>(Circle<Type>(V2d<Type>(x, y), r), br);
      else
        c_src = Source<Type>(Circle<Type>(V2d<Type>(x, y), r), br, ld);
    } else {
      Type x, y, a, b, ph;
      in >> x >> y >> a >> b >> ph;
      if (gr == Gradation::Uniform)
        c_src = Source<Type>(Ellipse<Type>(V2d<Type>(x, y), a, b, ph), br);
      else
        c_src = Source<Type>(Ellipse<Type>(V2d<Type>(x, y), a, b, ph), br, ld);
    }
    in >> c_round;
    Type td;
    in >> td;
    c_t = Time<Type>(td);
    in.close();
  }
};

template <typename Type>
std::array<V2d<Type>, 2> single_image(Mlens<Type> ml, Source<Type> src);

/// @brief Return Square that include all lenses and the source
/// @param[in] ml : multiple lens system
/// @param[in] src : the source star
/// @return Square<Type> include_all(Mlens<Type> ml, Source<Type> src) : Square
/// that include all lenses and the source
template <typename Type>
Square<Type> include_all(Mlens<Type> ml, Source<Type> src);

/**
 * template class FrTri
 * @brief The class for triangle fractal analysis
 * Type : float,  double,  long double
 *
 */
template <class Type> class FrTri {
private:
  /// @brief Two triangles to define area to search
  std::array<Triangle<Type>, 2> c_search;
  /// @brief "in" queue
  TMQue<Type> c_in;
  /// @brief "out" queue
  TMQue<Type> c_out;
  /// @brief "uncertain" queue
  TMQue<Type> c_uncertain;
  /// @brief Lens system
  Mlens<Type> c_ml;
  /// @brief Source
  Source<Type> c_src;
  /// @brief Parameters
  FrParam<Type> c_param;
  /// @brief Number of round
  int c_round = 0;

public:
  /// @brief Blank constructor
  FrTri() {
    auto tl0 = TriL(V2d<Type>(1.5, -1.5), V2d<Type>(1.5, 1.5));
    auto tl1 = TriL(V2d<Type>(-1.5, 1.5), V2d<Type>(-1.5, -1.5));
    auto tm0 = TriMap(tl0, c_ml);
    auto tm1 = TriMap(tl1, c_ml);
    c_uncertain = {tm0, tm1};
    while (c_round < c_param.rd_init) {
      c_uncertain.split(c_ml);
      c_round++;
    }
  }
  /// @brief Constructor, ml and src
  /// @param[in] ml : Lens system
  /// @param[in] src : Source
  FrTri(Mlens<Type> ml, Source<Type> src) : c_ml(ml), c_src(src) {
    c_ml.sort_i();
    auto sq = include_all(ml, src);
    c_search = sq.enlarge(c_param.fact).to_tri();
    TriL<Type> tl0(c_search[0].apex()[1], c_search[0].apex()[2]);
    TriL<Type> tl1(c_search[1].apex()[1], c_search[1].apex()[2]);
    c_uncertain = {TriMap<Type>(tl0, ml), TriMap<Type>(tl1, ml)};
    while (c_round < c_param.rd_init) {
      c_uncertain.split(c_ml);
      c_round++;
    }
  };
  /// @brief Constructor, ml, src, and fact
  /// @param[in] ml : Lens system
  /// @param[in] src : Source
  /// @param[in] fact : Factor of the region to analyze
  FrTri(Mlens<Type> ml, Source<Type> src, Type fact) : c_ml(ml), c_src(src) {
    c_ml.sort_i();
    c_param.fact = fact;
    auto sq = include_all(ml, src);
    c_search = sq.enlarge(c_param.fact).to_tri();
    TriL<Type> tl0(c_search[0].apex()[1], c_search[0].apex()[2]);
    TriL<Type> tl1(c_search[1].apex()[1], c_search[1].apex()[2]);
    c_uncertain = {TriMap<Type>(tl0, ml), TriMap<Type>(tl1, ml)};
    while (c_round < c_param.rd_init) {
      c_uncertain.split(c_ml);
      c_round++;
    }
  };
  /// @brief Constructor, ml, src, and param
  /// @param[in] ml : Lens system
  /// @param[in] src : Source
  /// @param[in] param : Parameters of the analysis
  FrTri(Mlens<Type> ml, Source<Type> src, FrParam<Type> param)
      : c_ml(ml), c_src(src), c_param(param) {
    c_ml.sort_i();
    auto sq = include_all(ml, src);
    c_search = sq.enlarge(c_param.fact).to_tri();
    TriL<Type> tl0(c_search[0].apex()[1], c_search[0].apex()[2]);
    TriL<Type> tl1(c_search[1].apex()[1], c_search[1].apex()[2]);
    c_uncertain = {TriMap<Type>(tl0, ml), TriMap<Type>(tl1, ml)};
    while (c_round < c_param.rd_init) {
      c_uncertain.split(c_ml);
      c_round++;
    }
  };
  /// @brief Constructor with snap
  /// @param[in] snap : snap shot
  FrTri(Snap<Type> snap)
      : c_ml(snap.ml()), c_src(snap.src()), c_param(snap.param()),
        c_round(snap.round()) {};
  /// @brief Constructor, unc
  /// @param[in] unc : "uncertain" queue
  FrTri(TMQue<Type> unc) : c_uncertain(unc) {};
  /// @brief Constructor, unc and param
  /// @param[in] unc : "uncertain" queue
  /// @param[in] param : Parameters for the analysis
  FrTri(TMQue<Type> unc, FrParam<Type> param)
      : c_uncertain(unc), c_param(param) {};
  /// @brief Destructor
  ~FrTri() {
    c_in.clear();
    c_out.clear();
    c_uncertain.clear();
    // std::cout << "cleared" << std::endl;
  };
  /// @brief Init queues
  /// @param[in] unc : "uncertain" queue
  void init(TMQue<Type> unc);
  /// @brief Move source and lenses
  /// @param[in] src : Source
  /// @param[in] ml : Lens system
  void move_to(Source<Type> src, Mlens<Type> ml);
  /// @brief Move the source
  /// @param[in] src : Source
  void move_to(Source<Type> src);
  /// brief Move the source to the position
  /// @param[in] ps : Position of the source
  void move_to(V2d<Type> ps);
  /// brief Move the source and the lenses to the positions
  /// @param[in] p_src : Position of the source
  /// @param[in] p_ml : Positions of the lenses
  void move_to(V2d<Type> p_src, std::vector<V2d<Type>> p_ml);
  /// @brief Split "uncertain" queue n-times without processing
  /// @param[in] n : Number of split processes
  /// @param[in] ml : Lens system
  void split_u(int n, Mlens<Type> &ml);
  /// @brief Return search area to find images
  /// @return Sqare<Type> search() : search area to find images
  Square<Type> search() const;
  /// @brief Reurn number of round
  /// @return Number of round
  int round() const;
  /// @brief Returns short side length of tril
  /// @return Type d() : length of short side of tril
  Type d() const;
  /// @brief Reurn "in" queue
  /// @return TMQue<Type> in() :  "in" queue
  TMQue<Type> in() const;
  /// @brief Reurn "out" queue
  /// @return TMQue<Type> out() :  "out" queue
  TMQue<Type> out() const;
  /// @brief Reurn "uncertain" queue
  /// @return TMQue<Type> uncertain() :  "uncertain" queue
  TMQue<Type> uncertain() const;
  /// @brief Reurn analysis patameters
  /// @return FrParam<Type> param() :  analysis parameters
  FrParam<Type> param() const;
  /// @brief Reurn the lens system
  /// @return Mlens<Type> ml() : the lens system
  Mlens<Type> ml() const;
  /// @brief Reurn the source
  /// @return Source<Type> src() : the source
  Source<Type> src() const;
  /// @brief Reurn the total area of "in"
  /// @return Type area_in() : Total area of "in" queue
  Snap<Type> snap() const { return Snap<Type>(c_ml, c_src, c_param, c_round); };
  /// @brief Return square include all
  /// @return square include all of the elements
  Square<Type> include() const;
  /// @brief Reurn the total area of "in"
  /// @return Type area_in() : Total area of "in" queue
  Type area_in() const;
  /// @brief Reurn the total area of "uncertain"
  /// @return Type area_uncertain() : Total area of "uncertain" queue
  Type area_uncertain() const;
  /// @brief Return a probable value of area obtained from "in" and "uncertain"
  /// @return Type area_iu() : area_in *+ area_uncertain / 2
  Type area_iu() const;
  /// @brief Reurn the total area of sensible number
  /// @return Type area_sn() : Total area of sensible number
  Type area_sn() const;
  /// @brief Reurn the total brightness
  /// @return Brightness<Type> brightness() : Total brightness
  Brightness<Type> brightness() const;
  /// @brief Reurn the total brightness of "in"
  /// @return Brightness<Type> br_in() : Total brightness of "in"
  Brightness<Type> br_in() const;
  /// @brief Reurn the total brightness of "out"
  /// @return Brightness<Type> br_out() : Total brightness of "out"
  Brightness<Type> br_out() const;
  /// @brief Reurn the total brightness of "uncertain"
  /// @return Brightness<Type> br_uncertain() : Total brightness of "uncertain"
  Brightness<Type> br_uncertain() const;
  /// @brief Return summation of the center of the images weighted by the area
  /// @return V2d<Type> sump_area() : summation of the center of the images
  /// weighted by the area
  V2d<Type> sump_area() const;
  /// @brief Return summation of the center of the images weighted by the
  /// brightness
  /// @return V2d<Type> sump_area() : summation of the center of the images
  /// weighted by the brightness
  V2d<Type> sump_br() const;
  /// @brief Return the image centroid weighted by the area
  /// @return V2d<Type> ic_area() : the image centroid weighted by the area
  V2d<Type> ic_area() const; // Image centroid in area
  /// @brief Return the image centroid weighted by the brightness
  /// @return V2d<Type> ic_area() : the image centroid weighted by the
  /// brightness
  V2d<Type> ic_br() const; // Image centroid in brightness
  /// @brief Return the astrometric shift of the images centroid weighted by the
  /// area
  /// @return V2d<Type> as_area() : the astrometric shift of the images centroid
  /// weighted by the area
  V2d<Type> as_area() const; // Astrometric shift in area
  /// @brief Return the astrometric shift of the images centroid weighted by the
  /// brightness
  /// @return V2d<Type> as_area() : the astrometric shift of the images centroid
  /// weighted by the brightness
  V2d<Type> as_br() const; // Astrometric shift in brightness
  /// @brief Process lens plane algorithm
  /// @return EnProcess process_l() : Continue or Break
  EnProcess process_l();
  // EnProcess process_il(size_t i);
  /// @brief Process source plane algorithm
  void process_s();
  /// @brief Process front TriMap on the uncertain queue
  void process_front();
  /// @brief Process single round of uncertain queue both lens and source planes
  void process();
  /// @brief Process n rounds
  /// @param[in] n : Number of rounds
  void process(int n);
  /// @brief Repeat process until the size of in queue exceed in_size
  /// @param[in] in_size : minimum size of in queue
  void process_while(int in_size);
  /// @brief Repeat process until the ratio of the areas of uncertain and in
  /// queues become less than ratio_unc_in)
  /// @param[in] ratio_unc_in : upper limit of the ratio
  void process_while(Type ratio_unc_in);
  /// @brief Process specifying relative precision, ml, src
  /// @param[in] prec : relative precision
  /// @param[in] &ml : lens system
  /// @param[in] &src : source
  void process(Type prec, Mlens<Type> &ml, Source<Type> &src);
  /// @brief Returns benchmark result
  /// @return Type ctime() : CPU time
  Type ctime() const {
    return static_cast<Type>(c_param.end - c_param.begin) / CLOCKS_PER_SEC;
  }
  /// @brief Clear queues
  void clear();
  /// @brief Returns converted FrTri<float>
  /// @return to_f() : converted FrTri<float>
  FrTri<float> to_f() const {
    return FrTri<float>(this->ml().to_f(), this->src().to_f(),
                        this->param().to_f());
  }
  /// @brief Returns converted FrTri<double>
  /// @return to_d() : converted FrTri<double>
  FrTri<double> to_d() const {
    return FrTri<double>(this->ml().to_d(), this->src().to_d(),
                         this->param().to_d());
  }
  /// @brief Returns converted FrTri<long double>
  /// @return to_l() : converted FrTri<long double>
  FrTri<long double> to_l() const {
    return FrTri<long double>(this->ml().to_l(), this->src().to_l(),
                              this->param().to_l());
  }
  /// @brief Print all parameters
  void print() const;
};

// template <typename Type> using Corner = std::array<V2d<Type>, 4>;
// template <typename Type>
// Corner<Type> Corner_0 = {V2d_0<Type>, V2d_0<Type>, V2d_0<Type>, V2d_0<Type>};
//
// template <typename Type> using Inner = std::array<V2d<Type>, 5>;
// template <typename Type>
// Inner<Type> Inner_0 = {V2d_0<Type>, V2d_0<Type>, V2d_0<Type>, V2d_0<Type>,
//                        V2d_0<Type>};
//
// template <class Type> class Quad {
// private:
//   Corner<Type> c_corner = Corner_0<Type>;
//
// public:
//   Quad()
//       : c_corner({V2d<Type>(-1.5, -1.5), V2d<Type>(1.5, -1.5),
//                   V2d<Type>(1.5, 1.5), V2d<Type>(-1.5, 1.5)}){};
//   Quad(V2d<Type> min, Type d)
//       : c_corner({min, min + V2d<Type>(d, 0.0), min + V2d<Type>(d, d),
//                   min + V2d<Type>(0.0, d)}){};
//   Quad(Corner<Type> cn) : c_corner(cn){};
//   Type d() const;
//   Type area() const;
//   bool adjoin(Quad<Type> qd) const;
//   Zone zone(Circle<Type> cl) const;
//   Zone zone(Circle<Type> cl, Type mrg) const;
//   Zone zone(Ellipse<Type> el) const;
//   Zone zone(Ellipse<Type> el, Type mrg) const;
//   bool border(Circle<Type> cl) const;
//   bool too_large(Circle<Type> cl) const;
//   bool too_large(Circle<Type> cl, Type mrgn) const;
//   bool in(V2d<Type> p) const;
//   bool in(Circle<Type> cl, Type mrg) const;
//   bool in(Circle<Type> cl) const;
//   Corner<Type> corner() const;
//   std::vector<Type> angles() const;
//   bool all_p() const;
//   bool all_n() const;
//   bool mixed() const;
//   V2d<Type> min() const;
//   V2d<Type> max() const;
//   Inner<Type> newpl() const;
//   Inner<Type> newps(Inner<Type> npl, Mlens<Type> ml) const;
//   std::vector<Quad<Type>> split() const;
//   std::vector<Quad<Type>> split(Inner<Type> newps) const;
//   std::array<TriL<Type>, 2> tril() const;
//   Quad<float> to_f() const;
//   Quad<double> to_d() const;
//   Quad<long double> to_l() const;
//   void print() const;
// };
//
// template <typename Type> Quad<Type> Quad_0 = Quad<Type>(Corner_0<Type>);
//
// template <class Type> class QuadMap {
// private:
//   Quad<Type> c_ql = Quad_0<Type>;
//   Quad<Type> c_qs = Quad_0<Type>;
//   Mlens<Type> c_ml = Mlens_0<Type>;
//   Source<Type> c_src;
//
// public:
//   QuadMap(Quad<Type> ql, Mlens<Type> ml, Source<Type> src)
//       : c_ql(ql), c_ml(ml), c_src(src) {
//     Corner<Type> av = Corner_0<Type>;
//     for (size_t i = 0; i < c_ql.corner().size(); i++)
//       av[i] = c_ml.beta(c_ql.corner()[i]);
//     c_qs = Quad<Type>(av);
//   };
//   QuadMap(Quad<Type> ql, Quad<Type> qs, Mlens<Type> ml, Source<Type> src)
//       : c_ql(ql), c_qs(qs), c_ml(ml), c_src(src){};
//   Mlens<Type> ml() const;
//   Source<Type> src() const;
//   Quad<Type> ql() const;
//   Quad<Type> qs() const;
//   Zone zone_s() const;
//   Zone zone_s(Type mrg) const;
//   Zone zone_l(int i) const;
//   std::vector<Zone> zone_l(Type mrg) const;
//   Zone zone_l(int i, Type mrg) const;
//   bool too_large_l(int i, Type drl) const;
//   bool too_large_s(Type drs) const;
//   std::vector<QuadMap<Type>> split() const;
//   QuadMap<float> to_f() const;
//   QuadMap<double> to_d() const;
//   QuadMap<long double> to_l() const;
//   void print() const;
// };
//
// template <typename Type>
// QuadMap<Type> QuadMap_0 =
//     QuadMap<Type>(Quad_0<Type>, Mlens_NaN<Type>, Source_0<Type>);
//
// template <class Type> class QMQue : public std::deque<QuadMap<Type>> {
//   using size_type = typename std::deque<QuadMap<Type>>::size_type;
//   using std::deque<QuadMap<Type>>::deque;
//   std::queue<QuadMap<Type>> c_qmqp;
//
// public:
//   QMQue(QuadMap<Type> qm) : std::deque<QuadMap<Type>>({qm}){};
//   QMQue(const QMQue<Type> &qmu) : std::deque<QuadMap<Type>>(qmu){};
//   QuadMap<Type> &operator[](size_type i);
//   const QuadMap<Type> &operator[](size_type i) const;
//   Type area() const;
//   QMQue<Type> all_p() const;
//   QMQue<Type> all_n() const;
//   QMQue<Type> mixed() const;
//   QMQue<Type> border() const;
//   void qmqp_init();
//   void qmqp_pop();
//   void qmqp_rotate();
//   size_t qmqp_size() const;
//   QMQue<Type> group();
//   std::vector<QMQue<Type>> groups();
//   std::queue<QuadMap<Type>> qmqp();
//   QMQue<Type> qmqp_to_qmque();
//   QuadMap<Type> qmqp_front() const;
//   QuadMap<Type> qmqp_back() const;
//   void print() const;
// };
//
// // template <typename Type> QMQue<Type> QMQue_0 =
// // QMQue<Type>(QuadMap_0<Type>);
//
// template <class Type> class FrQuad {
// private:
//   QMQue<Type> c_in;
//   QMQue<Type> c_out;
//   QMQue<Type> c_uncertain;
//   FrParam<Type> c_param = {0.5, 0.5, 0.5, 1.0};
//
// public:
//   FrQuad() : c_param({0.5, 0.5, 0.5, 1.0}){};
//   FrQuad(QMQue<Type> unc) : c_uncertain(unc), c_param({0.5, 0.5,
//   0.5, 1.0}){}; FrQuad(QMQue<Type> unc, FrParam<Type> param)
//       : c_uncertain(unc), c_param(param){};
//   QMQue<Type> in() const;
//   QMQue<Type> out() const;
//   QMQue<Type> uncertain() const;
//   FrParam<Type> param() const;
//   Type area_in() const;
//   Type area_uncertain() const;
//   EnProcess process_l();
//   void process_s();
//   void process_front();
//   void process();
//   void process(int n);
//   void process(Type prec);
//   void print() const;
// };
//
// //
// // Suppress instantations for float
// //
// extern template class Quad<float>;
// extern template class QuadMap<float>;
// extern template class QMQue<float>;
// extern template class FrQuad<float>;

extern template class TriL<float>;
extern template struct Brightness<float>;
extern template class TriS<float>;
extern template class TriMap<float>;
extern template class TMQue<float>;
extern template struct FrParam<float>;
extern template class Snap<float>;
extern template class FrTri<float>;

// using fQuad = Quad<float>;
// using fQuadMap = QuadMap<float>;
// using fQMQue = QMQue<float>;
// using fFrQuad = FrQuad<float>;

/// @brief Alias of TriL<float>
using fTriL = TriL<float>;
/// @brief Alias of TriS<float>
using fTriS = TriS<float>;
/// @brief Alias of TriMap<float>
using fTriMap = TriMap<float>;
/// @brief Alias of TMQue<float>
using fTMQue = TMQue<float>;
/// brief Alias of Snap<float>
using fSnap = Snap<float>;
/// @brief Alias of FrTri<float>
using fFrTri = FrTri<float>;

//
// For double
//
// extern template class Quad<double>;
// extern template class QuadMap<double>;
// extern template class QMQue<double>;
// extern template class FrQuad<double>;

extern template class TriL<double>;
extern template struct Brightness<double>;
extern template class TriS<double>;
extern template class TriMap<double>;
extern template class TMQue<double>;
extern template struct FrParam<double>;
extern template class Snap<double>;
extern template class FrTri<double>;

// using dQuad = Quad<double>;
// using dQuadMap = QuadMap<double>;
// using dQMQue = QMQue<double>;
// using dFrQuad = FrQuad<double>;

/// @brief Alias of TriL<double>
using dTriL = TriL<double>;
/// @brief Alias of TriS<double>
using dTriS = TriS<double>;
/// @brief Alias of TriMap<double>
using dTriMap = TriMap<double>;
/// @brief Alias of TMQue<double>
using dTMQue = TMQue<double>;
/// brief Alias of Snap<double>
using dSnap = Snap<double>;
/// @brief Alias of FrTri<double>
using dFrTri = FrTri<double>;

//
// For long double
//
// extern template class Quad<long double>;
// extern template class QuadMap<long double>;
// extern template class QMQue<long double>;
// extern template class FrQuad<long double>;

extern template class TriL<long double>;
extern template struct Brightness<long double>;
extern template class TriS<long double>;
extern template class TriMap<long double>;
extern template class TMQue<long double>;
extern template struct FrParam<long double>;
extern template class Snap<long double>;
extern template class FrTri<long double>;

// using lQuad = Quad<long double>;
// using lQuadMap = QuadMap<long double>;
// using lQMQue = QMQue<long double>;
// using lFrQuad = FrQuad<long double>;

/// @brief Alias of TriL<long double>
using lTriL = TriL<long double>;
/// @brief Alias of TriS<long double>
using lTriS = TriS<long double>;
/// @brief Alias of TriMap<long double>
using lTriMap = TriMap<long double>;
/// @brief Alias of TMQue<long double>
using lTMQue = TMQue<long double>;
/// brief Alias of Snap<long double>
using lSnap = Snap<long double>;
/// @brief Alias of FrTri<long double>
using lFrTri = FrTri<long double>;

#endif
