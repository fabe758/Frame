/**
 * @file Lcurve.h
 *
 * @brief header file for functions, classes to handle light curve
 * @author Fumio Abe
 *
 */

// define _INC_Lcurve_h_ to load only once
#ifndef _INC_Lcurve_h_
#define _INC_Lcurve_h_

#include "Fractal.h"
#include "Motion.h"
// #include <ctime>

template <typename Type> struct LcParam {
  bool record = false;
  std::vector<std::vector<int>> in_size;
  std::vector<std::vector<int>> unc_size;
  std::vector<std::vector<int>> out_size;
  std::vector<std::vector<Type>> in_area;
  std::vector<std::vector<Type>> unc_area;
  std::vector<std::vector<Type>> area_sn;
  std::vector<std::vector<int>> round;
  bool bench = false;
  long start;
  long stop;
  LcParam<float> to_f() const {
    LcParam<float> prm;
    prm.record = this->record;
    prm.in_size = this->in_size;
    prm.unc_size = this->unc_size;
    prm.out_size = this->out_size;
    // std::vector<std::vector<float>> ia(
    //     this->in_area.size(), std::vector<float>(this->in_area[0].size()));
    // std::cout << "ok, 1 : size = " << ia.size() << ", size = " <<
    // ia[0].size()
    //           << std::endl;
    // for (size_t i1 = 0; i1 < this->in_area.size(); i1++)
    //   for (size_t i2 = 0; i2 < this->in_area[i1].size(); i2++)
    //     ia[i1][i2] = static_cast<float>(this->in_area[i1][i2]);
    // prm.in_area = ia;
    // std::vector<std::vector<float>> ua(
    //     this->unc_area.size(), std::vector<float>(this->unc_area[0].size()));
    // for (size_t i1 = 0; i1 < this->unc_area.size(); i1++)
    //   for (size_t i2 = 0; i2 < this->unc_area[i1].size(); i2++)
    //     ua[i1][i2] = static_cast<float>(this->unc_area[i1][i2]);
    // prm.unc_area = ua;
    // std::vector<std::vector<float>> as(
    //     this->area_sn.size(), std::vector<float>(this->area_sn[0].size()));
    // for (size_t i1 = 0; i1 < this->area_sn.size(); i1++)
    //   for (size_t i2 = 0; i2 < this->area_sn[i1].size(); i2++)
    //     as[i1][i2] = static_cast<float>(this->area_sn[i1][i2]);
    // prm.area_sn = as;
    return prm;
  }
  LcParam<double> to_d() const {
    LcParam<double> prm;
    prm.record = this->record;
    prm.in_size = this->in_size;
    prm.unc_size = this->unc_size;
    prm.out_size = this->out_size;
    // std::vector<std::vector<double>> ia(
    //     this->in_area.size(), std::vector<double>(this->in_area[0].size()));
    // for (size_t i1 = 0; i1 < this->in_area.size(); i1++)
    //   for (size_t i2 = 0; i2 < this->in_area[i1].size(); i2++)
    //     ia[i1][i2] = static_cast<double>(this->in_area[i1][i2]);
    // prm.in_area = ia;
    // std::vector<std::vector<double>> ua(
    //     this->unc_area.size(),
    //     std::vector<double>(this->unc_area[0].size()));
    // for (size_t i1 = 0; i1 < this->unc_area.size(); i1++)
    //   for (size_t i2 = 0; i2 < this->unc_area[i1].size(); i2++)
    //     ua[i1][i2] = static_cast<double>(this->unc_area[i1][i2]);
    // prm.unc_area = ua;
    // std::vector<std::vector<double>> as(
    //     this->area_sn.size(), std::vector<double>(this->area_sn[0].size()));
    // for (size_t i1 = 0; i1 < this->area_sn.size(); i1++)
    //   for (size_t i2 = 0; i2 < this->area_sn[i1].size(); i2++)
    //     as[i1][i2] = static_cast<double>(this->area_sn[i1][i2]);
    // prm.area_sn = as;
    return prm;
  }
  LcParam<long double> to_l() const {
    LcParam<long double> prm;
    prm.record = this->record;
    prm.in_size = this->in_size;
    prm.unc_size = this->unc_size;
    prm.out_size = this->out_size;
    // std::vector<std::vector<long double>> ia(
    //     this->in_area.size(),
    //     std::vector<long double>(this->in_area[0].size()));
    // for (size_t i1 = 0; i1 < this->in_area.size(); i1++)
    //   for (size_t i2 = 0; i2 < this->in_area[i1].size(); i2++)
    //     ia[i1][i2] = static_cast<long double>(this->in_area[i1][i2]);
    // prm.in_area = ia;
    // std::vector<std::vector<long double>> ua(
    //     this->unc_area.size(),
    //     std::vector<long double>(this->unc_area[0].size()));
    // for (size_t i1 = 0; i1 < this->unc_area.size(); i1++)
    //   for (size_t i2 = 0; i2 < this->unc_area[i1].size(); i2++)
    //     ua[i1][i2] = static_cast<long double>(this->unc_area[i1][i2]);
    // prm.unc_area = ua;
    // std::vector<std::vector<long double>> as(
    //     this->area_sn.size(),
    //     std::vector<long double>(this->area_sn[0].size()));
    // for (size_t i1 = 0; i1 < this->area_sn.size(); i1++)
    //   for (size_t i2 = 0; i2 < this->area_sn[i1].size(); i2++)
    //     as[i1][i2] = static_cast<long double>(this->area_sn[i1][i2]);
    // prm.area_sn = as;
    return prm;
  }
};

/**
 * template class Lcurve
 * @brief Light curve
 *
 * Type : float,  double,  long double
 * Aliases : fLcurve (Lcurve<float),
 *           dLcurve (Lcurve<double),
 *           lLcurve (Lcurve<loaaa double)
 *
 */
template <class Type> class Lcurve {
private:
  /// @brief Source parameters
  std::vector<SrcMotion<Type>> c_srcm;
  /// @brief Multiple lens parameters
  MlMotion<Type> c_mlm;
  /// @brief A series of time
  std::vector<Time<Type>> c_t;
  /// @brief A series of magnification
  std::vector<std::vector<Type>> c_mag;
  std::vector<std::vector<Type>> c_mag2;
  /// @brief A series of astrometric shift
  std::vector<V2d<Type>> c_as;
  std::vector<V2d<Type>> c_as2;
  /// @brief Fractal parameters
  std::vector<std::vector<FrTri<Type>>> c_frtri;
  /// @brief Analysis parameters
  FrParam<Type> c_param;
  LcParam<Type> c_lcp;

public:
  /// @brief Constructor
  /// @param[in] srcm : multiple sorce motions
  /// @param[in] mlm : multiple lens mptions
  /// @param[in] t : a series of time
  Lcurve(std::vector<SrcMotion<Type>> srcm, MlMotion<Type> mlm,
         std::vector<Time<Type>> t)
      : c_srcm(srcm), c_mlm(mlm), c_t(t) {};
  /// @brief Constructor
  /// @param[in] srcm : multiple sorce motions
  /// @param[in] mlm : multiple lens mptions
  /// @param[in] t : a series of time
  /// @param[in] param : analysis parameters
  Lcurve(std::vector<SrcMotion<Type>> srcm, MlMotion<Type> mlm,
         std::vector<Time<Type>> t, FrParam<Type> param)
      : c_srcm(srcm), c_mlm(mlm), c_t(t), c_param(param) {};
  Lcurve(std::vector<SrcMotion<Type>> srcm, MlMotion<Type> mlm,
         std::vector<Time<Type>> t, FrParam<Type> param, LcParam<Type> lcp)
      : c_srcm(srcm), c_mlm(mlm), c_t(t), c_param(param), c_lcp(lcp) {};
  /// @brief Constructor
  /// @param[in] srcm : single sorce motions
  /// @param[in] ml : multiple lens parameters
  /// @param[in] t : a series of time
  Lcurve(SrcMotion<Type> srcm, Mlens<Type> ml, std::vector<Time<Type>> t)
      : c_srcm({srcm}), c_mlm(MlMotion<Type>(ml)), c_t(t) {};
  /// @brief Constructor
  /// @param[in] srcm : single sorce motions
  /// @param[in] ml : multiple lens parameters
  /// @param[in] t : a series of time
  /// @param[in] param : analysis parameters
  Lcurve(SrcMotion<Type> srcm, Mlens<Type> ml, std::vector<Time<Type>> t,
         FrParam<Type> param)
      : c_srcm({srcm}), c_mlm(MlMotion<Type>(ml)), c_t(t), c_param(param) {};
  /// @brief Constructor
  /// @param[in] src : static single sorce
  /// @param[in] mlm : multiple lens mptions
  /// @param[in] t : a series of time
  Lcurve(Source<Type> src, MlMotion<Type> mlm, std::vector<Time<Type>> t)
      : c_srcm({SrcMotion<Type>(src)}), c_mlm(mlm), c_t(t) {};
  /// @brief Destructor : to make memory free
  ~Lcurve() {
    this->clear();
    // c_mlm.clear();
    // c_t.clear();
    // c_as.clear();
    // for (size_t is = 0; is < c_srcm.size(); is++) {
    //   c_mag[is].clear();
    //   c_frtri[is].clear();
    // }
    // c_mag.clear();
    // c_frtri.clear();
    // c_srcm.clear();
  };
  void set_record(bool record) { c_lcp.record = record; };
  void set_bench(bool bench) { c_lcp.bench = bench; };
  /// @brief Return the series of time
  /// @return std::vector<Time<Type>> time() : the series of time
  std::vector<Time<Type>> time() const { return c_t; };
  /// @brief Return the source motions
  /// @return std::vector<SrcMotion<Type>> srcm() : source motions
  std::vector<SrcMotion<Type>> srcm() const { return c_srcm; };
  /// @brief Return multiple lens motion
  /// @return MlMotion<Type> mlm() : multiple lens motion
  MlMotion<Type> mlm() const { return c_mlm; };
  std::vector<std::vector<FrTri<Type>>> frtri() const { return c_frtri; };
  Source<Type> src(int it) const;
  Source<Type> src(Time<Type> t) const;
  Mlens<Type> ml(int it) const;
  Mlens<Type> ml(Time<Type> t) const;
  FrTri<Type> frtri(int it) const;
  FrTri<Type> frtri(Time<Type> t) const;
  /// @brief Clculate for srcm[is] and t[it] for np times and return result
  /// @param[in] is : index of source
  /// @param[in] it : index of time
  /// @param[in] round : number of round of calculation
  /// @return FrTri<Type> calc(int is, int it, int np) : result
  FrTri<Type> calc(int is, int it, int round) {
    FrTri<Type> frtri(c_mlm.ml(c_t[it]), c_srcm[is].src(c_t[it]), c_param);
    frtri.process(round);
    return frtri;
  };
  /// @brief Clculate
  /// @param[in] round : number of round of calculation
  void process(int round);
  void process_while(int in_size);
  void process_while(Type ratio_unc_in);
  void process_while2(int in_size);
  void process_while2(Type ratio_unc_in);
  std::vector<std::vector<Type>> mags2() const;
  std::vector<Type> mag2() const;
  std::vector<V2d<Type>> as2() const;
  /// @brief Return magnifications
  /// @return std::vector<std::vector<Type>> mag() : series of magnifications
  std::vector<std::vector<Type>> mags() const;
  std::vector<Type> mag() const;
  /// @brief Return astrometric shifts
  /// @return std::vector<V2d<Type>> as() : series of astrometric shifts
  std::vector<V2d<Type>> as() const;
  /// @brief Return snap shots between min and max it
  /// @param[in] min : lower end of it
  /// @param[in] max : upper end of it
  /// @return std::vector<std::vector<Snap<Type>>> snap(size_t min, size_t max)
  /// : selected snap shots
  std::vector<std::vector<Snap<Type>>> snap(size_t min, size_t max) const;
  /// @brief Return selected snap shots selected by the function select
  /// @return std::vector<std::vector<Snap<Type>>>
  /// snap(std::function<bool(Time<Type>) : selected snap shots
  std::vector<std::vector<Snap<Type>>>
  snap(std::function<bool(Time<Type>)> select) const;
  /// @brief Return selected snap shots selected by the function select
  /// @return  : selected snap shots
  std::vector<std::vector<Snap<Type>>>
  snap(std::function<bool(FrTri<Type>)> select) const;
  /// @brief Clear c_frti
  void clear() {
    c_mlm.clear();
    c_t.clear();
    c_as.clear();
    c_as2.clear();
    for (size_t is = 0; is < c_srcm.size(); is++) {
      c_mag[is].clear();
      c_mag2[is].clear();
      c_frtri[is].clear();
    }
    c_mag.clear();
    c_frtri.clear();
    c_srcm.clear();
    c_mag.clear();
    c_mag2.clear();
    c_frtri.clear();
    c_srcm.clear();
  };
  FrParam<Type> param() const { return c_param; };
  LcParam<Type> lcp() const { return c_lcp; };
  Type ctime() const;
  Lcurve<float> to_f() const {
    std::vector<SrcMotion<float>> vs;
    for (auto sm : this->srcm())
      vs.push_back(sm.to_f());
    std::vector<Time<float>> vt;
    for (auto tm : this->time())
      vt.push_back(tm.to_f());
    return Lcurve<float>(vs, this->mlm().to_f(), vt, this->param().to_f());
  }
  Lcurve<double> to_d() const {
    std::vector<SrcMotion<double>> vs;
    for (auto sm : this->srcm())
      vs.push_back(sm.to_d());
    std::vector<Time<double>> vt;
    for (auto tm : this->time())
      vt.push_back(tm.to_d());
    return Lcurve<double>(vs, this->mlm().to_d(), vt, this->param().to_d(),
                          this->lcp().to_d());
  }
  Lcurve<long double> to_l() const {
    std::vector<SrcMotion<long double>> vs;
    for (auto sm : this->srcm())
      vs.push_back(sm.to_l());
    std::vector<Time<long double>> vt;
    for (auto tm : this->time())
      vt.push_back(tm.to_l());
    return Lcurve<long double>(vs, this->mlm().to_l(), vt, this->param().to_l(),
                               this->lcp().to_l());
  }
};

//
//  Suppress instantiations of templates of classes
//
//   For float
//
extern template struct LcParam<float>;
extern template class Lcurve<float>;

/// @brief Alias of Lcurve<float>
using fLcurve = Lcurve<float>;

//
//  For double
//
extern template struct LcParam<double>;
extern template class Lcurve<double>;

/// @brief Alias of Lcurve<double>
using dLcurve = Lcurve<double>;

//
//  For long double
//
extern template class Lcurve<long double>;

/// @brief Alias of Lcurve<long double>
extern template struct LcParam<long double>;
using lLcurve = Lcurve<long double>;

#endif
