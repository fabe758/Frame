/**
 * @file Fractal.cpp
 +
 * @brief Class library to handle fractal algorith
 * @author Fumio Abe
 * @date 25 April,  2023
 *
 * @details class library for self-similar division of squares
 * in the lensing plane and corresponding corner points in the
 * source plane (QuadMap). Queues to handle QuadMaps.
 * @note
 *
 */

/**
 * include headers  for this project
 * all c++ headers are loaded in Geom.h
 * Geom.h,  Lens.h,  and Source.h are loaded in Fractal.h
 *
 */
#include "Fractal.h"

/**
 * std::string zmap(Zone zone);
 * brief Function to return a string that represent the zone (used to print)
 * Return enum class Zone
 * @return std::string zmap() : a string that represent the zone
 * @param[in] zone : Zone parameters
 *
 */
std::string zmap(Zone zone) {
  if (zone == Zone::in)
    return "Zone::in";
  if (zone == Zone::uncertain)
    return "Zone::uncertain";
  return "Zone::out";
}

template <typename Type> V2d<Type> TriL<Type>::right() const {
  return this->c_apex[0];
}

template <typename Type> std::array<V2d<Type>, 3> TriL<Type>::apex() const {
  return this->c_apex;
}

template <typename Type> Type TriL<Type>::d() const {
  return this->sides()[0].length();
}

template <typename Type>
bool TriL<Type>::too_large(Circle<Type> cl, Type &drl) const {
  return this->d() > (cl.r() * drl);
}

template <class Type> Type TriL<Type>::area() const {
  return 0.5 * (this->c_apex[1] - this->c_apex[0]).n2();
}

template <class Type> Zone TriL<Type>::zone(Circle<Type> cl) const {
  return this->zone(cl, 0.2);
}

template <class Type> Zone TriL<Type>::zone(Circle<Type> cl, Type mrg) const {
  Type rmrg = (this->c_apex[0] - this->c_apex[1]).norm() * mrg;
  auto cns = this->c_apex;
  auto sds = this->sides();
  V2d<Type> maxap =
      (*std::max_element(cns.begin(), cns.end(), [&cl](auto a, auto b) {
        return (a - cl.cent()).norm() < (b - cl.cent()).norm();
      }));
  Type maxd = (maxap - cl.cent()).norm();
  if (maxd - cl.r() < -rmrg)
    return Zone::in;
  auto mnd = (std::min_element(sds.begin(), sds.end(), [&cl](auto a, auto b) {
    return a.distance(cl.cent()) < b.distance(cl.cent());
  }));
  Type mind = mnd->distance(cl.cent());
  if ((mind - cl.r()) > rmrg)
    return Zone::out;
  return Zone::uncertain;
}

template <typename Type>
Zone TriL<Type>::zone(Ellipse<Type> el, Type mrg) const {
  V2d<Type> eqvap0 = el.eqvCirc01(this->c_apex[0]);
  V2d<Type> eqvap1 = el.eqvCirc01(this->c_apex[1]);
  TriL<Type> eqvtril(eqvap0, eqvap1);
  Circle<Type> circ(V2d<Type>(0.0, 0.0), 1.0);
  return eqvtril.zone(circ, mrg);
}

template <class Type> bool TriL<Type>::border(Circle<Type> cl) const {
  auto cns = this->c_apex;
  V2d<Type> maxap =
      (*std::max_element(cns.begin(), cns.end(), [&cl](auto a, auto b) {
        return (a - cl.cent()).norm() < (b - cl.cent()).norm();
      }));
  Type maxd = (maxap - cl.cent()).norm();
  V2d<Type> minap =
      (*std::min_element(cns.begin(), cns.end(), [&cl](auto a, auto b) {
        return (a - cl.cent()).norm() < (b - cl.cent()).norm();
      }));
  Type mind = (minap - cl.cent()).norm();
  if (maxd > cl.r() && mind < cl.r())
    return true;
  return false;
}

template <class Type> bool TriL<Type>::adjoin(TriL<Type> tril) const {
  auto cns0 = this->c_apex;
  auto cns1 = tril.apex();
  for (size_t i = 0; i < cns0.size(); i++) {
    if (cns0[i].between(cns1[0], cns1[1]))
      return true;
    if (cns0[i].between(cns1[1], cns1[2]))
      return true;
    if (cns0[i].between(cns1[2], cns1[0]))
      return true;
  }
  for (size_t i = 0; i < cns1.size(); i++) {
    if (cns1[i].between(cns0[0], cns0[1]))
      return true;
    if (cns1[i].between(cns0[1], cns0[2]))
      return true;
    if (cns1[i].between(cns0[2], cns0[0]))
      return true;
  }
  return false;
}

template <class Type> std::array<TriL<Type>, 2> TriL<Type>::split() const {
  auto cns = this->c_apex;
  return {TriL((cns[1] + cns[2]) * 0.5, cns[2]),
          TriL((cns[1] + cns[2]) * 0.5, cns[0])};
}

template <class Type> TriL<float> TriL<Type>::to_f() const {
  return TriL<float>(this->c_apex[0].to_f(), this->c_apex[1].to_f());
}

template <class Type> TriL<double> TriL<Type>::to_d() const {
  return TriL<double>(this->c_apex[0].to_d(), this->c_apex[1].to_d());
}

template <class Type> TriL<long double> TriL<Type>::to_l() const {
  return TriL<long double>(this->c_apex[0].to_l(), this->c_apex[1].to_l());
}

template <class Type> std::array<V2d<Type>, 3> TriS<Type>::apex() const {
  return this->c_apex;
}

template <typename Type> Type TriS<Type>::maxd() const {
  auto sds = this->sides();
  return std::max_element(
             sds.begin(), sds.end(),
             [](auto a, auto b) { return a.length() < b.length(); })
      ->length();
}

template <class Type> Zone TriS<Type>::zone(Circle<Type> cl) const {
  return this->zone(cl, 0.2);
}

template <class Type> Zone TriS<Type>::zone(Circle<Type> cl, Type mrg) const {
  Type rmrg = (this->c_apex[0] - this->c_apex[1]).norm() * mrg;
  auto cns = this->c_apex;
  auto sds = this->sides();
  V2d<Type> maxap =
      (*std::max_element(cns.begin(), cns.end(), [&cl](auto a, auto b) {
        return (a - cl.cent()).norm() < (b - cl.cent()).norm();
      }));
  Type maxd = (maxap - cl.cent()).norm();
  if (maxd - cl.r() < -rmrg)
    return Zone::in;

  V2d<Type> minap =
      (*std::min_element(sds.begin(), sds.end(), [&cl](auto a, auto b) {
        return a.distance(cl.cent()) < b.distance(cl.cent());
      })).closest(cl.cent());
  Type mind = (minap - cl.cent()).norm();
  if ((mind - cl.r()) > rmrg)
    return Zone::out;

  return Zone::uncertain;
}

template <class Type> Zone TriS<Type>::zone(Ellipse<Type> el) const {
  return this->zone(el, 0.2);
}

template <class Type> Zone TriS<Type>::zone(Ellipse<Type> el, Type mrg) const {
  TriS<Type> tris({el.eqvCirc01(this->apex()[0]), el.eqvCirc01(this->apex()[1]),
                   el.eqvCirc01(this->apex()[2])});
  Circle<Type> circ(V2d<Type>(0.0, 0.0), 1.0);
  return tris.zone(circ, mrg);
}

template <class Type> bool TriS<Type>::border(Circle<Type> cl) const {
  auto cns = this->c_apex;
  V2d<Type> maxap =
      (*std::max_element(cns.begin(), cns.end(), [&cl](auto a, auto b) {
        return (a - cl.cent()).norm() < (b - cl.cent()).norm();
      }));
  Type maxd = (maxap - cl.cent()).norm();
  V2d<Type> minap =
      (*std::min_element(cns.begin(), cns.end(), [&cl](auto a, auto b) {
        return (a - cl.cent()).norm() < (b - cl.cent()).norm();
      }));
  Type mind = (minap - cl.cent()).norm();
  if (maxd > cl.r() && mind < cl.r())
    return true;
  return false;
}

template <class Type> bool TriS<Type>::border(Ellipse<Type> el) const {
  TriS<Type> tris({el.eqvCirc01(this->apex()[0]), el.eqvCirc01(this->apex()[1]),
                   el.eqvCirc01(this->apex()[2])});
  Circle<Type> circ(V2d<Type>(0.0, 0.0), 1.0);
  return tris.border(circ);
}

template <typename Type>
std::vector<Type> TriS<Type>::sb(Source<Type> src) const {
  return {src.sb(this->c_apex[0]), src.sb(this->c_apex[1]),
          src.sb(this->c_apex[2])};
}

template <typename Type>
Brightness<Type> TriS<Type>::brightness(Source<Type> src) const {
  Type mean = this->area() *
              std::accumulate(this->sb(src).begin(), this->sb(src).end(), 0.0) /
              this->sb(src).size();
  Type sqs = this->area() * this->area() *
             std::inner_product(this->sb(src).begin(), this->sb(src).end(),
                                this->sb(src).begin(), 0.0) /
             this->sb(src).size();
  return {mean, std::sqrt(sqs - mean * mean)};
}

template <class Type>
std::array<TriS<Type>, 2> TriS<Type>::split(V2d<Type> right) const {
  return {TriS<Type>({right, this->c_apex[2], this->c_apex[0]}),
          TriS<Type>({right, this->c_apex[0], this->c_apex[1]})};
}

template <class Type> TriS<float> TriS<Type>::to_f() const {
  auto aps = this->apex();
  return TriS<float>({aps[0].to_f(), aps[1].to_f(), aps[2].to_f()});
}

template <class Type> TriS<double> TriS<Type>::to_d() const {
  auto aps = this->apex();
  return TriS<double>({aps[0].to_d(), aps[1].to_d(), aps[2].to_d()});
}

template <class Type> TriS<long double> TriS<Type>::to_l() const {
  auto aps = this->apex();
  return TriS<long double>({aps[0].to_l(), aps[1].to_l(), aps[2].to_l()});
}

template <class Type> TriL<Type> TriMap<Type>::tril() const { return c_tril; }

template <class Type> TriS<Type> TriMap<Type>::tris() const { return c_tris; }

template <typename Type> Type TriMap<Type>::mu() const {
  return c_tril.area() / c_tris.area();
}

template <class Type>
Zone TriMap<Type>::zone(Source<Type> &src, Type mrg, Mlens<Type> &ml) const {
  if (src.shape() == Shape::Circular) {
    for (auto lens : ml)
      if (c_tril.in(lens.p()))
        if (c_tris.zone(src.circ(), mrg) == Zone::in)
          return Zone::uncertain;
    return c_tris.zone(src.circ(), mrg);
    return Zone::uncertain;
  } else if (src.shape() == Shape::Elliptical) {
    for (auto lens : ml)
      if (c_tril.in(lens.p()))
        if (c_tris.zone(src.ell(), mrg) == Zone::in)
          return Zone::uncertain;
    return c_tris.zone(src.ell(), mrg);
    return Zone::uncertain;
  }
  return Zone::uncertain;
}

template <class Type>
Zone TriMap<Type>::zone(Source<Type> &src, Mlens<Type> &ml) const {
  return this->zone(src, 0.2, ml);
}

template <typename Type> Type TriMap<Type>::area_sn(Source<Type> src) const {
  if (src.shape() == Shape::Circular)
    return c_tril.area() * src.circ().ratio_ia(c_tris);
  return c_tril.area() * src.ell().ratio_ia(c_tris);
}

template <class Type>
std::array<TriMap<Type>, 2> TriMap<Type>::split(Mlens<Type> &ml) const {
  auto spl = c_tril.split();
  auto sps = c_tris.split(ml.beta(spl[0].right()));
  std::array<TriMap<Type>, 2> tms{TriMap(spl[0], sps[0]),
                                  TriMap(spl[1], sps[1])};
  return tms;
}

template <class Type>
std::vector<TriMap<Type>> TriMap<Type>::split(int n, Mlens<Type> &ml) const {
  std::vector<TriMap<Type>> vtm;
  for (int i = 0; i < n; i++) {
    std::array<TriMap<Type>, 2> atm = this->split(ml);
    vtm.push_back(atm[0]);
    vtm.push_back(atm[1]);
  }
  return vtm;
}

template <typename Type> void TriMap<Type>::print() const {
  std::cout << "TriL : " << std::endl;
  c_tril.print();
  std::cout << "TriS : " << std::endl;
  c_tris.print();
}

template <class Type> TriMap<float> TriMap<Type>::to_f() const {
  return TriMap<float>(c_tril.to_f(), c_tris.to_f());
}

template <class Type> TriMap<double> TriMap<Type>::to_d() const {
  return TriMap<double>(c_tril.to_d(), c_tris.to_d());
}

template <class Type> TriMap<long double> TriMap<Type>::to_l() const {
  return TriMap<long double>(c_tril.to_l(), c_tris.to_l());
}

template <class Type> TriMap<Type> &TMQue<Type>::operator[](size_type i) {
  return this->at(i);
}

template <class Type>
const TriMap<Type> &TMQue<Type>::operator[](size_type i) const {
  return this->at(i);
}

template <typename Type> void TMQue<Type>::operator+=(TMQue<Type> tm) {
  for (size_t i = 0; i < tm.size(); i++)
    this->push_back(tm[i]);
}
template <typename Type> V2d<Type> TMQue<Type>::min() const {
  return V2d<Type>(
      (*std::min_element(this->begin(), this->end(),
                         [](auto a, auto b) {
                           return a.tril().min().x() < b.tril().min().x();
                         }))
          .tril()
          .min()
          .x(),
      (*std::min_element(this->begin(), this->end(),
                         [](auto a, auto b) {
                           return a.tril().min().y() < b.tril().min().y();
                         }))
          .tril()
          .min()
          .y());
}

template <typename Type> V2d<Type> TMQue<Type>::max() const {
  return V2d<Type>(
      (*std::max_element(this->begin(), this->end(),
                         [](auto a, auto b) {
                           return a.tril().max().x() < b.tril().max().x();
                         }))
          .tril()
          .max()
          .x(),
      (*std::max_element(this->begin(), this->end(),
                         [](auto a, auto b) {
                           return a.tril().max().y() < b.tril().max().y();
                         }))
          .tril()
          .max()
          .y());
}

template <typename Type> Square<Type> TMQue<Type>::include() const {
  V2d<Type> cent = (this->min() + this->max()) / 2.0;
  Type d = std::max(this->max().x() - cent.x(), this->max().y() - cent.y());
  return Square<Type>(cent, d);
}

template <class Type> Type TMQue<Type>::area() const {
  Type ar = 0.0;
  for (size_type i = 0; i < this->size(); i++)
    ar += this->at(i).tril().area();
  return ar;
}

template <typename Type> Type TMQue<Type>::area_sn(Source<Type> src) const {
  Type ar = 0.0;
  for (size_t i = 0; i < this->size(); i++)
    ar += this->at(i).area_sn(src);
  return ar;
}

template <class Type> TMQue<Type> TMQue<Type>::border(Source<Type> src) const {
  TMQue<Type> bd;
  for (size_type i = 0; i < this->size(); i++)
    if (this->at(i).tris().border(src.circ()))
      bd.push_back(this->at(i));
  return bd;
}

template <typename Type> TMQue<Type> TMQue<Type>::ccca(Type min_mu) const {
  TMQue<Type> tmq;
  for (size_t i = 0; i < this->size(); i++)
    if (this->at(i).mu() > min_mu)
      tmq.push_back(this->at(i));
  return tmq;
}

template <class Type> void TMQue<Type>::tmqp_init() {
  while (!c_tmqp.empty())
    c_tmqp.pop();
  if (!c_grps.empty())
    c_grps.clear();
  for (TriMap<Type> tm : *this)
    c_tmqp.push(tm);
}

template <class Type> void TMQue<Type>::tmqp_pop() { c_tmqp.pop(); }

template <class Type> void TMQue<Type>::tmqp_rotate() {
  c_tmqp.push(c_tmqp.front());
  c_tmqp.pop();
}

template <class Type> size_t TMQue<Type>::tmqp_size() const {
  return c_tmqp.size();
}

template <class Type> TMQue<Type> TMQue<Type>::group() {
  TMQue<Type> tmq;
  tmq.push_back(c_tmqp.front());
  c_tmqp.pop();
  size_t it = 0;
  while (it < tmq.size()) {
    size_t nc = c_tmqp.size();
    for (size_t ic = 0; ic < nc; ic++) {
      if (tmq[it].tril().adjoin(c_tmqp.front().tril())) {
        // std::cout << "adjoint found : " << std::endl;
        // c_tmqp.front().tril().print();
        tmq.push_back(c_tmqp.front());
        c_tmqp.pop();
      } else {
        c_tmqp.push(c_tmqp.front());
        c_tmqp.pop();
      }
    }
    it++;
  }
  return tmq;
}

template <class Type> std::vector<TMQue<Type>> TMQue<Type>::groups() {
  this->tmqp_init();
  while (c_tmqp.size() > 0) {
    TMQue<Type> grp = this->group();
    c_grps.push_back(grp);
  }
  return c_grps;
}

template <class Type> std::vector<Square<Type>> TMQue<Type>::squares() {
  if (!c_grps.empty()) {
    c_grps = this->groups();
  }
  std::vector<Square<Type>> sqs;
  for (auto grp : c_grps)
    sqs.push_back(grp.include());
  return sqs;
}

template <class Type> std::queue<TriMap<Type>> TMQue<Type>::tmqp() {
  std::queue<TriMap<Type>> tmqp;
  size_t nq = c_tmqp.size();
  for (size_t i = 0; i < nq; i++) {
    tmqp.push(c_tmqp.front());
    c_tmqp.push(c_tmqp.front());
    c_tmqp.pop();
  }
  return tmqp;
}

template <class Type> TMQue<Type> TMQue<Type>::tmqp_to_tmque() {
  TMQue<Type> tmq;
  size_t nq = c_tmqp.size();
  for (size_t i = 0; i < nq; i++) {
    tmq.push_back(c_tmqp.front());
    c_tmqp.push(c_tmqp.front());
    c_tmqp.pop();
  }
  return tmq;
}

template <class Type> TriMap<Type> TMQue<Type>::tmqp_front() const {
  return c_tmqp.front();
}

template <class Type> TriMap<Type> TMQue<Type>::tmqp_back() const {
  return c_tmqp.back();
}

template <typename Type> void TMQue<Type>::split(Mlens<Type> &ml) {
  size_t n = this->size();
  for (size_t i = 0; i < n; i++) {
    auto atm = this->at(i).split(ml);
    this->push_back(atm[0]);
    this->push_back(atm[1]);
  }
  this->erase(this->begin(), this->begin() + n);
}

template <typename Type> void TMQue<Type>::split(int n, Mlens<Type> &ml) {
  for (int i = 0; i < n; i++)
    this->split(ml);
}

template <class Type> void TMQue<Type>::print() const {
  std::cout << "TMQue : size = " << this->size() << std::endl;
}

template <typename Type>
std::array<V2d<Type>, 2> single_image(Mlens<Type> ml, Source<Type> src) {
  V2d<Type> bt = (src.cent() - ml.cm()) / std::sqrt(ml.sumq());
  Type b = bt.norm();
  Type sqr = std::sqrt(b * b + 4.0);
  std::array<Type, 2> th({(b + sqr) / 2, (b - sqr) / 2});
  return std::array<V2d<Type>, 2>(
      {bt.unit() * th[0] + ml.cm(), bt.unit() * th[1] + ml.cm()});
}

template <typename Type>
Square<Type> include_all(Mlens<Type> ml, Source<Type> src) {
  std::vector<Square<Type>> vsq;
  if (src.shape() == Shape::Circular)
    vsq.push_back(src.circ().include());
  else
    vsq.push_back(src.ell().include());
  for (auto lens : ml)
    vsq.push_back(lens.circ().include());
  return include_all(vsq);
}

template <typename Type> void FrTri<Type>::init(TMQue<Type> unc) {
  this->clear();
  c_uncertain = unc;
}

template <typename Type>
void FrTri<Type>::move_to(Source<Type> src, Mlens<Type> ml) {
  c_src = src;
  c_ml = ml;
  auto sq = include_all(c_ml, c_src);
  c_search = sq.enlarge(c_param.fact).to_tri();
  TriL<Type> tl0(c_search[0].apex()[1], c_search[0].apex()[2]);
  TriL<Type> tl1(c_search[1].apex()[1], c_search[1].apex()[2]);
  c_uncertain = {TriMap<Type>(tl0, c_ml), TriMap<Type>(tl1, c_ml)};
}

template <typename Type> void FrTri<Type>::move_to(Source<Type> src) {
  c_src = src;
  auto sq = include_all(c_ml, c_src);
  c_search = sq.enlarge(c_param.fact).to_tri();
  TriL<Type> tl0(c_search[0].apex()[1], c_search[0].apex()[2]);
  TriL<Type> tl1(c_search[1].apex()[1], c_search[1].apex()[2]);
  c_uncertain = {TriMap<Type>(tl0, c_ml), TriMap<Type>(tl1, c_ml)};
}

template <typename Type> void FrTri<Type>::move_to(V2d<Type> ps) {
  c_src.move_to(ps);
  auto sq = include_all(c_ml, c_src);
  c_search = sq.enlarge(c_param.fact).to_tri();
  TriL<Type> tl0(c_search[0].apex()[1], c_search[0].apex()[2]);
  TriL<Type> tl1(c_search[1].apex()[1], c_search[1].apex()[2]);
  c_uncertain = {TriMap<Type>(tl0, c_ml), TriMap<Type>(tl1, c_ml)};
}

template <typename Type>
void FrTri<Type>::move_to(V2d<Type> p_src, std::vector<V2d<Type>> p_ml) {
  c_src.move_to(p_src);
  for (size_t i = 0; i < c_ml.size(); i++)
    c_ml[i].move_to(p_ml[i]);
  auto sq = include_all(c_ml, c_src);
  c_search = sq.enlarge(c_param.fact).to_tri();
  TriL<Type> tl0(c_search[0].apex()[1], c_search[0].apex()[2]);
  TriL<Type> tl1(c_search[1].apex()[1], c_search[1].apex()[2]);
  c_uncertain = {TriMap<Type>(tl0, c_ml), TriMap<Type>(tl1, c_ml)};
}

template <typename Type> void FrTri<Type>::split_u(int n, Mlens<Type> &ml) {
  c_uncertain.split(n, ml);
}

template <typename Type> Square<Type> FrTri<Type>::search() const {
  V2d<Type> cent = (c_search[0].apex()[0] + c_search[1].apex()[0]) * 0.5;
  Type d = std::abs(cent.x() - c_search[0].apex()[0].x());
  return Square<Type>(cent, d);
}

template <class Type> int FrTri<Type>::round() const { return c_round; }

template <typename Type> Type FrTri<Type>::d() const {
  return c_uncertain[0].tril().d();
}

template <class Type> TMQue<Type> FrTri<Type>::in() const { return c_in; }

template <class Type> TMQue<Type> FrTri<Type>::out() const { return c_out; }

template <class Type> TMQue<Type> FrTri<Type>::uncertain() const {
  return c_uncertain;
}

template <class Type> FrParam<Type> FrTri<Type>::param() const {
  return c_param;
}

template <typename Type> Mlens<Type> FrTri<Type>::ml() const { return c_ml; }

template <typename Type> Source<Type> FrTri<Type>::src() const { return c_src; }

template <typename Type> Square<Type> FrTri<Type>::include() const {
  return include_all(c_ml, c_src);
}

template <class Type> Type FrTri<Type>::area_in() const { return c_in.area(); }

template <class Type> Type FrTri<Type>::area_uncertain() const {
  return c_uncertain.area();
}

template <class Type> Type FrTri<Type>::area_iu() const {
  return this->area_iu() + this->area_uncertain() / 2.0;
}

template <typename Type> Type FrTri<Type>::area_sn() const {
  return c_in.area() + c_uncertain.area_sn(c_src);
}

template <typename Type> Brightness<Type> FrTri<Type>::brightness() const {
  Type br = std::accumulate(
      c_in.begin(), c_in.end(), 0.0, [&](Type acc, TriMap<Type> tm) {
        return acc + tm.tril().area() * tm.tris().brightness(c_src).mean;
      });
  br += std::accumulate(c_uncertain.begin(), c_uncertain.end(), 0.0,
                        [&](Type acc, TriMap<Type> tm) {
                          return acc + tm.tril().area() *
                                           tm.tris().brightness(c_src).mean;
                        });
  Type ms = std::accumulate(
      c_in.begin(), c_in.end(), 0.0, [&](Type acc, TriMap<Type> tm) {
        return acc + tm.tril().area() * tm.tris().brightness(c_src).rms *
                         tm.tril().area() * tm.tris().brightness(c_src).rms;
      });
  ms += std::accumulate(c_uncertain.begin(), c_uncertain.end(), 0.0,
                        [&](Type acc, TriMap<Type> tm) {
                          return acc + tm.tril().area() *
                                           tm.tris().brightness(c_src).rms *
                                           tm.tril().area() *
                                           tm.tris().brightness(c_src).rms;
                        });
  return {br, std::sqrt(ms)};
}

template <typename Type> Brightness<Type> FrTri<Type>::br_in() const {
  Type br = std::accumulate(
      c_in.begin(), c_in.end(), 0.0, [&](Type acc, TriMap<Type> tm) {
        return acc + tm.tril().area() * tm.tris().brightness(c_src).mean;
      });
  Type ms = std::accumulate(
      c_in.begin(), c_in.end(), 0.0, [&](Type acc, TriMap<Type> tm) {
        return acc + tm.tril().area() * tm.tris().brightness(c_src).rms *
                         tm.tril().area() * tm.tris().brightness(c_src).rms;
      });
  return {br, std::sqrt(ms)};
}

template <typename Type> Brightness<Type> FrTri<Type>::br_out() const {
  Type br = std::accumulate(
      c_out.begin(), c_out.end(), 0.0, [&](Type acc, TriMap<Type> tm) {
        return acc + tm.tril().area() * tm.tris().brightness(c_src).mean;
      });
  Type ms = std::accumulate(
      c_out.begin(), c_out.end(), 0.0, [&](Type acc, TriMap<Type> tm) {
        return acc + tm.tril().area() * tm.tris().brightness(c_src).rms *
                         tm.tril().area() * tm.tris().brightness(c_src).rms;
      });
  return {br, std::sqrt(ms)};
}

template <typename Type> Brightness<Type> FrTri<Type>::br_uncertain() const {
  Type br = std::accumulate(c_uncertain.begin(), c_uncertain.end(), 0.0,
                            [&](Type acc, TriMap<Type> tm) {
                              return acc + tm.tril().area() *
                                               tm.tris().brightness(c_src).mean;
                            });
  Type ms = std::accumulate(c_uncertain.begin(), c_uncertain.end(), 0.0,
                            [&](Type acc, TriMap<Type> tm) {
                              return acc + tm.tril().area() *
                                               tm.tris().brightness(c_src).rms *
                                               tm.tril().area() *
                                               tm.tris().brightness(c_src).rms;
                            });
  return {br, std::sqrt(ms)};
}

template <typename Type> V2d<Type> FrTri<Type>::sump_area() const {
  V2d<Type> smp = V2d_0<Type>;
  for (auto in : c_in)
    smp += in.tril().cent() * in.tril().area();
  for (auto unc : c_uncertain)
    smp += unc.tril().cent() * unc.area_sn(c_src);
  return smp;
}

template <typename Type> V2d<Type> FrTri<Type>::sump_br() const {
  V2d<Type> smp = V2d_0<Type>;
  for (auto in : c_in)
    smp += in.tril().cent() * in.tris().brightness(c_src).mean;
  for (auto unc : c_uncertain)
    smp += unc.tril().cent() * unc.tris().brightness(c_src).mean;
  return smp;
}

template <typename Type> V2d<Type> FrTri<Type>::ic_area() const {
  V2d<Type> cnt = V2d_0<Type>;
  for (auto in : c_in)
    cnt += in.tril().cent() * in.tril().area();
  for (auto unc : c_uncertain)
    cnt += unc.tril().cent() * unc.area_sn(c_src);
  return cnt / (c_in.area() + c_uncertain.area_sn(c_src));
}

template <typename Type> V2d<Type> FrTri<Type>::ic_br() const {
  V2d<Type> cnt = V2d_0<Type>;
  for (auto in : c_in)
    cnt += in.tril().cent() * in.tris().brightness(c_src).mean;
  for (auto unc : c_uncertain)
    cnt += unc.tril().cent() * unc.tris().brightness(c_src).mean;
  Type br = 0.0;
  for (auto in : c_in)
    br += in.tris().brightness(c_src).mean;
  for (auto unc : c_uncertain)
    br += unc.tris().brightness(c_src).mean;
  return cnt / br;
}

template <typename Type> V2d<Type> FrTri<Type>::as_area() const {
  return this->ic_area() - c_src.cent();
}

template <typename Type> V2d<Type> FrTri<Type>::as_br() const {
  return this->ic_br() - c_src.cent();
}

template <class Type> EnProcess FrTri<Type>::process_l() {
  if ((c_src.cent() - c_ml.cm()).norm() > c_param.sla_minbt) {
    Type d = c_src.circ().r();
    if (c_src.shape() == Shape::Elliptical)
      Type d = c_src.ell().ab().x();
    if (c_uncertain.front().tril().in(single_image(c_ml, c_src)[0]) &&
        c_uncertain.front().tril().d() > d * c_param.sla_fact_out) {
      std::array<TriMap<Type>, 2> vsp = c_uncertain.front().split(c_ml);
      c_uncertain.push_back(vsp[0]);
      c_uncertain.push_back(vsp[1]);
      return EnProcess::Break;
    }
  }
  for (size_t i = 0; i < c_ml.size(); i++) {
    switch (c_uncertain.front().tril().zone(c_ml[i].circ(), c_param.mrg_l)) {
    case Zone::uncertain:
      if (c_uncertain.front().tril().too_large(c_ml[i].circ(),
                                               c_param.max_drl)) {
        std::array<TriMap<Type>, 2> vsp = c_uncertain.front().split(c_ml);
        // for (TriMap<Type> sp : vsp)
        //   c_uncertain.push_back(sp);
        c_uncertain.push_back(vsp[0]);
        c_uncertain.push_back(vsp[1]);
        return EnProcess::Break;
      }
      break;
    case Zone::in:
      if (c_uncertain.front().tril().too_large(c_ml[i].circ(),
                                               c_param.max_drl)) {
        std::array<TriMap<Type>, 2> vsp = c_uncertain.front().split(c_ml);
        // for (TriMap<Type> sp : vsp)
        //   c_uncertain.push_back(sp);
        c_uncertain.push_back(vsp[0]);
        c_uncertain.push_back(vsp[1]);
        return EnProcess::Break;
      }
      if (c_src.shape() == Shape::Elliptical) {
        if (c_uncertain.front().tril().in(c_ml[i].p())) {
          if (c_uncertain.front().tris().in(c_src.ell()), c_param.in_fact) {
            if (c_param.acc_out)
              c_out.push_back(c_uncertain.front());
            return EnProcess::Break;
          } else {
            std::array<TriMap<Type>, 2> vsp = c_uncertain.front().split(c_ml);
            c_uncertain.push_back(vsp[0]);
            c_uncertain.push_back(vsp[1]);
            return EnProcess::Break;
          }
        }
      } else {
        if (c_uncertain.front().tril().in(c_ml[i].p())) {
          if (c_uncertain.front().tris().in(c_src.circ()), c_param.in_fact) {
            if (c_param.acc_out)
              c_out.push_back(c_uncertain.front());
            return EnProcess::Break;
          } else {
            std::array<TriMap<Type>, 2> vsp = c_uncertain.front().split(c_ml);
            c_uncertain.push_back(vsp[0]);
            c_uncertain.push_back(vsp[1]);
            return EnProcess::Break;
          }
        }
      }
    default:
      break;
    }
    return EnProcess::Continue;
  }
  return EnProcess::Continue;
}

template <class Type> void FrTri<Type>::process_s() {
  // switch (c_uncertain.front().tris().zone(src.circ(), c_param.max_drs,
  //                                         c_param.mxdr)) {
  // switch (c_uncertain.front().zone(c_src, c_ml)) {
  switch (c_uncertain.front().zone(c_src, c_param.mrg_s, c_ml)) {
  case Zone::in:
    c_in.push_back(c_uncertain.front());
    break;
  case Zone::out:
    if (c_param.acc_out)
      c_out.push_back(c_uncertain.front());
    break;
  case Zone::uncertain:
    std::array<TriMap<Type>, 2> vsp = c_uncertain.front().split(c_ml);
    // for (TriMap<Type> sp : vsp)
    //   c_uncertain.push_back(sp);
    c_uncertain.push_back(vsp[0]);
    c_uncertain.push_back(vsp[1]);
    break;
  }
}

template <class Type> void FrTri<Type>::process_front() {
  if (this->process_l() == EnProcess::Continue)
    this->process_s();
  c_uncertain.pop_front();
}

template <class Type> void FrTri<Type>::process() {
  // Type di = c_uncertain.front().tril().d();
  size_t sz = c_uncertain.size();
  // while (std::abs(c_uncertain.front().tril().d() - di) <=
  //        std::numeric_limits<Type>::epsilon())
  //   this->process_front();
  for (size_t i = 0; i < sz; i++)
    this->process_front();
  c_round++;
}

template <class Type> void FrTri<Type>::process(int n) {
  if (c_param.bench)
    c_param.begin = clock();
  for (int i = 0; i < n; i++)
    this->process();
  if (c_param.bench)
    c_param.end = clock();
}

template <class Type> void FrTri<Type>::process_while(int stop_in) {
  if (c_param.bench)
    c_param.begin = clock();
  while (c_in.size() < stop_in)
    this->process();
  if (c_param.bench)
    c_param.end = clock();
}

template <class Type> void FrTri<Type>::process_while(Type ratio_unc_in) {
  if (c_param.bench)
    c_param.begin = clock();
  while (c_uncertain.area() / c_in.area() > ratio_unc_in)
    this->process();
  if (c_param.bench)
    c_param.end = clock();
}

template <class Type>
void FrTri<Type>::process(Type prec, Mlens<Type> &ml, Source<Type> &src) {
  if (c_param.bench)
    c_param.begin = clock();
  Type una = std::numeric_limits<Type>::max();
  while (una > prec) {
    this->process();
    una = c_uncertain.area() / c_in.area();
  }
  if (c_param.bench)
    c_param.end = clock();
}

template <typename Type> void FrTri<Type>::clear() {
  c_in.clear();
  c_out.clear();
  c_uncertain.clear();
}

template <class Type> void FrTri<Type>::print() const {
  std::cout << "Area to search images :" << std::endl;
  this->search().print();
  std::cout << "size : in        = " << c_in.size() << std::endl;
  std::cout << "size : out       = " << c_out.size() << std::endl;
  std::cout << "size : uncertain = " << c_uncertain.size() << std::endl;
  if (c_uncertain.size() > 0)
    std::cout << "uncertain[0].tril().d() = " << c_uncertain[0].tril().d()
              << std::endl;
  std::cout << "FrParam :" << std::endl;
  std::cout << "mrg_l = " << c_param.mrg_l << ", mrg_s = " << c_param.mrg_s
            << std::endl;
  // std::cout << "max_drl = " << c_param.max_drl
  //           << ", max_drs = " << c_param.max_drs << std::endl;
  std::cout << "max_drl = " << c_param.max_drl << std::endl;
  std::cout << "fact = " << c_param.fact << std::endl;
  std::cout << "round = " << c_round << std::endl;
}

// /**
//  * template class Quad
//  * Type : float,  double,  long double
//  * Apex points of a square or a quadrilateral
//  * constructor Quad<Type>::Quad(V2d<Type> min, Type d)
//  * @param[in]  min : lower left apex of the square
//  * @param[in]  d : width of the square
//  *
//  */
//
// /**
//  * constructor Quad<Type>::Quad(Corner<Type> cn)
//  * @param[in] cn : four corner points
//  *
//  */
//
// /**
//  * method Quad<Type>::d()
//  * returns width of the square (valid only when Quad is
//  * a "right" square (on the lens plane)
//  * @return Type d() : width
//  *
//  */
// template <typename Type> Type Quad<Type>::d() const {
//   return c_corner[1].x() - c_corner[0].x();
// }
//
// /**
//  * method Quad<Type>::area()
//  * returns area of the square (valid only when Quad is
//  * a "right" square (on the lens plane)
//  * @return Type area() : area
//  *
//  */
// template <typename Type> Type Quad<Type>::area() const {
//   return this->d() * this->d();
// }
//
// template <typename Type> bool Quad<Type>::adjoin(Quad<Type> qd) const {
//   for (V2d<Type> p : qd.corner())
//     if (p.between(c_corner[0], c_corner[1]) ||
//         p.between(c_corner[1], c_corner[2]) ||
//         p.between(c_corner[2], c_corner[3]) ||
//         p.between(c_corner[3], c_corner[0]))
//       return true;
//   for (V2d<Type> p : c_corner)
//     if (p.between(qd.corner()[0], qd.corner()[1]) ||
//         p.between(qd.corner()[1], qd.corner()[2]) ||
//         p.between(qd.corner()[2], qd.corner()[3]) ||
//         p.between(qd.corner()[3], qd.corner()[0]))
//       return true;
//   return false;
// }
//
// /**
//  * method Quad<Type>::zone(Circle<Type> cl)
//  * returns zone of the points in the circle:
//  * Zone::in(i.e, securely in),  Zone::out(i.e,  securely out),
//  * or Zone::uncertain
//  * with default margin : 0.5
//  * @param[in] cl: a circle
//  * @return Zone zone() : zone of the points
//  *
//  */
// template <typename Type> Zone Quad<Type>::zone(Circle<Type> cl) const {
//   return this->zone(cl, 0.5);
// }
//
// /**
//  * method Quad<Type>::zone(Circle<Type> cl, Type mrg)
//  * returns zone of the points in the circle:
//  * Zone::in(i.e, securely in),  Zone::out(i.e,  securely out),
//  * or Zone::uncertain
//  * @param[in] cl: a circle
//  * @param[in] mrg : margin
//  * @return Zone zone() : zone of the points
//  *
//  */
// template <typename Type>
// Zone Quad<Type>::zone(Circle<Type> cl, Type mrg) const {
//   V2d<Type> vc(0.0, 0.0);
//   for (auto v : c_corner)
//     vc += v - cl.cent();
//
//   std::vector<Type> vt;
//   for (auto v : c_corner)
//     vt.push_back((v - cl.cent()) * (vc / vc.norm()));
//   Type mn = *std::min_element(vt.begin(), vt.end());
//   Type mx = *std::max_element(vt.begin(), vt.end());
//   std::vector<Type> vr;
//   for (auto v : c_corner)
//     vr.push_back((v - cl.cent()).norm());
//   Type mnr = *std::min_element(vr.begin(), vr.end());
//   Type mxr = *std::max_element(vr.begin(), vr.end());
//   if (mxr < (cl.r() - (mxr - mnr) * mrg) && (mxr - mnr) < 0.3 * cl.r()) {
//     return Zone::in;
//     //       if(this->in(cl)) return Zone::out;
//     //       else return Zone::in;
//   }
//   if (mn > (cl.r() + (mx - mn) * mrg * (mx - mn) / cl.r())) {
//     return Zone::out;
//     //        if(this->in(cl)) return Zone::out;
//     //       else return Zone::in;
//   }
//   return Zone::uncertain;
// }
//
// /**
//  * method Quad<Type>::zone(Ellipse<Type> el)
//  * returns zone of the points in the Ellipse:
//  * Zone::in(i.e, securely in),  Zone::out(i.e,  securely out),
//  * or Zone::uncertain
//  * @param[in] el: an ellipse
//  * @return Zone zone() : zone of the points
//  *
//  */
// template <typename Type> Zone Quad<Type>::zone(Ellipse<Type> el) const {
//   Corner<Type> cn2 = Corner_0<Type>;
//   std::transform(c_corner.begin(), c_corner.end(), cn2.begin(),
//                  [&el](V2d<Type> cn) {
//                    return (cn - el.cent()) / el.rho(cn - el.cent());
//                  });
//   Quad<Type> q2(cn2);
//   return q2.zone(Circle<Type>(V2d<Type>(0.0, 0.0), 1.0));
// }
//
// /**
//  * method Quad<Type>::zone(Ellipse<Type> el)
//  * returns zone of the points in the Ellipse in margin mrg:
//  * Zone::in(i.e, securely in),  Zone::out(i.e,  securely out),
//  * or Zone::uncertain
//  * @param[in] el: an ellipse
//  * @param[in] mrg : margin
//  * @return Zone zone() : zone of the points
//  *
//  */
// template <typename Type>
// Zone Quad<Type>::zone(Ellipse<Type> el, Type mrg) const {
//   Corner<Type> cn2 = Corner_0<Type>;
//   std::transform(c_corner.begin(), c_corner.end(), cn2.begin(),
//                  [&el](V2d<Type> cn) {
//                    return (cn - el.cent()) / el.rho(cn - el.cent());
//                  });
//   Quad<Type> q2(cn2);
//   return q2.zone(Circle<Type>(V2d<Type>(0.0, 0.0), 1.0), mrg);
// }
//
// template <typename Type> bool Quad<Type>::border(Circle<Type> cl) const {
//   if (cl.in(c_corner[0]) && cl.in(c_corner[1]) && cl.in(c_corner[2]) &&
//       cl.in(c_corner[3]))
//     return false;
//   if (!cl.in(c_corner[0]) && !cl.in(c_corner[1]) && !cl.in(c_corner[2]) &&
//       !cl.in(c_corner[3]))
//     return false;
//   return true;
// }
//
// /**
//  * method Quad<Type>::too_large(Circle<Type> cl)
//  * returns if the quadrilateral is too large compared to the circle
//  * with default factor 0.4
//  * @param[in] cl: a circle
//  * @return bool too_large() : true or false
//  *
//  */
// template <typename Type> bool Quad<Type>::too_large(Circle<Type> cl) const {
//   return this->too_large(cl, 0.4);
// }
//
// /**
//  * method Quad<Type>::too_large(Circle<Type> cl, Type fact)
//  * returns if the quadrilateral is too large compared to the circle
//  * @param[in] cl: a circle
//  * @param[in] fact: a factor to identify
//  * @return bool too_large() : true or false
//  *
//  */
// template <typename Type>
// bool Quad<Type>::too_large(Circle<Type> cl, Type fact) const {
//   return std::abs(this->d()) > cl.r() * fact;
// }
//
// template <typename Type> bool Quad<Type>::in(V2d<Type> p) const {
//   return (c_corner[1] - c_corner[0]) % (p - c_corner[0]) > 0.0 &&
//          (c_corner[2] - c_corner[1]) % (p - c_corner[1]) > 0.0 &&
//          (c_corner[3] - c_corner[2]) % (p - c_corner[2]) > 0.0 &&
//          (c_corner[0] - c_corner[3]) % (p - c_corner[3]) > 0.0;
// }
//
// template <typename Type> bool Quad<Type>::in(Circle<Type> cl, Type mrg) const
// {
//   if (this->in(cl.cent()))
//     return (cl.cent() - (c_corner[1] - c_corner[0]).unit()).norm() >
//                cl.r() * mrg &&
//            (cl.cent() - (c_corner[2] - c_corner[1]).unit()).norm() >
//                cl.r() * mrg &&
//            (cl.cent() - (c_corner[3] - c_corner[2]).unit()).norm() >
//                cl.r() * mrg &&
//            (cl.cent() - (c_corner[0] - c_corner[3]).unit()).norm() >
//                cl.r() * mrg;
//   else
//     return (c_corner[0] - cl.cent()).norm() < cl.r() / mrg &&
//            (c_corner[1] - cl.cent()).norm() < cl.r() / mrg &&
//            (c_corner[2] - cl.cent()).norm() < cl.r() / mrg &&
//            (c_corner[3] - cl.cent()).norm() < cl.r() / mrg;
//   return false;
// }
//
// template <typename Type> bool Quad<Type>::in(Circle<Type> cl) const {
//   return this->in(cl, 1.5);
// }
//
// /**
//  * method Quad<Type>::corner()
//  * returns corner points
//  * @return Corner<Type> corner() : corner points
//  *
//  */
// template <typename Type> Corner<Type> Quad<Type>::corner() const {
//   return c_corner;
// }
//
// template <typename Type> std::vector<Type> Quad<Type>::angles() const {
//   std::vector<Type> ang;
//   ang.push_back(std::asin((c_corner[1] - c_corner[0]).unit() %
//                           (c_corner[2] - c_corner[1]).unit()));
//   ang.push_back(std::asin((c_corner[2] - c_corner[1]).unit() %
//                           (c_corner[3] - c_corner[2]).unit()));
//   ang.push_back(std::asin((c_corner[3] - c_corner[2]).unit() %
//                           (c_corner[0] - c_corner[3]).unit()));
//   ang.push_back(std::asin((c_corner[0] - c_corner[3]).unit() %
//                           (c_corner[1] - c_corner[0]).unit()));
//   return ang;
// }
//
// template <typename Type> bool Quad<Type>::all_p() const {
//   std::vector<Type> ang = this->angles();
//   return ang[0] >= 0.0 && ang[1] >= 0.0 && ang[2] >= 0.0 && ang[3] >= 0.0;
// }
//
// template <typename Type> bool Quad<Type>::all_n() const {
//   std::vector<Type> ang = this->angles();
//   return ang[0] < 0.0 && ang[1] < 0.0 && ang[2] < 0.0 && ang[3] < 0.0;
// }
//
// template <typename Type> bool Quad<Type>::mixed() const {
//   return !this->all_p() && !this->all_n();
// }
//
// template <typename Type> V2d<Type> Quad<Type>::min() const {
//   return V2d<Type>(
//       std::min_element(c_corner.begin(), c_corner.end(),
//                        [](auto a, auto b) { return a.x() < b.x(); })[0]
//           .x(),
//       std::min_element(c_corner.begin(), c_corner.end(),
//                        [](auto a, auto b) { return a.y() < b.y(); })[0]
//           .y());
// }
//
// template <typename Type> V2d<Type> Quad<Type>::max() const {
//   return V2d<Type>(
//       std::max_element(c_corner.begin(), c_corner.end(),
//                        [](auto a, auto b) { return a.x() < b.x(); })[0]
//           .x(),
//       std::max_element(c_corner.begin(), c_corner.end(),
//                        [](auto a, auto b) { return a.y() < b.y(); })[0]
//           .y());
// }
//
// /**
//  * method Quad<Type>::newpl()
//  * returns  newly produced inner points to split
//  * on the lensing plane (valid only for lensing plane)
//  * @return Inner<Type> newpl() : inner points
//  *
//  */
// template <typename Type> Inner<Type> Quad<Type>::newpl() const {
//   Type d1 = this->d();
//   Type d2 = this->d() * 0.5;
//   return Inner<Type>{
//       c_corner[0] + V2d<Type>(d2, 0.0), c_corner[0] + V2d<Type>(d1, d2),
//       c_corner[0] + V2d<Type>(d2, d1), c_corner[0] + V2d<Type>(0.0, d2),
//       c_corner[0] + V2d<Type>(d2, d2)};
// }
//
// /**
//  * method Quad<Type>::newps(Inner<Type> npl, Mlens<Type> ml)
//  * returns calculated newly produced inner points to split
//  * on the source plane (valid only for source plane)
//  * @param[in] npl : inner points on the lensing plane
//  * @param[in] ml : lens system
//  * @return Inner<Type> newps() : corresponding inner points
//  * on the source plane
//  *
//  */
// template <typename Type>
// Inner<Type> Quad<Type>::newps(Inner<Type> npl, Mlens<Type> ml) const {
//   Inner<Type> nps = Inner_0<Type>;
//   for (int i = 0; i < static_cast<int>(npl.size()); i++)
//     nps[i] = ml.beta(npl[i]);
//   return nps;
// }
//
// /**
//  * method Quad<Type>::split()
//  * returns newly produced four squares
//  * on the lensing plane (valid only for lensing plane)
//  * @return std::vector<Quad<Type>> split() : produced
//  * four squares
//  *
//  */
// template <typename Type> std::vector<Quad<Type>> Quad<Type>::split() const {
//   Type d2 = this->d() * 0.5;
//   Inner<Type> newl = this->newpl();
//   std::vector<Quad<Type>> vq{Quad<Type>(c_corner[0], d2),
//                              Quad<Type>(newl[0], d2), Quad<Type>(newl[4],
//                              d2), Quad<Type>(newl[3], d2)};
//   return vq;
// }
//
// /**
//  * method Quad<Type>::split(Inner<Type> newps)
//  * returns newly produced four quadrilaterals
//  * on the source plane (valid only for source plane)
//  * @param[in] newps : inner points
//  * @return std::vector<Quad<Type>> split() : produced
//  * four squares
//  *
//  */
// template <typename Type>
// std::vector<Quad<Type>> Quad<Type>::split(Inner<Type> newps) const {
//   std::vector<Quad<Type>> vq{
//       Quad<Type>({c_corner[0], newps[0], newps[4], newps[3]}),
//       Quad<Type>({newps[0], c_corner[1], newps[1], newps[4]}),
//       Quad<Type>({newps[4], newps[1], c_corner[2], newps[2]}),
//       Quad({newps[3], newps[4], newps[2], c_corner[3]})};
//   return vq;
// }
//
// /**
//  * method Quad<Type>::tril()
//  * returns new two triangles
//  * on the lens plane
//  * @return std::vector<Quad<Type>> split() : produced
//  * four squares
//  *
//  */
// template <typename Type> std::array<TriL<Type>, 2> Quad<Type>::tril() const {
//   return {TriL(c_corner[1], c_corner[2]), TriL(c_corner[3], c_corner[0])};
// }
//
// /**
//  * method  Quad<Type>::to_f()
//  * returns converted Quad<Type> to Quad<float>
//  * @return Quad<float> Quad<Type>::to_f() :
//  * converted Quad<float>
//  *
//  */
// template <typename Type> Quad<float> Quad<Type>::to_f() const {
//   return Quad<float>({c_corner[0].to_f(), c_corner[1].to_f(),
//                       c_corner[2].to_f(), c_corner[3].to_f()});
// }
//
// /**
//  * method  Quad<Type>::to_d()
//  * returns converted Quad<Type> to Quad<double>
//  * @return Quad<float> Quad<Type>::to_d() :
//  * converted Quad<double>
//  *
//  */
// template <typename Type> Quad<double> Quad<Type>::to_d() const {
//   return Quad<double>({c_corner[0].to_d(), c_corner[1].to_d(),
//                        c_corner[2].to_d(), c_corner[3].to_d()});
// }
//
// /**
//  * method  Quad<Type>::to_l()
//  * returns converted Quad<Type> to Quad<long double>
//  * @return Quad<float> Quad<Type>::to_l() :
//  * converted Quad<long double>
//  *
//  */
// template <typename Type> Quad<long double> Quad<Type>::to_l() const {
//   return Quad<long double>({c_corner[0].to_l(), c_corner[1].to_l(),
//                             c_corner[2].to_l(), c_corner[3].to_l()});
// }
//
// /**
//  * method  Quad<Type>::print()
//  * prints parameters (corner points)
//  * @return none
//  *
//  */
// template <typename Type> void Quad<Type>::print() const {
//   std::cout << "Quad : corner :" << std::endl;
//   for (int i = 0; i < static_cast<int>(c_corner.size()); i++) {
//     std::cout << "i = " << i << " : ";
//     c_corner[i].print();
//   }
// }
//
// //
// // template class QuadMap square corner points and
// // corresponding points on the source plane
// //
// // template <typename Type>
// // QuadMap<Type>::QuadMap(Quad<Type> ql, Mlens<Type> ml, Source<Type> src) {
// //   c_ql = ql;
// //   c_ml = ml;
// //   c_srcs.push_back(src);
// //   Corner<Type> av = Corner_0<Type>;
// //   for (int i = 0; i < static_cast<int>(c_ql.corner().size()); i++)
// //     av[i] = c_ml.beta(c_ql.corner()[i]);
// //   c_qs = Quad<Type>(av);
// // }
//
// // template <typename Type>
// // QuadMap<Type>::QuadMap(Quad<Type> ql, Mlens<Type> ml,
// //                        std::vector<Source<Type>> srcs) {
// //   c_ql = ql;
// //   c_ml = ml;
// //   c_srcs = srcs;
// //   Corner<Type> av = Corner_0<Type>;
// //   for (int i = 0; i < static_cast<int>(c_ql.corner().size()); i++)
// //     av[i] = c_ml.beta(c_ql.corner()[i]);
// //   c_qs = Quad<Type>(av);
// // }
//
// // template <typename Type>
// // QuadMap<Type>::QuadMap(Quad<Type> ql, Quad<Type> qs, Mlens<Type> ml,
// //                        Source<Type> src) {
// //   c_ql = ql;
// //   c_qs = qs;
// //   c_ml = ml;
// //   c_srcs.push_back(src);
// // }
//
// // template <typename Type>
// // QuadMap<Type>::QuadMap(Quad<Type> ql, Quad<Type> qs, Mlens<Type> ml,
// //                        std::vector<Source<Type>> srcs) {
// //   c_ql = ql;
// //   c_qs = qs;
// //   c_ml = ml;
// //   c_srcs = srcs;
// // }
//
// template <typename Type> Mlens<Type> QuadMap<Type>::ml() const { return c_ml;
// }
//
// template <typename Type> Source<Type> QuadMap<Type>::src() const {
//   return c_src;
// }
//
// template <typename Type> Quad<Type> QuadMap<Type>::ql() const { return c_ql;
// }
//
// template <typename Type> Quad<Type> QuadMap<Type>::qs() const { return c_qs;
// }
//
// template <typename Type> Zone QuadMap<Type>::zone_s() const {
//   switch (c_src.shape()) {
//   case Circular:
//     return c_qs.zone(c_src.circ());
//     break;
//   case Elliptical:
//     return c_qs.zone(c_src.ell());
//     break;
//   }
//   return Zone::out;
// }
//
// template <typename Type> Zone QuadMap<Type>::zone_s(Type mrg) const {
//   switch (c_src.shape()) {
//   case Circular:
//     if (c_qs.zone(c_src.circ(), mrg) == Zone::uncertain)
//       return Zone::uncertain;
//     else if (c_qs.zone(c_src.circ(), mrg) == Zone::in) {
//       for (size_t i = 0; i < c_ml.size(); i++) {
//         if (this->zone_l(static_cast<int>(i)) == Zone::in) {
//           if (this->ql().in(c_ml[i].p()))
//             return Zone::out;
//         }
//       }
//       return Zone::in;
//     } else {
//       for (size_t i = 0; i < c_ml.size(); i++) {
//         if (this->zone_l(static_cast<int>(i)) == Zone::in) {
//           if (this->ql().in(c_ml[i].p()))
//             return Zone::in;
//         }
//       }
//       return Zone::out;
//     }
//     return Zone::uncertain;
//     break;
//   case Elliptical:
//     return c_qs.zone(c_src.ell(), mrg);
//     break;
//   }
//   return Zone::out;
// }
//
// template <typename Type>
// std::vector<Zone> QuadMap<Type>::zone_l(Type mrg) const {
//   std::vector<Zone> vz;
//   for (int i = 0; i < static_cast<int>(c_ml.size()); i++)
//     vz.push_back(this->zone_l(i, mrg));
//   return vz;
// }
//
// template <typename Type> Zone QuadMap<Type>::zone_l(int i, Type mrg) const {
//   return c_ql.zone(c_ml[i].circ(), mrg);
// }
//
// template <typename Type> Zone QuadMap<Type>::zone_l(int i) const {
//   return this->zone_l(i, 0.2);
// }
//
// template <typename Type>
// bool QuadMap<Type>::too_large_l(int i, Type drl) const {
//   if (c_ql.d() > drl * c_ml[i].re())
//     return true;
//   return false;
// }
//
// template <typename Type> bool QuadMap<Type>::too_large_s(Type drs) const {
//   if (c_src.shape() == Circular) {
//     if (c_qs.too_large(c_src.circ()))
//       return true;
//   }
//   return false;
// }
//
// template <typename Type>
// std::vector<QuadMap<Type>> QuadMap<Type>::split() const {
//   std::vector<QuadMap<Type>> vqm;
//   std::vector<Quad<Type>> vql = c_ql.split();
//   std::vector<Quad<Type>> vqs = c_qs.split(c_qs.newps(c_ql.newpl(), c_ml));
//   for (size_t i = 0; i < vql.size(); i++)
//     vqm.push_back(QuadMap<Type>(vql[i], vqs[i], c_ml, c_src));
//   return vqm;
// }
//
// template <typename Type> QuadMap<float> QuadMap<Type>::to_f() const {
//   return QuadMap<float>(c_ql.to_f(), c_ml.to_f(), c_src.to_f());
// }
//
// template <typename Type> QuadMap<double> QuadMap<Type>::to_d() const {
//   return QuadMap<double>(c_ql.to_d(), c_ml.to_d(), c_src.to_d());
// }
//
// template <typename Type> QuadMap<long double> QuadMap<Type>::to_l() const {
//   return QuadMap<long double>(c_ql.to_l(), c_ml.to_l(), c_src.to_l());
// }
//
// template <typename Type> void QuadMap<Type>::print() const {
//   std::cout << "QuadMap :" << std::endl;
//   std::cout << "ml : " << std::endl;
//   c_ml.print();
//   std::cout << "src = " << std::endl;
//   c_src.print();
//   std::cout << "ql : " << std::endl;
//   c_ql.print();
//   std::cout << "qs : " << std::endl;
//   c_qs.print();
// }
//
// // template <typename Type> QMQue<Type>::QMQue(QuadMap<Type> qm) {
// //   this->push_back(qm);
// // }
//
// // template <typename Type> QMQue<Type>::QMQue(const QMQue<Type> &qmu) {
// //   *this = qmu;
// // }
// template <typename Type> QuadMap<Type> &QMQue<Type>::operator[](size_type i)
// {
//   return this->at(i);
// }
//
// template <typename Type>
// const QuadMap<Type> &QMQue<Type>::operator[](size_type i) const {
//   return this->at(i);
// }
//
// template <typename Type> Type QMQue<Type>::area() const {
//   Type ar = 0.0;
//   for (QuadMap<Type> qm : *this)
//     ar += qm.ql().area();
//   return ar;
// }
//
// template <typename Type> QMQue<Type> QMQue<Type>::all_p() const {
//   QMQue<Type> qmq;
//   for (size_t i = 0; i < this->size(); i++)
//     if (this->at(i).qs().all_p())
//       qmq.push_back(this->at(i));
//   return qmq;
// }
//
// template <typename Type> QMQue<Type> QMQue<Type>::all_n() const {
//   QMQue<Type> qmq;
//   for (size_t i = 0; i < this->size(); i++)
//     if (this->at(i).qs().all_n())
//       qmq.push_back(this->at(i));
//   return qmq;
// }
//
// template <typename Type> QMQue<Type> QMQue<Type>::mixed() const {
//   QMQue<Type> qmq;
//   for (size_t i = 0; i < this->size(); i++)
//     if (!this->at(i).qs().all_p() && !this->at(i).qs().all_n())
//       qmq.push_back(this->at(i));
//   return qmq;
// }
//
// template <typename Type> QMQue<Type> QMQue<Type>::border() const {
//   QMQue<Type> bd;
//   for (size_t i = 0; i < this->size(); i++)
//     if (this->at(i).qs().border(this->at(i).src().circ()))
//       bd.push_back(this->at(i));
//   return bd;
// }
//
// template <class Type> void QMQue<Type>::qmqp_init() {
//   for (QuadMap<Type> qm : *this)
//     c_qmqp.push(qm);
// }
//
// template <class Type> void QMQue<Type>::qmqp_pop() { c_qmqp.pop(); }
//
// template <class Type> void QMQue<Type>::qmqp_rotate() {
//   c_qmqp.push(c_qmqp.front());
//   c_qmqp.pop();
// }
//
// template <class Type> size_t QMQue<Type>::qmqp_size() const {
//   return c_qmqp.size();
// }
//
// template <class Type> QMQue<Type> QMQue<Type>::group() {
//   QMQue<Type> qmq;
//   qmq.push_back(c_qmqp.front());
//   c_qmqp.pop();
//   size_t iq = 0;
//   while (iq < qmq.size()) {
//     size_t nc = c_qmqp.size();
//     for (size_t ic = 0; ic < nc; ic++) {
//       if (qmq[iq].ql().adjoin(c_qmqp.front().ql())) {
//         qmq.push_back(c_qmqp.front());
//         c_qmqp.pop();
//       }
//       c_qmqp.push(c_qmqp.front());
//       c_qmqp.pop();
//     }
//     iq++;
//   }
//   return qmq;
// }
//
// template <class Type> std::vector<QMQue<Type>> QMQue<Type>::groups() {
//   std::vector<QMQue<Type>> grps;
//   this->qmqp_init();
//   QMQue<Type> grp;
//   while (c_qmqp.size() > 0) {
//     grp = this->group();
//     grps.push_back(grp);
//   }
//   return grps;
// }
//
// template <class Type> std::queue<QuadMap<Type>> QMQue<Type>::qmqp() {
//   std::queue<QuadMap<Type>> qmqp;
//   size_t nq = c_qmqp.size();
//   for (size_t i = 0; i < nq; i++) {
//     qmqp.push(c_qmqp.front());
//     c_qmqp.push(c_qmqp.front());
//     c_qmqp.pop();
//   }
//   return qmqp;
// }
//
// template <class Type> QMQue<Type> QMQue<Type>::qmqp_to_qmque() {
//   QMQue<Type> qmq;
//   size_t nq = c_qmqp.size();
//   for (size_t i = 0; i < nq; i++) {
//     qmq.push_back(c_qmqp.front());
//     c_qmqp.push(c_qmqp.front());
//     c_qmqp.pop();
//   }
//   return qmq;
// }
//
// template <class Type> QuadMap<Type> QMQue<Type>::qmqp_front() const {
//   return c_qmqp.front();
// }
//
// template <class Type> QuadMap<Type> QMQue<Type>::qmqp_back() const {
//   return c_qmqp.back();
// }
//
// template <typename Type> void QMQue<Type>::print() const {}
//
// // template <typename Type>
// // FrQuad<Type>::FrQuad(QMQue<Type> unc, FrParam<Type> param) {
// //   c_uncertain = unc;
// //   c_param = param;
// // }
//
// template <typename Type> QMQue<Type> FrQuad<Type>::in() const { return c_in;
// }
//
// template <typename Type> QMQue<Type> FrQuad<Type>::out() const { return
// c_out; }
//
// template <typename Type> QMQue<Type> FrQuad<Type>::uncertain() const {
//   return c_uncertain;
// }
//
// template <typename Type> FrParam<Type> FrQuad<Type>::param() const {
//   return c_param;
// }
//
// template <typename Type> Type FrQuad<Type>::area_in() const {
//   Type area = 0.0;
//   for (QuadMap<Type> in : c_in)
//     area += in.ql().area();
//   return area;
// }
//
// template <typename Type> Type FrQuad<Type>::area_uncertain() const {
//   Type area = 0.0;
//   for (QuadMap<Type> unc : c_uncertain)
//     area += unc.ql().area();
//   return area;
// }
//
// template <typename Type> EnProcess FrQuad<Type>::process_l() {
//   for (size_t i = 0; i < c_uncertain.front().ml().size(); i++) {
//     switch (c_uncertain.front().zone_l(static_cast<int>(i))) {
//     case Zone::uncertain:
//       if (c_uncertain.front().too_large_l(i, c_param.max_drl)) {
//         std::vector<QuadMap<Type>> vsp = c_uncertain.front().split();
//         for (QuadMap<Type> sp : vsp)
//           c_uncertain.push_back(sp);
//         c_uncertain.pop_front();
//         return EnProcess::Break;
//       }
//       break;
//     default:
//       return EnProcess::Continue;
//       break;
//     }
//   }
//   return EnProcess::Continue;
// }
//
// template <typename Type> void FrQuad<Type>::process_s() {
//   switch (c_uncertain.front().zone_s()) {
//   case Zone::in:
//     c_in.push_back(c_uncertain.front());
//     c_uncertain.pop_front();
//     break;
//   case Zone::out:
//     c_out.push_back(c_uncertain.front());
//     c_uncertain.pop_front();
//     break;
//   case Zone::uncertain:
//     std::vector<QuadMap<Type>> vsp = c_uncertain.front().split();
//     for (QuadMap<Type> sp : vsp)
//       c_uncertain.push_back(sp);
//     c_uncertain.pop_front();
//     break;
//   }
// }
//
// template <typename Type> void FrQuad<Type>::process_front() {
//   if (this->process_l() == EnProcess::Continue)
//     this->process_s();
// }
//
// template <typename Type> void FrQuad<Type>::process() {
//   Type di = c_uncertain.front().ql().d();
//   while (c_uncertain.front().ql().d() == di)
//     this->process_front();
// }
//
// template <typename Type> void FrQuad<Type>::process(int n) {
//   for (int i = 0; i < n; i++)
//     this->process_front();
// }
//
// template <typename Type> void FrQuad<Type>::process(Type prec) {
//   Type una = std::numeric_limits<Type>::max();
//   while (una > prec) {
//     this->process_front();
//     una = c_uncertain.area() / c_in.area();
//   }
// }
//
// template <typename Type> void FrQuad<Type>::print() const {}

//
// Exticit instantations for float
//
// template class Quad<float>;
// template class QuadMap<float>;
// template class QMQue<float>;
// template class FrQuad<float>;

template class TriL<float>;
template struct Brightness<float>;
template class TriS<float>;
template class TriMap<float>;
template class TMQue<float>;
template struct FrParam<float>;
template class Snap<float>;
template class FrTri<float>;

template std::array<V2d<float>, 2> single_image(Mlens<float> ml,
                                                Source<float> src);
template Square<float> include_all(Mlens<float> ml, Source<float> src);

//
//
// For double
//
// template class Quad<double>;
// template class QuadMap<double>;
// template class QMQue<double>;
// template class FrQuad<double>;

template class TriL<double>;
template struct Brightness<double>;
template class TriS<double>;
template class TriMap<double>;
template class TMQue<double>;
template struct FrParam<double>;
template class Snap<double>;
template class FrTri<double>;

template std::array<V2d<double>, 2> single_image(Mlens<double> ml,
                                                 Source<double> src);
template Square<double> include_all(Mlens<double> ml, Source<double> src);

//
// For long double
//
// template class Quad<long double>;
// template class QuadMap<long double>;
// template class QMQue<long double>;
// template class FrQuad<long double>;

template class TriL<long double>;
template struct Brightness<long double>;
template class TriS<long double>;
template class TriMap<long double>;
template class TMQue<long double>;
template struct FrParam<long double>;
template class Snap<long double>;
template class FrTri<long double>;

template std::array<V2d<long double>, 2> single_image(Mlens<long double> ml,
                                                      Source<long double> src);
template Square<long double> include_all(Mlens<long double> ml,
                                         Source<long double> src);
