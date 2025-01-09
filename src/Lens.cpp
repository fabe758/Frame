/*
 *   Lens.cpp
 *
 *
 *   Created by F. Abe
 *
 */

/**
 * @file Lens.cpp
 +
 * @brief Class library to handle lens system
 * @author Fumio Abe
 * @date 23 April,  2023
 *
 * @details Single lens and multiple lens system
 * @note
 *
 */

/**
 * include headers  for this project
 *
 */
#include "Lens.h"
#include "Geom.h"
#include <algorithm>

// template <class Type> Lens<Type>::Lens(){};

template <class Type> Type Lens<Type>::q() const { return c_q; }

template <class Type> Type Lens<Type>::re() const { return std::sqrt(c_q); }

template <class Type> Circle<Type> Lens<Type>::circ() const {
  return Circle<Type>(c_p, this->re());
}

template <class Type> V2d<Type> Lens<Type>::p() const { return c_p; }

template <typename Type> V2d<Type> Lens<Type>::max(Type fre) const {
  return V2d<Type>(c_p.x() + fre * re(), c_p.y() + fre * re());
}

template <typename Type> V2d<Type> Lens<Type>::min(Type fre) const {
  return V2d<Type>(c_p.x() - fre * re(), c_p.y() - fre * re());
}

template <typename Type> Rectangle<Type> Lens<Type>::rec(Type fre) const {
  return Rectangle<Type>(this->min(fre), this->max(fre));
}

template <typename Type> Square<Type> Lens<Type>::sq(Type fre) const {
  return Square<Type>(this->p(), this->re() * fre);
}

template <typename Type> Type Lens<Type>::g(V2d<Type> th, int alpha) const {
  if (th == this->p())
    return std::numeric_limits<Type>::max() * 0.01; // Ignore singularity
  else
    return std::pow(((th - this->p()).n2()), -alpha);
}

template <typename Type> Lens<float> Lens<Type>::to_f() const {
  return Lens<float>(static_cast<float>(c_q), c_p.to_f());
}

template <typename Type> Lens<double> Lens<Type>::to_d() const {
  return Lens<double>(static_cast<double>(c_q), c_p.to_d());
}

template <typename Type> Lens<long double> Lens<Type>::to_l() const {
  return Lens<long double>(static_cast<long double>(c_q), c_p.to_l());
}

template <class Type> void Lens<Type>::print() const {
  std::cout << "q = " << c_q << ", p : ";
  c_p.print();
}

template <typename Type> void Lens<Type>::move_to(V2d<Type> p) { c_p = p; }

template <class Type> bool Lens<Type>::operator==(Lens<Type> l2) const {
  return l2.q() == c_q && l2.p() == c_p;
}

template <typename Type> Lens<Type> &Mlens<Type>::operator[](size_type i) {
  return this->at(i);
}

//
// class Mlens
//
template <typename Type>
const Lens<Type> &Mlens<Type>::operator[](size_type i) const {
  return this->at(i);
}

template <typename Type> void Mlens<Type>::sort_i() {
  std::sort(this->begin(), this->end(),
            [](auto &a, auto &b) { return a.q() < b.q(); });
}

template <typename Type> void Mlens<Type>::sort_d() {
  std::sort(this->begin(), this->end(),
            [](auto &a, auto &b) { return a.q() > b.q(); });
}

template <typename Type> void Mlens<Type>::origin(V2d<Type> org) {
  for (int i = 0; i < static_cast<int>(this->size()); i++)
    (*this)[i] = Lens<Type>((*this)[i].q(), (*this)[i].p() - org);
}

template <typename Type> V2d<Type> Mlens<Type>::cm() const {
  V2d<Type> cent(0.0, 0.0);
  for (int i = 0; i < static_cast<int>(this->size()); i++)
    cent += (*this)[i].p() * (*this)[i].q();
  return cent / sumq();
}

template <typename Type> Type Mlens<Type>::sumq() const {
  return std::accumulate(this->begin(), this->end(), 0.0,
                         [](Type a, auto b) { return a + b.q(); });
}

template <typename Type> void Mlens<Type>::normq() {
  Type Q = sumq();
  for (int i = 0; i < static_cast<int>(this->size()); i++)
    (*this)[i] = Lens<Type>((*this)[i].q() / Q, (*this)[i].p());
}

template <typename Type> Lens<Type> Mlens<Type>::maxq() const {
  return *std::max_element(this->begin(), this->end(),
                           [](auto a, auto b) { return a.q() < b.q(); });
}

template <typename Type> Lens<Type> Mlens<Type>::minq() const {
  return *std::min_element(this->begin(), this->end(),
                           [](auto a, auto b) { return a.q() < b.q(); });
}

template <typename Type> Lens<Type> Mlens<Type>::closest(V2d<Type> p) const {
  return *std::min_element(this->begin(), this->end(), [&p](auto a, auto b) {
    return (a.p() - p).norm() < (b.p() - p).norm();
  });
}

template <typename Type> V2d<Type> Mlens<Type>::max(Type fre) const {
  Lens<Type> maxx =
      *std::max_element(this->begin(), this->end(), [&fre](auto &a, auto &b) {
        return a.max(fre).x() < b.max(fre).x();
      });
  Lens<Type> maxy =
      *std::max_element(this->begin(), this->end(), [&fre](auto &a, auto &b) {
        return a.max(fre).y() < b.max(fre).y();
      });
  return V2d<Type>(maxx.max(fre).x(), maxy.max(fre).y());
}

template <typename Type> V2d<Type> Mlens<Type>::min(Type fre) const {
  Lens<Type> minx =
      *std::min_element(this->begin(), this->end(), [&fre](auto &a, auto &b) {
        return a.min(fre).x() < b.min(fre).x();
      });
  Lens<Type> miny =
      *std::min_element(this->begin(), this->end(), [&fre](auto &a, auto &b) {
        return a.min(fre).y() < b.min(fre).y();
      });
  return V2d<Type>(minx.min(fre).x(), miny.min(fre).y());
}

template <typename Type> Square<Type> Mlens<Type>::sq(Type fre) {
  return Square<Type>((this->min(fre) + this->max(fre)) / 2.0,
                      (this->max(fre) - this->min(fre)).x() / 2.0);
}

template <typename Type> M2d<Type> Mlens<Type>::J(V2d<Type> th) const {
  M2d<Type> j(1.0, 0.0, 0.0, 1.0);
  for (int i = 0; i < static_cast<int>(this->size()); i++)
    j += M2d<Type>(-(*this)[i].g(th, 1) +
                       2.0 * std::pow(th.x() - (*this)[i].p().x(), 2) *
                           (*this)[i].g(th, 2),
                   2.0 * (th.x() - (*this)[i].p().x()) *
                       (th.y() - (*this)[i].p().y()) * (*this)[i].g(th, 2),
                   2.0 * (th.x() - (*this)[i].p().x()) *
                       (th.y() - (*this)[i].p().y()) * (*this)[i].g(th, 2),
                   -(*this)[i].g(th, 1) +
                       2.0 * std::pow(th.y() - (*this)[i].p().y(), 2) *
                           (*this)[i].g(th, 2)) *
         (*this)[i].q();
  return j;
}

template <typename Type> Type Mlens<Type>::mu(V2d<Type> th) const {
  return 1.0 / J(th).det();
}

template <typename Type> V2d<Type> Mlens<Type>::beta(V2d<Type> th) const {
  V2d<Type> bt = th;
  for (Lens<Type> lens : *this)
    bt -= (th - lens.p()) / (th - lens.p()).n2() * lens.q();
  return bt;
}

template <typename Type>
Type Mlens<Type>::db(V2d<Type> th, V2d<Type> bt) const {
  return (this->beta(th) - bt).norm();
}

template <typename Type> void Mlens<Type>::operator+=(Lens<Type> l) {
  this->push_back(l);
}

template <typename Type> Mlens<float> Mlens<Type>::to_f() const {
  Mlens<float> ml(this->size());
  std::transform(this->begin(), this->end(), ml.begin(),
                 [](auto x) { return x.to_f(); });
  return ml;
}

template <typename Type> Mlens<double> Mlens<Type>::to_d() const {
  Mlens<double> ml(this->size());
  std::transform(this->begin(), this->end(), ml.begin(),
                 [](auto x) { return x.to_d(); });
  return ml;
}

template <typename Type> Mlens<long double> Mlens<Type>::to_l() const {
  Mlens<long double> ml(this->size());
  std::transform(this->begin(), this->end(), ml.begin(),
                 [](auto x) { return x.to_l(); });
  return ml;
}

template <typename Type> void Mlens<Type>::print() const {
  std::cout << "size = " << this->size() << std::endl;
  for (int i = 0; i < static_cast<int>(this->size()); i++)
    this->at(i).print();
}

template <typename Type> void Mlens<Type>::move_to(std::vector<V2d<Type>> vp) {
  for (int i = 0; i < static_cast<int>(this->size()); i++)
    this->at(i).move_to(vp[i]);
}

//
// Functions
//

//
// Explicit instantation of templates of classes
//
// For float
//
template class Lens<float>;
template class Mlens<float>;

//
template class Lens<double>;
template class Mlens<double>;

//
// For long double
//
template class Lens<long double>;
template class Mlens<long double>;

//
// Explicit instantation of templates of Functions
//
//  For Lens conv
//

//
// For Mlens conv
//
