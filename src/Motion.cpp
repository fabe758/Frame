/**
 * @file Motion.cpp
 +
 * @brief Class library to handle motions of sources and lenses
 * @author Fumio Abe
 * @date 6 January,  2024
 *
 * @details class library for Linear, parallax, and other motions
 * @note
 *
 */

/**
 * include headers  for this project
 * all c++ headers are loaded in Geom.h
 * Geom.h,  Lens.h,  Source.h, and Fractal.h are loaded in Motion.h
 *
 */
#include "Motion.h"

template <typename Type> V2d<Type> Linear<Type>::p0() const { return c_p0; };

template <typename Type> Time<Type> Linear<Type>::t0() const { return c_t0; };

template <typename Type> V2d<Type> Linear<Type>::v() const { return c_v; };

template <typename Type> Type Linear<Type>::u0() const { return c_p0.norm(); };

template <typename Type> Angle<Type> Linear<Type>::phi() const {
  return std::acos(c_v.x() / c_v.norm());
};

template <typename Type> Time<Type> Linear<Type>::te() const {
  return Time<Type>(1.0 / c_v.norm());
};

template <typename Type> V2d<Type> Linear<Type>::p(Time<Type> t) const {
  return c_p0 + c_v * (t - c_t0).day();
};

template <typename Type> Type Linear<Type>::u(Time<Type> t) const {
  return std::sqrt(this->p(t).x() * this->p(t).x() +
                   this->p(t).y() * this->p(t).y());
};

template <typename Type> void Linear<Type>::print() const {
  std::cout << "p0 :" << std::endl;
  c_p0.print();
  std::cout << "t0 = " << c_t0.day() << std::endl;
  std::cout << "v :" << std::endl;
  c_v.print();
};

template <typename Type> Angle<Type> Kepler<Type>::apa_i(Time<Type> t) const {
  Type l = this->ma(t).rad();
  Type e = c_obt.e;
  Type u = l;
  u += (e - e * e * e / 8.0 + e * e * e * e * e / 192.0 -
        e * e * e * e * e * e * e / 9216.0) *
       std::sin(l);
  u += (e * e / 2.0 - e * e * e * e / 6.0 + e * e * e * e * e * e / 48.0) *
       std::sin(2.0 * l);
  u += (3.0 * e * e * e / 8.0 - 27.0 * e * e * e * e * e / 128.0 +
        243.0 * e * e * e * e * e * e * e / 5120.0) *
       std::sin(3.0 * l);
  u += (e * e * e * e / 3.0 - 4.0 * e * e * e * e * e * e / 15.0) *
       std::sin(4.0 * l);
  u += (125.0 * e * e * e * e * e / 384.0 -
        3125.0 * e * e * e * e * e * e * e / 9216.0) *
       std::sin(5.0 * l);
  u += 27.0 * e * e * e * e * e * e * e / 80.0 * std::sin(6.0 * l);
  u += 16807 * e * e * e * e * e * e * e / 46080 * std::sin(7.0 * l);
  return u;
};

template <typename Type> Kepler<Type>::Kepler(std::string obj) {
  Time<Type> T;
  Time<Type> tp;
  if (obj == "Moon") {
  } else {
    T.from_year(1.000017);
    Time<Type> tp = JD245<Type>(2024, 01, 03, 00, 39, 0.0);
    c_obt = {V2d_0<Type>,      Length<Type>(1.00), 0.01670,         T, tp,
             Angle<Type>(0.0), Angle<Type>(0.0),   Angle<Type>(0.0)};
  };
};

template <typename Type> V2d<Type> SrcMotion<Type>::p(Time<Type> t) const {
  return c_f(t);
}

template <typename Type>
Type SrcMotion<Type>::u(Time<Type> t, V2d<Type> p0) const {
  return (this->p(t) - p0).norm();
};

template <typename Type> Type SrcMotion<Type>::u(Time<Type> t) const {
  return this->p(t).norm();
};

template <typename Type>
std::vector<V2d<Type>> SrcMotion<Type>::p(std::vector<Time<Type>> t) const {
  std::vector<V2d<Type>> p;
  for (size_t i = 0; i < t.size(); i++)
    p.push_back(this->p(t[i]));
  return p;
};

template <typename Type>
std::vector<Type> SrcMotion<Type>::u(std::vector<Time<Type>> t,
                                     V2d<Type> p0) const {
  std::vector<Type> u;
  for (size_t i = 0; i < t.size(); i++)
    u.push_back(this->u(t[i], p0));
  return u;
};

template <typename Type>
std::vector<Type> SrcMotion<Type>::u(std::vector<Time<Type>> t) const {
  std::vector<Type> u;
  for (size_t i = 0; i < t.size(); i++)
    u.push_back(this->u(t[i]));
  return u;
};

template <typename Type> SrcMotion<float> SrcMotion<Type>::to_f() const {
  return SrcMotion<float>(c_src.to_f(), [this](Time<float> t) {
    auto p = this->c_f(Time<Type>(static_cast<Type>(t.day())));
    return p.to_f();
  });
};

template <typename Type> SrcMotion<double> SrcMotion<Type>::to_d() const {
  return SrcMotion<double>(c_src.to_d(), [this](Time<double> t) {
    auto p = this->c_f(Time<Type>(static_cast<Type>(t.day())));
    return p.to_d();
  });
};

template <typename Type> SrcMotion<long double> SrcMotion<Type>::to_l() const {
  return SrcMotion<long double>(c_src.to_l(), [this](Time<long double> t) {
    auto p = this->c_f(Time<Type>(static_cast<Type>(t.day())));
    return p.to_l();
  });
};

template <typename Type> void SrcMotion<Type>::print() const {
  std::cout << "Source : " << std::endl;
  c_src.print();
};

template <typename Type> LensMotion<float> LensMotion<Type>::to_f() const {
  return LensMotion<float>(c_l.to_f(), [this](Time<float> t) {
    auto p = this->c_f(Time<Type>(static_cast<Type>(t.day())));
    return p.to_f();
  });
};

template <typename Type> LensMotion<double> LensMotion<Type>::to_d() const {
  return LensMotion<double>(c_l.to_d(), [this](Time<double> t) {
    auto p = this->c_f(Time<Type>(static_cast<Type>(t.day())));
    return p.to_d();
  });
};

template <typename Type>
LensMotion<long double> LensMotion<Type>::to_l() const {
  return LensMotion<long double>(c_l.to_l(), [this](Time<long double> t) {
    auto p = this->c_f(Time<Type>(static_cast<Type>(t.day())));
    return p.to_l();
  });
};

template <typename Type> void LensMotion<Type>::print() const {
  std::cout << "Lens : " << std::endl;
  c_l.print();
};

template <typename Type>
LensMotion<Type> &MlMotion<Type>::operator[](size_type i) {
  return this->at(i);
}

template <typename Type>
const LensMotion<Type> &MlMotion<Type>::operator[](size_type i) const {
  return this->at(i);
}

//
// Exticit instantations for float
//
template class Linear<float>;
template class Length<float>;
template class V2dU<float>;
template class V3dU<float>;
template class Time<float>;
template class JD245<float>;
template class EllipseU<float>;
template struct Orbit<float>;
template class Kepler<float>;
// template class ObtEarth<float>;
template class PiE<float>;
template class Binary<float>;
template class SrcMotion<float>;
template class LensMotion<float>;
template class MlMotion<float>;

//
// For double
//
template class Linear<double>;
template class Length<double>;
template class V2dU<double>;
template class V3dU<double>;
template class Time<double>;
template class JD245<double>;
template class EllipseU<double>;
template struct Orbit<double>;
template class Kepler<double>;
// template class ObtEarth<double>;
template class PiE<double>;
template class Binary<double>;
template class SrcMotion<double>;
template class LensMotion<double>;
template class MlMotion<double>;

//
// For long double
//
template class Linear<long double>;
template class Length<long double>;
template class V2dU<long double>;
template class V3dU<long double>;
template class Time<long double>;
template class JD245<long double>;
template class EllipseU<long double>;
template struct Orbit<long double>;
template class Kepler<long double>;
// template class ObtEarth<long double>;
template class PiE<long double>;
template class Binary<long double>;
template class SrcMotion<long double>;
template class LensMotion<long double>;
template class MlMotion<long double>;
