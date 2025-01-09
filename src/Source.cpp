/*
 *   Source.cpp
 *
 *
 *   Created by F. Abe on 08/12/21.
 *
 */

#include "Source.h"
#include <algorithm>

std::string shmap(enum Shape sh) {
  if (sh == Shape::Elliptical)
    return "Shape::Elliptical";
  return "Shape::Circular";
}

std::string grmap(enum Gradation gr) {
  if (gr == Gradation::LimbDark)
    return "Gradation::LimbDark";
  return "Gradation::Uniform";
}

//
// class Source
//
// template <typename Type> Source<Type>::Source(Circle<Type> cl)
// {
// 	c_shape = Circular;
// 	c_grad = Uniform;
// 	c_circ = cl;
// 	c_br = 1.0;
// 	c_sfb = c_br / (M_PI * c_circ.r() * c_circ.r());
// }

// template <typename Type> Source<Type>::Source(Circle<Type> cl, Type br) {
//   c_shape = Circular;
//   c_grad = Uniform;
//   c_circ = cl;
//   c_br = br;
//   c_sfb = c_br / (M_PI * c_circ.r() * c_circ.r());
// }

// template <typename Type>
// Source<Type>::Source(Circle<Type> cl, Type br, Type ld) {
//   c_shape = Circular;
//   c_grad = LimbDark;
//   c_circ = cl;
//   c_br = br;
//   c_ld = ld;
//   c_sfb = c_br / (M_PI * c_circ.r() * c_circ.r());
//   c_sb = [this](V2d<Type> p) {
//     return c_br / (1.0 - c_ld / 3.0) *
//            (1.0 - c_ld * (1.0 - std::sqrt(1.0 - p.n2())));
//   };
// }

// template <typename Type>
// Source<Type>::Source(Circle<Type> cl, std::function<Type(V2d<Type>)> f) {
//   c_shape = Circular;
//   c_grad = LimbDark;
//   c_circ = cl;
//   c_sfb = c_br / c_circ.area();
//   c_sb = f;
// }
//
// template <typename Type> Source<Type>::Source(Ellipse<Type> ell) {
//   c_shape = Elliptical;
//   c_grad = Uniform;
//   c_ell = ell;
//   c_br = 1.0;
//   c_sfb = c_br / c_ell.area();
// }
//
// template <typename Type> Source<Type>::Source(Ellipse<Type> ell, Type br) {
//   c_shape = Elliptical;
//   c_grad = Uniform;
//   c_ell = ell;
//   c_br = br;
//   c_sfb = c_br / c_ell.area();
// }
//
// template <typename Type>
// Source<Type>::Source(Ellipse<Type> ell, Type br, Type ld) {
//   c_shape = Elliptical;
//   c_grad = LimbDark;
//   c_ell = ell;
//   c_ld = ld;
//   c_br = br;
//   c_sfb = c_br / c_ell.area();
//   c_sb = [this](V2d<Type> p) {
//     return 1.0 / ((1.0 - 1.0 / c_ld) *
//                   (1.0 - c_ld * (1.0 - std::sqrt(1.0 - c_ell.dist(p) *
//                                                            c_ell.dist(p)))));
//   };
// }
//
// template <typename Type>
// Source<Type>::Source(Ellipse<Type> ell, std::function<Type(V2d<Type>)> f) {
//   c_shape = Elliptical;
//   c_grad = LimbDark;
//   c_ell = ell;
//   c_sfb = c_br / (M_PI / (c_ell.ab().x() * c_ell.ab().y()));
//   c_sb = f;
// }

template <typename Type> Shape Source<Type>::shape() const { return c_shape; }

template <typename Type> Gradation Source<Type>::gradation() const {
  return c_grad;
}

template <typename Type> Circle<Type> Source<Type>::circ() const {
  return c_circ;
}

template <typename Type> Ellipse<Type> Source<Type>::ell() const {
  return c_ell;
}

template <typename Type> V2d<Type> Source<Type>::cent() const {
  switch (c_shape) {
  case Circular:
    return c_circ.cent();
    break;
  case Elliptical:
    return c_ell.cent();
    break;
  default:
    return c_circ.cent();
    break;
  }
  return V2d_NaN<Type>;
}

template <typename Type> Type Source<Type>::rho(V2d<Type> p) {
  switch (c_shape) {
  case Circular:
    return (p - c_circ.cent()).norm() / c_circ.r();
    break;
  case Elliptical:
    return std::numeric_limits<Type>::quiet_NaN();
    break;
  default:
    return std::numeric_limits<Type>::quiet_NaN();
    break;
  }
  return Numeric_NaN<Type>;
}

template <typename Type> Type Source<Type>::sb(V2d<Type> p) const {
  if (c_shape == Shape::Circular) {
    if (c_circ.in(p)) {
      if (c_grad == Gradation::Uniform)
        return c_sfb;
      else
        return c_sb((p - this->cent()) / this->circ().r());
    }
    return 0.0;
  } else {
    if (c_ell.in(p)) {
      if (c_grad == Gradation::Uniform)
        return c_sfb;
      else
        return c_sb(p);
    }
    return 0.0;
  }
}

template <typename Type> Type Source<Type>::area() const {
  if (c_shape == Shape::Circular)
    return c_circ.area();
  return c_ell.area();
}

template <typename Type> Type Source<Type>::brightness() const { return c_br; }

template <typename Type> Source<float> Source<Type>::to_f() const {
  switch (c_shape) {
  case Circular:
    switch (c_grad) {
    case Uniform:
      return Source<float>(c_circ.to_f(), static_cast<float>(c_br));
      break;
    case LimbDark:
      if (c_sb == nullptr)
        return Source<float>(c_circ.to_f(), static_cast<float>(c_br),
                             static_cast<float>(c_ld));
      else
        return Source<float>(c_circ.to_f(), [this](V2d<float> p) {
          return static_cast<float>(c_sb(conv<Type, float>(p)));
        });
      break;
    default:
      //                    return Numeric_NaN<float>;
      break;
    }
    break;
  case Elliptical:
    switch (c_grad) {
    case Uniform:
      return Source<float>(c_ell.to_f(), c_br);
      break;
    case LimbDark:
      if (c_sb == nullptr)
        return Source<float>(c_ell.to_f(), static_cast<float>(c_br),
                             static_cast<float>(c_ld));
      break;
    default:
      //                    return std::numeric_limits<Type>::quiet_NaN();
      break;
    }
    break;
  default:
    switch (c_grad) {
    case Uniform:
      //                    return std::numeric_limits<Type>::quiet_NaN();
      break;
    case LimbDark:
      //                    return std::numeric_limits<Type>::quiet_NaN();
      break;
    default:
      //                    return std::numeric_limits<Type>::quiet_NaN();
      break;
    }
    //            return std::numeric_limits<Type>::quiet_NaN();
    break;
  }
  return Source_NaN<float>;
}

template <typename Type> Source<double> Source<Type>::to_d() const {
  switch (c_shape) {
  case Circular:
    switch (c_grad) {
    case Uniform:
      return Source<double>(c_circ.to_d(), static_cast<double>(c_br));
      break;
    case LimbDark:
      if (c_sb == nullptr)
        return Source<double>(c_circ.to_d(), static_cast<double>(c_br),
                              static_cast<double>(c_ld));
      else
        return Source<double>(c_circ.to_d(), [this](V2d<double> p) {
          return static_cast<double>(c_sb(conv<Type, double>(p)));
        });
      break;
    default:
      //                    return std::numeric_limits<Type>::quiet_NaN();
      break;
    }
    break;
  case Elliptical:
    switch (c_grad) {
    case Uniform:
      return Source<double>(c_ell.to_d(), static_cast<double>(c_br));
      break;
    case LimbDark:
      if (c_sb == nullptr)
        return Source<double>(c_ell.to_d(), static_cast<double>(c_br),
                              static_cast<double>(c_ld));
      break;
    default:
      //                    return std::numeric_limits<Type>::quiet_NaN();
      break;
    }
    break;
  default:
    switch (c_grad) {
    case Uniform:
      //                    return std::numeric_limits<Type>::quiet_NaN();
      break;
    case LimbDark:
      //                    return std::numeric_limits<Type>::quiet_NaN();
      break;
    default:
      //                    return std::numeric_limits<Type>::quiet_NaN();
      break;
    }
    //            return std::numeric_limits<Type>::quiet_NaN();
    break;
  }
  return Source_NaN<double>;
}

template <typename Type> Source<long double> Source<Type>::to_l() const {
  switch (c_shape) {
  case Circular:
    switch (c_grad) {
    case Uniform:
      return Source<long double>(c_circ.to_l(), static_cast<long double>(c_br));
      break;
    case LimbDark:
      if (c_sb == nullptr)
        return Source<long double>(c_circ.to_l(),
                                   static_cast<long double>(c_br),
                                   static_cast<long double>(c_ld));
      else
        return Source<long double>(c_circ.to_l(), [this](V2d<long double> p) {
          return static_cast<long double>(c_sb(conv<Type, long double>(p)));
        });
      break;
    default:
      //                    return std::numeric_limits<Type>::quiet_NaN();
      break;
    }
    break;
  case Elliptical:
    switch (c_grad) {
    case Uniform:
      return Source<long double>(c_ell.to_l(), c_br);
      break;
    case LimbDark:
      if (c_sb == nullptr)
        return Source<long double>(c_ell.to_l(), static_cast<long double>(c_br),
                                   static_cast<long double>(c_ld));
      break;
    default:
      //                    return std::numeric_limits<long
      //                    double>::quiet_NaN();
      break;
    }
    break;
  default:
    switch (c_grad) {
    case Uniform:
      //                    return std::numeric_limits<long
      //                    double>::quiet_NaN();
      break;
    case LimbDark:
      //                    return std::numeric_limits<long
      //                    double>::quiet_NaN();
      break;
    default:
      //                    return std::numeric_limits<long
      //                    double>::quiet_NaN();
      break;
    }
    //            return std::numeric_limits<long double>::quiet_NaN();
    break;
  }
  return Source_NaN<long double>;
}

template <typename Type> void Source<Type>::print() const {
  std::cout << "Shape : " << shmap(c_shape) << std::endl;
  if (c_shape == Circular)
    c_circ.print();
  else
    c_ell.print();
  std::cout << "c_br (total brightness) = " << c_br << std::endl;
  std::cout << "c_sfb (Surface brightness) = " << c_sfb << std::endl;
  std::cout << "Gradation = " << grmap(c_grad) << std::endl;
  if (c_grad == LimbDark)
    std::cout << "c_ld (Limb darkening coefficient) = " << c_ld << std::endl;
}

template <typename Type> void Source<Type>::move_to(V2d<Type> p) {
  switch (c_shape) {
  case Circular:
    c_circ = Circle<Type>(p, c_circ.r());
    break;
  case Elliptical:
    c_ell = Ellipse<Type>(p, c_ell.ab().x(), c_ell.ab().y());
    break;
  default:
    break;
  }
}

//
// Functions
//

//
// Explicit instantation of templates of classes
//
// For float
//
template class Source<float>;
// template class Msource<float>;

//
// For double
//
template class Source<double>;
// template class Msource<double>;

//
// For long double
//
template class Source<long double>;
// template class Msource<long double>;
