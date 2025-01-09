/**
 * @file Source.h
 *
 * @brief header file for Source star
 * @author Fumio Abe
 *
 */

// define _INC_Source_h_ to load only once
#ifndef _INC_Source_h_
#define _INC_Source_h_

// include Lens.h
#include "Lens.h"

/// @brief enum for Shape : Circular or Elliptical
enum Shape { Circular, Elliptical };

/// @brief enum for Gradation : Uniform or LimbDark
enum Gradation { Uniform, LimbDark };

/// @brief Convert Shape to string
std::string shmap(enum Shape sh);

/// @brief Convert Gradation to string
std::string grmap(enum Gradation gr);

/**
 * template class Source
 * @brief class for source star
 *
 * Type : float,  double,  long double
 * Aliases : fSource (Source<float),
 *           dSource (Source<double),
 *           lSource (Source<long double)
 *
 */
template <class Type> class Source {
private:
  /// @brief Shape of the source star
  Shape c_shape;
  /// @brief Gradation of the surface brightness
  Gradation c_grad;
  /// @brief Total brightness
  Type c_br;
  /// @brief Surface brightness (if Gradation::Uniform)
  Type c_sfb;
  /// @brief Limb darkening coefficient
  Type c_ld;
  /// @brief Circular shape of the star (if Shape::Circular)
  Circle<Type> c_circ;
  /// @brief Elliptical shape of the star (if Shape::Elliptical)
  Ellipse<Type> c_ell;
  /// @brief 2d profile of the surface brightness
  /// @param[in] V2d<Type> : 2d position on the star
  /// @return Type : Surface brightness
  std::function<Type(V2d<Type>)> c_sb;

public:
  /// @brief Blank constructor : Circular, Uniform, total 1.0, radius 0.01, at
  /// (0.1, 0.0)
  Source()
      : c_shape(Circular), c_grad(Uniform), c_br(1.0), c_sfb(0.0), c_ld(0.0),
        c_circ(Circle<Type>(V2d<Type>(0.1, 0.0), 0.01)), c_ell(Ellipse_0<Type>),
        c_sb(nullptr) {
    c_sfb = c_br / (M_PI * c_circ.r() * c_circ.r());
  };
  Source(Type rho)
      : c_shape(Circular), c_grad(Uniform), c_br(1.0),
        c_sfb(c_br / (M_PI * rho * rho)), c_ld(0.0),
        c_circ(Circle<Type>(V2d<Type>(0.0, 0.0), rho)), c_ell(Ellipse_0<Type>),
        c_sb(nullptr) {
    c_sfb = c_br / (M_PI * c_circ.r() * c_circ.r());
  };
  /// @brief Constructor : Circular, Uniform, total 1.0, Circle<Type> cl
  /// @param[in] cl : circular shape of the star
  Source(Circle<Type> cl)
      : c_shape(Circular), c_grad(Uniform), c_br(1.0),
        c_sfb(c_br / (M_PI * cl.r() * cl.r())), c_ld(0.0), c_circ(cl),
        c_ell(Ellipse_0<Type>), c_sb(nullptr) {
    c_sfb = c_br / (M_PI * c_circ.r() * c_circ.r());
  };
  /// @brief Constructor : Circular, Uniform, total br, Circle<Type> cl
  /// @param[in] cl : circular shape of the star
  /// @param[in] br : total brightness
  Source(Circle<Type> cl, Type br)
      : c_shape(Circular), c_grad(Uniform), c_br(br),
        c_sfb(br / (M_PI * cl.r() * cl.r())), c_circ(cl),
        c_ell(Ellipse_0<Type>), c_sb(nullptr) {
    c_sfb = c_br / (M_PI * c_circ.r() * c_circ.r());
  };
  /// @brief Constructor : Circular, LimbDark, total br, Circle<Type> cl
  /// @param[in] cl : circular shape of the star
  /// @param[in] br : total brightness
  /// @param[in] ld : limb darkening coefficient (Bozza 2010)
  Source(Circle<Type> cl, Type br, Type ld)
      : c_shape(Circular), c_grad(LimbDark), c_br(br),
        c_sfb(br / (M_PI * cl.r() * cl.r())), c_ld(ld), c_circ(cl),
        c_ell(Ellipse_0<Type>), c_sb([&](V2d<Type> p) {
          return 1.0 / (1.0 - c_ld / 3.0) *
                 (1.0 - c_ld * (1.0 - std::sqrt(1.0 - p * p)));
        }) {
    c_sfb = c_br / (M_PI * c_circ.r() * c_circ.r());
  };
  /// @brief Constructor : Circular, non Uniform, Circle<Type> cl
  /// @param[in] cl : circular shape of the star
  /// @param[in] f : 2d function that represent the surface brightness
  Source(Circle<Type> cl, std::function<Type(V2d<Type>)> f)
      : c_shape(Circular), c_grad(LimbDark), c_br(1.0),
        c_sfb(c_br / (M_PI * cl.r() * cl.r())), c_circ(cl),
        c_ell(Ellipse_0<Type>), c_sb(f) {
    c_sfb = c_br / (M_PI * c_circ.r() * c_circ.r());
  };
  /// @brief Constructor : Elliptical, Uniform, total 1.0, Ellipse<Type> ell
  /// @param[in] ell : Elliptical shape of the star
  Source(Ellipse<Type> ell)
      : c_shape(Elliptical), c_grad(Uniform), c_br(1.0),
        c_sfb(c_br / (M_PI * ell.ab().x() * ell.ab().y())), c_ld(0.0),
        c_circ(Circle_0<Type>), c_ell(ell), c_sb(nullptr) {
    c_sfb = c_br / (M_PI * c_circ.r() * c_circ.r());
  };
  /// @brief Constructor : Elliptical, Uniform, total br, Ellipse<Type> ell
  /// @param[in] ell : Elliptical shape of the star
  /// @param[in] br : total brightness
  Source(Ellipse<Type> ell, Type br)
      : c_shape(Elliptical), c_grad(Uniform), c_br(br),
        c_sfb(c_br / (M_PI * ell.ab().x() * ell.ab().y())), c_ld(0.0),
        c_circ(Circle_0<Type>), c_ell(ell), c_sb(nullptr) {
    c_sfb = c_br / (M_PI * c_circ.r() * c_circ.r());
  };
  /// @brief Constructor : Elliptical, LimbDark, total br, Ellipse<Type> ell
  /// @param[in] ell : Elliptical shape of the star
  /// @param[in] br : total brightness
  /// @param[in] ld : limb darkening coefficient (Bozza 2010)
  Source(Ellipse<Type> ell, Type br, Type ld)
      : c_shape(Elliptical), c_grad(LimbDark), c_br(br),
        c_sfb(c_br / (M_PI * ell.ab().x() * ell.ab().y())), c_ld(ld),
        c_circ(Circle_0<Type>), c_ell(ell), c_sb([&](V2d<Type> p) {
          return 1.0 / (1.0 - c_ld / 3.0) *
                 (1.0 - c_ld * (1.0 - std::sqrt(1.0 - p * p)));
        }) {};
  /// @brief Constructor : Elliptical, non Uniform, Ellipse<Type> ell
  /// @param[in] ell : Elliptical shape of the star
  /// @param[in] f : 2d function that represent the surface brightness
  Source(Ellipse<Type> ell, std::function<Type(V2d<Type>)> f)
      : c_shape(Elliptical), c_grad(LimbDark), c_br(1.0),
        c_sfb(c_br / (M_PI * ell.ab().x() * ell.ab().y())),
        c_circ(Circle_0<Type>), c_ell(ell), c_sb(f) {
    c_sfb = c_br / (M_PI * c_circ.r() * c_circ.r());
  };
  /// @brief Return the Shape (Circular or Elliptical)
  /// @return Shape shape() :  Shape of the star (Circular or Elliptical)
  Shape shape() const;
  /// @brief Return the Gradation (Uniform or LimbDark)
  /// @return Gradation gradation() : Gradation (Uniform or LimbDark)
  Gradation gradation() const;
  /// @brief Return the circular shape
  /// @return Circle<Type> circ() :  circular shape
  Circle<Type> circ() const;
  /// @brief Return the elliptical shape
  /// @return Ellipse<Type> ell() :  elliptical shape
  Ellipse<Type> ell() const;
  /// @brief Return the center of the star
  /// @return V2d<Type> cent() : center of the star
  V2d<Type> cent() const;
  /// @brief Return the circular or elliptical normalized distance from the
  /// center
  /// @param[in] p : 2d position in the star
  /// @return Type rho(V2d<Type> p) : Circular or Elliptical "normalized"
  /// distance from the center
  Type rho(V2d<Type> p);
  /// @brief Return the surface brightness of the star at p
  /// @param[in] p : 2d position in the star
  /// @return Type sb(V2d<Type> p) :  surface brightness
  Type sb(V2d<Type> p) const;
  /// @brief Return the surface brightness of the star at p
  /// @return Type sb(V2d<Type> p) :  surface brightness
  Type area() const;
  /// @brief Return Limb Darkening coefficient
  /// @return limbdark() : Return Limb Darkening coefficient
  Type limbdark() const { return c_ld; };
  /// @brief Return the total brightness of the star
  /// @return Type brightness() :  total brightness
  Type brightness() const;
  /// @brief Returns Source<float> (32-bit floating point)
  /// @return Source<float> to_f() : converted Source to float
  Source<float> to_f() const;
  /// @brief Returns Source<double> (64-bit floating point)
  /// @return Source<double> to_d() : converted Source to double
  Source<double> to_d() const;
  /// @brief Returns Source<long double> (128-bit floating point)
  /// @return Source<long double> to_l() : converted Source to long double
  Source<long double> to_l() const;
  /// @brief Print parameters
  void print() const;
  /// @brief Move the center to p
  /// @param[in] p : the point to move to
  void move_to(V2d<Type> p);
};

/// brief Global constant Source (Circlular, Uniform, all 0.0)
template <typename Type>
const Source<Type> Source_0 = Source<Type>(Circle_0<Type>);

/// brief Global constant Source exception
template <typename Type> const Source<Type> Source_NaN(Circle_NaN<Type>);

//
//  Function prototypes
//

//
// Suppress instantation of templates of classes
//
// For float
//
extern template class Source<float>;

/// @brief  using fSource = Source<float>;
using fSource = Source<float>;

//
// For double
//
extern template class Source<double>;

/// @brief  using dSource = Source<double>;
using dSource = Source<double>;

//
// For long double
//
extern template class Source<long double>;

/// @brief  using lSource = Source<long double>;
using lSource = Source<long double>;

#endif
