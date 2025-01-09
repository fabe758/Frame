/**
 * @file Lens.h
 *
 * @brief header file for functions, classes to handle lens system
 * @author Fumio Abe
 *
 */

// define _INC_Lens_h_ to load only once
#ifndef _INC_Lens_h_
#define _INC_Lens_h_
// include Geom.h
#include "Geom.h"

/**
 * template class Lens
 * @brief Class for a lens
 *
 * Type : float,  double,  long double
 *
 */
template <class Type> class Lens {
private:
  /// @brief Mass ratio
  Type c_q;
  /// @brief Center of the lens
  V2d<Type> c_p;

public:
  /// @brief Blank constructor (all zero)
  Lens() : c_q(1.0), c_p(V2d<Type>(0.0, 0.0)) {};
  /// @brief Constructor specifying mass ratio q at (0.0, 0.0)
  /// @param[in] q : mass ratio of the lens
  Lens(Type q) : c_q(q), c_p(V2d<Type>(0.0, 0.0)) {};
  /// @brief Constructor specifying mass ratio q and the center p
  /// @param[in] q : mass ratio of the lens
  /// @param[in] p : center of the lens
  Lens(Type q, V2d<Type> p) : c_q(q), c_p(p) {};
  /// @brief Return the mass ratio
  /// @return Type q() : the mass ratio
  Type q() const;
  /// @brief Return the Einstein radious
  /// @return Type re() : the Einstein radious
  Type re() const;
  /// @brief Return the Einstein ring
  /// @return Circle<Type> circle() : the Einstein ring
  Circle<Type> circ() const;
  /// @brief Return the center of the lens
  /// @return V2d<Type> p() : the center of the lens
  V2d<Type> p() const;
  /// @brief Return the upper right corner of the fre times large minimum suare
  /// @return V2d<Type> min() : the lower left corner
  /// @param[in] fre : a factor to multiply
  V2d<Type> max(Type fre) const;
  /// @brief Return the lower left corner of the fre times large minimum suare
  /// @return V2d<Type> min() : the lower left corner
  /// @param[in] fre : a factor to multiply
  V2d<Type> min(Type fre) const;
  /// @brief Return the fre times large minimum suare
  /// @param[in] fre : a factor to multiply
  /// @return Rectangle<Type> rec() : the fre times large minimum suare
  Rectangle<Type> rec(Type fre) const;
  /// @brief Return the fre times large minimum suare
  /// @param[in] fre : a factor to multiply
  /// @return Square<Type> sq() : the fre times large minimum suare
  Square<Type> sq(Type fre) const;
  /// @brief Return the function g for simple calculation
  /// @return Type g() : the function to use calculation of derivative
  /// @param[in] th : a position on the lens plane
  /// @param[in] alpha : the power
  Type g(V2d<Type> th, int alpha) const;
  /// @brief Return the converted Lens to float
  /// @return Lens<float> to_f() : the converted Lens to Lens<float> (32-bit
  /// floating point)
  Lens<float> to_f() const;
  /// @brief Return the converted Lens to double
  /// @return Lens<double> to_d() : the converted Lens to Lens<double> (64-bit
  /// floating point)
  Lens<double> to_d() const;
  /// @brief Return the converted Lens to long double
  /// @return Lens<long double> to_l() : the converted Lens to Lens<long double>
  /// (128-bit floating point)
  Lens<long double> to_l() const;
  /// @brief Move the center to the point p
  /// @param[in] p : the position to move the center to
  void move_to(V2d<Type> p);
  /// @brief Print the lens parameters
  void print() const;
  /// @brief Return if the lens is equal to l2
  /// @param[in] l2 : the lens to compare
  /// @return bool operator==() : equal (true) or not equal (false)
  bool operator==(Lens<Type> l2) const;
};

/// @brief Lens with all zero parameters
template <typename Type> const Lens<Type> Lens_0(0.0, V2d_0<Type>);

/**
 * template class Mlens
 * @brief Class for multiple lens system
 *
 * Type : float,  double,  long double
 *
 */
template <class Type> class Mlens : public std::vector<Lens<Type>> {
  using size_type = typename std::vector<Lens<Type>>::size_type;

public:
  /// @brief Use vector of Lens
  using std::vector<Lens<Type>>::vector;
  /// @brief Return specific Lens
  /// @param[in] i : index of the lens
  /// @return Lens<Type> &operator[] : The Lens
  Lens<Type> &operator[](size_type i);
  /// @brief Return specific Lens
  /// @param[in] i : index of the lens
  /// @return Lens<Type> &operator[] : The Lens
  const Lens<Type> &operator[](size_type i) const;
  /// @brief Blank constructor
  // Mlens() : std::vector<Lens<Type>>({Lens<Type>()}) {};
  /// @brief Constructor to initialize by a lens
  /// @param[in] lens : a lens
  Mlens(Lens<Type> lens) : std::vector<Lens<Type>>({lens}) {};
  /// @brief Constructor to initialize only mass ratios
  /// @param[in] vq : mass ratios
  Mlens(std::vector<Type> vq) {
    for (Lens<Type> q : vq)
      this->push_back(Lens(q));
  };
  /// @brief Constructor to initialize by ml
  /// @param[in] ml : multiple lenses
  Mlens(std::vector<Lens<Type>> ml) : std::vector<Lens<Type>>(ml) {};
  /// @brief Sort by the increasing order of q (mass ratio)
  void sort_i();
  /// @brief Sort by the increasing order of q (mass ratio)
  void sort_d();
  /// @brief Change origin of the coordinates to org
  /// @param[in] org : origin
  void origin(V2d<Type> org);
  /// @brief Return the center of mass of the lenses
  /// @return V2d<Type> cm() : center of mass
  V2d<Type> cm() const;
  /// @brief Return summation of the mass ratios
  /// @return Type sumq() : summation of the q (mass ratio)
  Type sumq() const;
  /// @brief Make summation of the mass ratio to be 1.0
  void normq();
  /// @brief Return the maximum mass ratio lens
  /// @return Lens<Type> maxq() : Largest mass ratio lens
  Lens<Type> maxq() const;
  /// @brief Return the minimum mass ratio lens
  /// @return Lens<Type> minq() : Smallest mass ratio lens
  Lens<Type> minq() const;
  /// @brief Return the closest lens from the point p
  /// @param[in] p : a point
  /// @return Lens<Type> closest() : the close3st lens from the point p
  Lens<Type> closest(V2d<Type> p) const;
  /// @brief Return upper right corner of the fre times large lens square
  /// @param[in] fre : a factor to enlarge
  /// @return V2d<Type> max(Type fre) : Upper right corner of the lens square
  V2d<Type> max(Type fre) const;
  /// @brief Return lower left corner of the fre times large lens square
  /// @param[in] fre : a factor to enlarge
  /// @return V2d<Type> min(Type fre) : Lower left corner of the lens square
  V2d<Type> min(Type fre) const;
  /// @brief Return the fre times large lens square
  /// @param[in] fre : a factor to enlarge
  /// @return Square<Type> rec(Type fre) : the fre times large lens square
  Square<Type> sq(Type fre);
  /// @brief Return the Jacobian matrix at th (a position on the lens plane)
  /// @param[in] th : a position on the lens plane
  /// @return M2d<Type> J(V2d<Type> th) : the Jacobian matrix
  M2d<Type> J(V2d<Type> th) const;
  /// @brief Return the magnification at th (a position on the lens plane)
  /// @param[in] th : a position on the lens plane
  /// @return mu(V2d<Type> th) : the magnification at th
  Type mu(V2d<Type> th) const;
  /// @brief Return the mapped position on the source plane by the lens system
  /// @param[in] th : a position on the lens plane
  /// @return V2d<Type> beta(V2d<Type> th) : the mapped position on the source
  /// plane
  V2d<Type> beta(V2d<Type> th) const;
  /// @brief Return the distance between the mapped position and a position on
  /// the source plane
  /// @param[in] th : a position on the lens plane
  /// @param[in] bt : a position on the source plane
  /// @return V2d<Type> db() : the distance between the mapped position and bt
  Type db(V2d<Type> th, V2d<Type> bt) const;
  /// @brief Add a new lens (does not work in PyROOT)
  void operator+=(Lens<Type> l);
  /// @brief Convert to Mlens<float> (32-bit floating point)
  /// @return Mlens<float> to_f() : converted Mlens to float
  Mlens<float> to_f() const;
  /// @brief Convert to Mlens<double> (64-bit floating point)
  /// @return Mlens<double> to_d() : converted Mlens to double
  Mlens<double> to_d() const;
  /// @brief Convert to Mlens<long double> (128-bit floating point)
  /// @return Mlens<long double> to_l() : converted Mlens to long double
  Mlens<long double> to_l() const;
  /// @brief Print parameters
  void print() const;
  /// @brief Move positions of the lenses
  /// @param[in] vp : the positions of the lelses to move to
  void move_to(std::vector<V2d<Type>> vp);
};

/// @brief Global null Mlens constant
template <typename Type> const Mlens<Type> Mlens_0 = Mlens<Type>(Lens_0<Type>);

/// @brief Numerical limit for Lens
template <typename Type>
const Lens<Type> Lens_NaN(Numeric_NaN<Type>, V2d_NaN<Type>);

/// @brief Numerical limit for Mlens
template <typename Type> const Mlens<Type> Mlens_NaN{Lens_NaN<Type>};

//
// Suppress instantation of templates of classes
//
// For float
//
extern template class Lens<float>;
extern template class Mlens<float>;

/// @brief using fLens = Lens<float>;
using fLens = Lens<float>;
/// @brief using fMlens = Mlens<float>;
using fMlens = Mlens<float>;

//
// For double
//
extern template class Lens<double>;
extern template class Mlens<double>;

/// @brief using dLens = Lens<double>;
using dLens = Lens<double>;
/// @brief using dMlens = Mlens<double>;
using dMlens = Mlens<double>;

//
// For long double
//
extern template class Lens<long double>;
extern template class Mlens<long double>;

/// @brief using lLens = Lens<long double>;
using lLens = Lens<long double>;
/// @brief using lMlens = Mlens<long double>;
using lMlens = Mlens<long double>;

//
//  Function prototypes
//

#endif
