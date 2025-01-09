/**
 * @file Geom.cpp
 +
 * @brief Class library to handle geometries
 * @author Fumio Abe
 * @details tow dimensional vector,  matrix, rectangle, Circle,
 * Ellipse,  etc.
 *
 */

/**
 * include headers  for this project
 * all c++ headers are loaded in Geom.h
 *
 */
#include "Geom.h" 

/**
 * method V2d<Type>::x()
 * @brief returns x coordinate
 * @return Type x() : x coordinate
 *
 */
template <typename Type> Type V2d<Type>::x() const { return c_x; }

/**
 * method V2d<Type>::y()
 * @brief returns y coordinate
 * @return Type y() : y coordinate
 *
 */
template <typename Type> Type V2d<Type>::y() const { return c_y; }

/**
 * method V2d<Type>::norm()
 * @brief returns norm of the vector
 * @return Type norm() : norm of the vector
 *
 */
template <typename Type> Type V2d<Type>::norm() const {
  return std::sqrt(c_x * c_x + c_y * c_y);
}

/**
 * method V2d<Type>::norm()
 * @brief returns norm of the vector
 * @return Type norm() : square norm of the vector
 *
 */
template <typename Type> Type V2d<Type>::n2() const {
  return c_x * c_x + c_y * c_y;
}

/**
 * method V2d<Type>::unit()
 * @brief returns unit vector with same direction
 * @return Type unit() : unit vector with same direction
 *
 */
template <typename Type> V2d<Type> V2d<Type>::unit() const {
  return *this / this->norm();
}

/**
 * method V2d<Type>::online()
 * @brief returns if the point is on the line of v1 and v2 or not
 * @return Type online() :the point is on the line of v1 and v2
 * @param[in] v1 : a 2d vector v1
 * @param[in] v2 : a 2d vector v2
 *
 */
template <typename Type> bool V2d<Type>::online(V2d<Type> v1, V2d<Type> v2) {
  bool b1 = v1 == *this;
  bool b2 = v2 == *this;
  if (b1 || b2)
    return true;
  return this->between(v1, v2) ||
         ((v1 - *this).unit() * (v2 - *this).unit()) - 1.0 <=
             4.0*std::numeric_limits<Type>::epsilon();
}

/**
 * method V2d<Type>::between()
 * @brief returns if the point is between v1 and v2 or not
 * @return Type online() :the point is between v1 and v2
 * @param[in] v1 : a 2d vector v1
 * @param[in] v2 : a 2d vector v2
 *
 */
template <typename Type> bool V2d<Type>::between(V2d<Type> v1, V2d<Type> v2) {
  bool b1 = v1 == *this;
  bool b2 = v2 == *this;
  if (b1 || b2)
    return true;
  return ((v1 - *this).unit() * (v2 - *this).unit()) + 1.0 <=
         4.0*std::numeric_limits<Type>::epsilon();
}

/**
 * operator V2d<Type>::operator+(V2d<Type> vv)
 * @brief returns the vector added by vv
 * @return V2d<Type> : vector added by vv
 * @param[in] vv : a 2d vector vv
 *
 */
template <typename Type> V2d<Type> V2d<Type>::operator+(V2d<Type> vv) const {
  return V2d(c_x + vv.x(), c_y + vv.y());
}

/**
 * operator V2d<Type>::operator-(V2d<Type> vv)
 * @brief returns the vector subtracted by vv
 * @return V2d<Type> : vector subtracted by vv
 * @param[in] vv : a 2d vector vv
 *
 */
template <typename Type> V2d<Type> V2d<Type>::operator-(V2d<Type> vv) const {
  return V2d(c_x - vv.x(), c_y - vv.y());
}

/**
 * operator V2d<Type>::operator*(V2d<Type> vv)
 * @brief returns inner product of the vector and vv
 * @return Type : inner product with vv
 * @param[in] vv : a 2d vector vv
 *
 */
template <typename Type> Type V2d<Type>::operator*(V2d<Type> vv) const {
  return c_x * vv.x() + c_y * vv.y();
}

/**
 * operator V2d<Type>::operator*(Type dd)
 * @brief returns the vector multiplied by a number dd
 * @param[in] dd : a 2d vector vv
 * @return V2d<Type> : vector multiplied by dd
 *
 */
template <typename Type> V2d<Type> V2d<Type>::operator*(Type dd) const {
  return V2d(c_x * dd, c_y * dd);
}

/**
 * operator V2d<Type>::operator*(M2d<Type> mm)
 * @brief returns the vector multiplied by a matrix mm
 * @param[in] mm : a 2d matrix mm
 * @return V2d<Type> : vector multiplied by mm
 *
 */
template <typename Type> V2d<Type> V2d<Type>::operator*(M2d<Type> mm) const {
  return V2d<Type>(c_x * mm.xx() + c_y * mm.yx(),
                   c_x * mm.xy() + c_y * mm.yy());
}

/**
 * operator V2d<Type>::operator/(Type dd)
 * @brief returns the vector devided by a number dd
 * @param[in] dd : a number dd
 * @return V2d<Type> : vector devided by dd
 *
 */
template <typename Type> V2d<Type> V2d<Type>::operator/(Type dd) const {
  return V2d(c_x / dd, c_y / dd);
}

/**
 * operator V2d<Type>::operator+=(V2d<Type> vv)
 * @brief adds a vector vv
 *
 */
template <typename Type> void V2d<Type>::operator+=(V2d<Type> vv) {
  c_x += vv.x();
  c_y += vv.y();
}

/**
 * operator V2d<Type>::operator-=(V2d<Type> vv)
 * @brief subtracts a vector vv
 *
 */
template <typename Type> void V2d<Type>::operator-=(V2d<Type> vv) {
  c_x -= vv.x();
  c_y -= vv.y();
}

/**
 * operator V2d<Type>::operator==(V2d<Type> vv)
 * @brief check if the vector is equal to a vector vv or not
 * @return bool : true if equal,  faulse if not
 *
 */
template <typename Type> bool V2d<Type>::operator==(V2d<Type> vv) const {
  return std::abs(vv.x() - c_x) <= 4.0*std::numeric_limits<Type>::epsilon() &&
         std::abs(vv.y() - c_y) <= 4.0*std::numeric_limits<Type>::epsilon();
}

/**
 * operator V2d<Type>::operator^(Type ph)
 * @brief returns the vector rotated by an angle ph (in radian)
 * @return V2d<Type> : vector rotated by ph
 *
 */
template <typename Type> V2d<Type> V2d<Type>::operator^(Type ph) const {
  M2d<Type> m(std::cos(ph), -std::sin(ph), std::sin(ph), std::cos(ph));
  return m * *this;
}

/**
 * operator V2d<Type>::operator<(V2d<Type> vv)
 * @brief returns the "vector product" with a vector vv (only z component)
 * @return Type : vector product with vv
 *
 */
template <typename Type> Type V2d<Type>::operator%(V2d<Type> vv) const {
  return c_x * vv.y() - c_y * vv.x();
}

/**
 * method V2d<Type>::to_f()
 * @brief returns converted vector to V2d<float> (32-bit floating point)
 * @return V2d<float> : converted vector to V2d<float>
 *
 */
template <typename Type> V2d<float> V2d<Type>::to_f() const {
  return V2d<float>(static_cast<float>(c_x), static_cast<float>(c_y));
}

/**
 * method V2d<Type>::to_d()
 * @brief returns converted vector to V2d<double> (64-bit floating point)
 * @return V2d<double> : converted vector to V2d<double>
 *
 */
template <typename Type> V2d<double> V2d<Type>::to_d() const {
  return V2d<double>(static_cast<double>(c_x), static_cast<double>(c_y));
}

/**
 * method V2d<Type>::to_l()
 * @brief returns converted vector to V2d<long double> (128-bit floating point)
 * @return V2d<long double> : converted vector to V2d<long double>
 *
 */
template <typename Type> V2d<long double> V2d<Type>::to_l() const {
  return V2d<long double>(static_cast<long double>(c_x),
                          static_cast<long double>(c_y));
}

/**
 * method V2d<Type>::print()
 * @brief prints  x and y components
 *
 */
template <typename Type> void V2d<Type>::print() const {
  std::cout << "x = " << c_x << ",   y = " << c_y << std::endl;
}

/**
 * template class Seg2d
 * @brief Type : float,  double, long double
 * @brief Aliases : fSeg2d (Seg2d<float>),dSeg2d (Seg2d<double>),lSeg2d (Seg2d<long double>),
 * @brief Class to express 2d line Segment
 */
/**
 * method std::array<V2d<Type>, 2> Seg2d<Type>::ends()
 * @brief returns end points of the segment
 * @return method std::array<V2d<Type>, 2> ends() : points of the segment
 *
 */
template <class Type> std::array<V2d<Type>, 2> Seg2d<Type>::ends() const {
  return c_ends;
}

/**
 * method bool Seg2d<Type>::on(V2d<Type> p)
 * @brief returns if point p is on the segment or not
 * @return bool on() : the point is on the segment or not (true or false) 
 * @param[in] p : a point on the plane
 *
 */
template <class Type> bool Seg2d<Type>::on(V2d<Type> p) const {
  return p.between(c_ends[0], c_ends[1]);
}

/**
 * method Type Seg2d<Type>::length()
 * @brief returns the length of the segment
 * @return Type Seg2d<Type>::length() : the length of the segment
 *
 */
template <class Type> Type Seg2d<Type>::length() const {
  return (c_ends[1] - c_ends[0]).norm();
}

/**
 * method Type Seg2d<Type>::alpha(V2d<Type> p)
 * @brief returns the length between end[0] and foot point of point p
 * @return Type Seg2d<Type>::alpha(V2d<Type> p) : the length between the end[0] and foot point of point p
 * @param[in] p : a point on the plane
 *
 */
template <class Type> Type Seg2d<Type>::alpha(V2d<Type> p) const {
  return -1.0 * ((c_ends[0] - p) * (c_ends[0] - c_ends[1]).unit());
}

/**
 * method V2d<Type> Seg2d<Type>::foot(V2d<Type> p)
 * @brief returns the foot point of the point p
 * @return V2d<Type> Seg2d<Type>::foot() : the foot point
 * @param[in] p : a point on the plane
 *
 */
template <class Type> V2d<Type> Seg2d<Type>::foot(V2d<Type> p) const {
  return c_ends[0] + (c_ends[0] - c_ends[1]).unit() * this->alpha(p);
}

/**
 * method V2d<Type> Seg2d<Type>::closest(V2d<Type> p)
 * @brief returns the closest point on the segment to the point p
 * @return V2d<Type> Seg2d<Type>::foot() : the foot point
 * @param[in] p : a point on the plane
 *
 */
template <class Type> V2d<Type> Seg2d<Type>::closest(V2d<Type> p) const {
  if (!this->on(this->foot(p)))
    return *std::min_element(
        c_ends.begin(), c_ends.end(),
        [&p](auto a, auto b) { return (a - p).norm() < (b - p).norm(); });
  return this->foot(p);
}

/**
 * method Type Seg2d<Type>::distance(V2d<Type> p)
 * @brief returns the closest point on the segment to the point p
 * @return Type Seg2d<Type>::distance() : the distance of the point p
 * @param[in] p : a point on the plane
 *
 */
template <class Type> Type Seg2d<Type>::distance(V2d<Type> p) const {
  return (this->closest(p) - p).norm();
}

/**
 * method Type Seg2d<Type>::angle(V2d<Type> p)
 * @brief returns the opening angle of the segment seen from the point p
 * @return Type Seg2d<Type>::angle() : the opening angle
 * @param[in] p : a point on the plane
 *
 */
template <typename Type> Type Seg2d<Type>::angle(V2d<Type> p) const {
  return std::abs(std::asin((c_ends[0] - p).unit() % (c_ends[1] - p).unit()));
}

/**
 * method Seg2d<float> Seg2d<Type>::to_f()
 * @brief returns converted segment to Seg2d<float> (32-bit floating point)
 * @return Seg2d<float> : converted segment to Seg2d<float>
 *
 */
template <class Type> Seg2d<float> Seg2d<Type>::to_f() const {
  return Seg2d<float>({c_ends[0].to_f(), c_ends[1].to_f()});
}

/**
 * method Seg2d<double> Seg2d<Type>::to_d()
 * @brief returns converted segment to Seg2d<double> (64-bit floating point)
 * @return Seg2d<double> : converted segment to Seg2d<double>
 *
 */
template <class Type> Seg2d<double> Seg2d<Type>::to_d() const {
  return Seg2d<double>({c_ends[0].to_d(), c_ends[1].to_d()});
}

/**
 * method Seg2d<long double> Seg2d<Type>::to_l()
 * @brief returns converted segment to Seg2d<long double> (128-bit floating point)
 * @return Seg2d<long double> : converted segment to Seg2d<long double>
 *
 */
template <class Type> Seg2d<long double> Seg2d<Type>::to_l() const {
  return Seg2d<long double>({c_ends[0].to_l(), c_ends[1].to_l()});
}

/**
 * method Seg2d<Type>::print()
 * @brief prints end points
 *
 */
template <class Type> void Seg2d<Type>::print() {
  std::cout << "End point 0 : " << std::endl;
  c_ends[0].print();
  std::cout << "End point 1 : " << std::endl;
  c_ends[1].print();
}


/**
 * method M2d<Type>::xx()
 * @brief returns xx component
 * @return Type : xx component
 *
 */
template <typename Type> Type M2d<Type>::xx() const { return c_xx; }

/**
 * method M2d<Type>::xy()
 * @brief returns xy component
 * @return Type : xy component
 *
 */
template <typename Type> Type M2d<Type>::xy() const { return c_xy; }

/**
 * method M2d<Type>::yx()
 * @brief returns yx component
 * @return Type : yx component
 *
 */
template <typename Type> Type M2d<Type>::yx() const { return c_yx; }

/**
 * method M2d<Type>::yy()
 * @brief returns yy component
 * @return Type : yy component
 *
 */
template <typename Type> Type M2d<Type>::yy() const { return c_yy; }

/**
 * operator M2d<Type>::operator+(M2d<Type> mm)
 * @brief returns the matrix added by mm
 * @return M2d<Type> : matrix added by mm
 *
 */
template <typename Type> M2d<Type> M2d<Type>::operator+(M2d<Type> mm) const {
  return M2d<Type>(c_xx + mm.xx(), c_xy + mm.xy(), c_yx + mm.yx(),
                   c_yy + mm.yy());
}

/**
 * operator M2d<Type>::operator-(M2d<Type> mm)
 * @brief returns the matrix subtracted by mm
 * @return M2d<Type> : matrix subtracted by mm
 *
 */
template <typename Type> M2d<Type> M2d<Type>::operator-(M2d<Type> mm) const {
  return M2d<Type>(c_xx - mm.xx(), c_xy - mm.xy(), c_yx - mm.yx(),
                   c_yy - mm.yy());
}

/**
 * operator M2d<Type>::operator*(M2d<Type> mm)
 * @brief returns the matrix mutiplied by mm
 * @return M2d<Type> : matrix multiplied by mm
 *
 */
template <typename Type> M2d<Type> M2d<Type>::operator*(M2d<Type> mm) const {
  return M2d<Type>(
      c_xx * mm.xx() + c_xy * mm.yx(), c_xx * mm.xy() + c_xy * mm.yy(),
      c_yx * mm.xx() + c_yy * mm.yx(), c_yx * mm.xy() + c_yy * mm.yy());
}

/**
 * operator M2d<Type>::operator*(V2d<Type> vv)
 * @brief returns the mutiplied matrix by a vector vv
 * @return V2d<Type> : multiplied by vv
 *
 */
template <typename Type> V2d<Type> M2d<Type>::operator*(V2d<Type> vv) const {
  return V2d<Type>(c_xx * vv.x() + c_xy * vv.y(),
                   c_yx * vv.x() + c_yy * vv.y());
}

/**
 * operator M2d<Type>::operator*(Type ss)
 * @brief returns the matrix mutiplied by a number ss
 * @return V2d<Type> : the matrix multiplied by ss
 *
 */
template <typename Type> M2d<Type> M2d<Type>::operator*(Type ss) const {
  return M2d<Type>(c_xx * ss, c_xy * ss, c_yx * ss, c_yy * ss);
}

/**
 * operator void M2d<Type>::operator+=(M2d<Type> mm)
 * @brief accumulates by the matrix mm
 *
 */
template <typename Type> void M2d<Type>::operator+=(M2d<Type> mm) {
  c_xx += mm.xx();
  c_xy += mm.xy();
  c_yx += mm.yx();
  c_yy += mm.yy();
}

/**
 * operator M2d<Type>::operator!()
 * @brief returns the transposed matrix 
 * @return M2d<Type> : the transposed matrix
 *
 */
template <typename Type> M2d<Type> M2d<Type>::operator!() const {
  return M2d<Type>(c_xx, c_yx, c_xy, c_yy);
}

/**
 * method Type M2d<Type>::det()
 * @brief returns the determinant of the matrix
 * @return Type M2d<Type>::det() : the determinant of the matrix
 *
 */
template <typename Type> Type M2d<Type>::det() const {
  return c_xx * c_yy - c_xy * c_yx;
}

/**
 * method M2d<Type> M2d<Type>::inv()
 * @brief returns the inverse matrix
 * @return M2d<Type> M2d<Type>::inv() : the inverse matrix
 *
 */
template <typename Type> M2d<Type> M2d<Type>::inv() const {
  if (std::abs(det()) <= std::numeric_limits<Type>::epsilon()) {
    constexpr Type de = std::numeric_limits<Type>::min();
    return M2d<Type>(c_yy / de, -c_xy / de, -c_yx / de, c_xx / de);
  } else
    return M2d(c_yy / det(), -c_xy / det(), -c_yx / det(), c_xx / det());
}

/**
 * method std::array<Type, 2> M2d<Type>::lambda()
 * @brief returns eigen values
 * @return std::array<Type, 2> lambda() : eigenvalues
 *
 */
template <typename Type> std::array<Type, 2> M2d<Type>::lambda() const {
  Type l0 = c_xx + c_yy;
  Type l1 = std::sqrt(pow(l0, 2) + 4.0 * (c_xy * c_yx - c_xy * c_yx));
  return std::array<Type, 2>{l0 + l1, l0 - l1};
}

/**
 * method M2d<float> M2d<Type>::to_f()
 * @brief returns converted matrix to M2d<float> (32-bit floating point)
 * @return M2d<float> to_f() : converted matrix to M2d<float>
 *
 */
template <typename Type> M2d<float> M2d<Type>::to_f() const {
  return M2d<float>(static_cast<float>(c_xx), static_cast<float>(c_xy),
                    static_cast<float>(c_yx), static_cast<float>(c_yy));
}

/**
 * method M2d<double> M2d<Type>::to_d()
 * @brief returns converted matrix to M2d<double> (64-bit floating point)
 * @return M2d<double> to_d() : converted matrix to M2d<double>
 *
 */
template <typename Type> M2d<double> M2d<Type>::to_d() const {
  return M2d<double>(static_cast<double>(c_xx), static_cast<double>(c_xy),
                     static_cast<double>(c_yx), static_cast<double>(c_yy));
}

/**
 * method M2d<long double> M2d<Type>::to_l()
 * @brief returns converted matrix to M2d<long double> (128-bit floating point)
 * @return M2d<double> to_l() : converted matrix to M2d<long double>
 *
 */
template <typename Type> M2d<long double> M2d<Type>::to_l() const {
  return M2d<long double>(
      static_cast<long double>(c_xx), static_cast<long double>(c_xy),
      static_cast<long double>(c_yx), static_cast<long double>(c_yy));
}

/**
 * method M2d<Type>::print()
 * @brief prints matrix
 *
 */
template <typename Type> void M2d<Type>::print() const {
  std::cout << "           x               y  " << std::endl;
  std::cout << "x " << c_xx << "    " << c_xy << std::endl;
  std::cout << "y " << c_yx << "    " << c_yy << std::endl;
}

/**
 * method std::array<V2d<Type>, 3> Triangle<Type>::apex()
 * @brief Returns three apexes
 * @return std::array<V2d<Type>, 3> Triangle<Type>::apex() : three apexes
 *
 */
template <typename Type> std::array<V2d<Type>, 3> Triangle<Type>::apex() const {
  return c_apex;
}

template <typename Type> V2d<Type> Triangle<Type>::min() const {
  return V2d<Type>((*std::min_element(c_apex.begin(), c_apex.end(), 
          [](auto a, auto b){return a.x() < b.x();})).x(),
      (*std::min_element(c_apex.begin(), c_apex.end(), 
                         [](auto a, auto b){return a.y() < b.y();})).y());
}

template <typename Type> V2d<Type> Triangle<Type>::max() const {
  return V2d<Type>((*std::max_element(c_apex.begin(), c_apex.end(), 
          [](auto a, auto b){return a.x() < b.x();})).x(),
      (*std::max_element(c_apex.begin(), c_apex.end(), 
                         [](auto a, auto b){return a.y() < b.y();})).y());
}

/**
 * method std::array<Seg2d<Type>, 3> Triangle<Type>::sides()
 * @brief Returns sides of the triangle
 * @return std::array<Seg2d<Type>, 3> sides() : sides of the triangle
 *
 */
template <typename Type>
std::array<Seg2d<Type>, 3> Triangle<Type>::sides() const {
  return {Seg2d<Type>({c_apex[0], c_apex[1]}),
          Seg2d<Type>({c_apex[1], c_apex[2]}),
          Seg2d<Type>({c_apex[2], c_apex[0]})};
}

/**
 * method V2d<Type> Triangle<Type>::cent()
 * @brief Returns the  center of the triangle
 * @return V2d<Type> cent() : center of the triangle
 *
 */
template <typename Type> V2d<Type> Triangle<Type>::cent() const {
  return (c_apex[0] + c_apex[1] + c_apex[2]) / 3.0;
}

/**
 * method Type Triangle<Type>::area()
 * @brief Returns the area of the triangle
 * @return Type area() : area of the triangle
 *
 */
template <typename Type> Type Triangle<Type>::area() const {
  return 0.5 *
         std::abs(
             (this->sides()[0].ends()[1] - this->sides()[0].ends()[0]).x() *
                 (this->sides()[0].ends()[1] + this->sides()[0].ends()[0]).y() +
             (this->sides()[1].ends()[1] - this->sides()[1].ends()[0]).x() *
                 (this->sides()[1].ends()[1] + this->sides()[1].ends()[0]).y() +
             (this->sides()[2].ends()[1] - this->sides()[2].ends()[0]).x() *
                 (this->sides()[2].ends()[1] + this->sides()[2].ends()[0]).y());
}

/**
 * method Type Triangle<Type>::closestd(V2d<Type> p)
 * @brief Returns the closest distance to the point p
 * @return Type losestd(V2d<Type> p) : closest distance to p
 * @param[in] p : a point on the plane
 *
 */
template <typename Type> Type Triangle<Type>::closestd(V2d<Type> p) const {
  auto sds = this->sides();
  return std::min_element(
             sds.begin(), sds.end(),
             [&p](auto a, auto b) { return a.distance(p) < b.distance(p); })
      ->distance(p);
}

/**
 * method bool Triangle<Type>::in(V2d<Type> p)
 * @brief Returns if the point p is in the triangle or not
 * @return bool in(V2d<Type> p) : point p is in (true : in, false : out)
 * @param[in] p : a point on the plane
 *
 */
template <typename Type>
bool Triangle<Type>::in(V2d<Type> p) const // Is p in?
{
  if ((this->sides()[0].ends()[1] - p) % (this->sides()[0].ends()[0] - p) >
          0.0 &&
      (this->sides()[1].ends()[1] - p) % (this->sides()[1].ends()[0] - p) >
          0.0 &&
      (this->sides()[2].ends()[1] - p) % (this->sides()[2].ends()[0] - p) > 0.0)
    return true;
  if ((this->sides()[0].ends()[1] - p) % (this->sides()[0].ends()[0] - p) <
          0.0 &&
      (this->sides()[1].ends()[1] - p) % (this->sides()[1].ends()[0] - p) <
          0.0 &&
      (this->sides()[2].ends()[1] - p) % (this->sides()[2].ends()[0] - p) < 0.0)
    return true;
  return false;
}

template <typename Type>
bool Triangle<Type>::in(Circle<Type> cl) const {
  if(!this->in(cl.cent())) return false;
  std::vector<Type> dst;
  for(auto sd : this->sides()) dst.push_back(sd.distance(cl.cent()));
  return dst[0] > cl.r() && dst[1] > cl.r() && dst[2] > cl.r();
}

template <typename Type>
bool Triangle<Type>::in(Circle<Type> cl, Type fact) const {
  Circle<Type> clf(cl.cent(), cl.r() * fact);
  return this->in(clf);
}

template <typename Type>
bool Triangle<Type>::in(Ellipse<Type> el) const{
  std::array<V2d<Type>, 3> eqvap;
  for(size_t i=0; i<c_apex.size(); i++) eqvap[i] = el.eqvCirc01(c_apex[i]);
  Triangle<Type> tri(eqvap);
  Circle<Type> circ(V2d<Type>(0.0, 0.0), 1.0);
  return tri.in(circ);
}

template <typename Type>
bool Triangle<Type>::in(Ellipse<Type> el, Type fact) const{
  Ellipse<Type> elf(el.cent(), el.ab().x(), el.ab().y(), el.ph());
  return this->in(elf);
}

/**
 * method bool Triangle<Type>::insideof(Circle<Type> cl)
 * brief Returns if the triangle is inside of the circle
 * @return bool insideof(Circle<Type> cl) : the triangle is inside of 
 * the circle (true : inside, false : not inside)
 * @param[in] cl : a circle on the plane
 *
 */
template <typename Type>
bool Triangle<Type>::insideof(
    Circle<Type> cl) const // Is this triangle inside of th cl?
{
  return (c_apex[0] - cl.cent()).norm() < cl.r() &&
         (c_apex[1] - cl.cent()).norm() < cl.r() &&
         (c_apex[2] - cl.cent()).norm() < cl.r();
}

/**
 * method Triangle<float> Triangle<Type>::to_f()
 * @brief returns converted triangle to Triangle<float> (32-bit floating point)
 * @return Triangle<float> to_f() : converted Triangle to Triangle<float>
 *
 */
template <typename Type> Triangle<float> Triangle<Type>::to_f() const {
  return Triangle<float>(
      {c_apex[0].to_f(), c_apex[1].to_f(), c_apex[2].to_f()});
}

/**
 * method Triangle<double> Triangle<Type>::to_d()
 * @brief returns converted triangle to Triangle<double> (64-bit floating point)
 * @return Triangle<double> to_d() : converted Triangle to Triangle<double>
 *
 */
template <typename Type> Triangle<double> Triangle<Type>::to_d() const {
  return Triangle<double>(
      {c_apex[0].to_d(), c_apex[1].to_d(), c_apex[2].to_d()});
}

/**
 * method Triangle<long double> Triangle<Type>::to_l()
 * @brief returns converted triangle to Triangle<long double> (128-bit floating point)
 * @return Triangle<long double> to_l() : converted Triangle to Triangle<long double>
 *
 */
template <typename Type> Triangle<long double> Triangle<Type>::to_l() const {
  return Triangle<long double>(
      {c_apex[0].to_l(), c_apex[1].to_l(), c_apex[2].to_l()});
}

/**
 * method void Triangle<Type>::print()
 * @brief prints triangle
 *
 */
template <typename Type> void Triangle<Type>::print() const {
  for (size_t i = 0; i < c_apex.size(); i++) {
    std::cout << "i = " << i << ", ";
    c_apex[i].print();
  }
}

//
//  Methods for Rectangle
//
// template <typename Type>
// Rectangle<Type>::Rectangle(V2d<Type> min, V2d<Type> max) {
//   c_min = min;
//   c_max = max;
// }

/**
 * method V2d<Type> Rectangle<Type>::min()
 * @brief Returns lower left corner
 * @return V2d<Type> min() : lower left corner 
 *
 */
template <typename Type> V2d<Type> Rectangle<Type>::min() const {
  return c_min;
}

/**
 * method V2d<Type> Rectangle<Type>::max()
 * @brief Returns upper right corner
 * @return V2d<Type> max() : upper right corner 
 *
 */
template <typename Type> V2d<Type> Rectangle<Type>::max() const {
  return c_max;
}

/**
 * method V2d<Type> Rectangle<Type>::cent()
 * @brief Returns center of the rectangle
 * @return V2d<Type> cent() : center
 *
 */
template <typename Type> V2d<Type> Rectangle<Type>::cent() const {
  return (c_min + c_max) / 2.0;
}

/**
 * method Type Rectangle<Type>::d()
 * @brief Returns width
 * @return Type d() : width
 *
 */
template <typename Type> Type Rectangle<Type>::d() const {
  return c_max.x() - c_min.x();
}

/**
 * method bool Rectangle<Type>::in()
 * brief Returns if point p is in or not
 * @return bool in() : p is in (true : yes, false : no)
 * @param[in] p : a point on the plane
 *
 */
template <typename Type> bool Rectangle<Type>::in(V2d<Type> p) const {
  return p.x() > c_min.x() && p.x() < c_max.x() && p.y() > c_min.y() &&
         p.y() < c_max.y();
}

/**
 * method Type Rectangle<Type>::area()
 * @brief Returns area
 * @return Type area() : area
 *
 */
template <typename Type> Type Rectangle<Type>::area() const {
  return (c_max.x() - c_min.x()) * (c_max.y() - c_min.y());
}

/**
 * method std::array<V2d<Type>, 4> Rectangle<Type>::corners()
 * @brief Returns corners
 * @return std::array<V2d<Type>, 4> corners() : four corners
 *
 */
template <typename Type>
std::array<V2d<Type>, 4> Rectangle<Type>::corners() const {
  return {c_min, V2d<Type>(c_max.x(), c_min.y()), c_max,
          V2d<Type>(c_min.x(), c_max.y())};
}

/**
 * method std::array<Seg2d<Type>, 4> Rectangle<Type>::segs()
 * @brief Returns segments
 * @return std::array<Seg2d<Type>, 4> segs() : four segments
 *
 */
template <typename Type>
std::array<Seg2d<Type>, 4> Rectangle<Type>::segs() const {
  return {Seg2d<Type>({this->corners()[0], this->corners()[1]}),
          Seg2d<Type>({this->corners()[1], this->corners()[2]}),
          Seg2d<Type>({this->corners()[2], this->corners()[3]}),
          Seg2d<Type>({this->corners()[3], this->corners()[0]})};
}

/**
 * method Rectangle<float> Rectangle<Type>::to_f()
 * @brief returns converted rectangle to Rectangle<float> (32-bit floating point)
 * @return Rectangle<float> to_f() : converted Rectangle to Rectangle<float>
 *
 */
template <typename Type> Rectangle<float> Rectangle<Type>::to_f() const {
  return Rectangle<float>(c_min.to_f(), c_max.to_f());
}

/**
 * method Rectangle<double> Rectangle<Type>::to_d()
 * @brief returns converted rectangle to Rectangle<double> (64-bit floating point)
 * @return Rectangle<double> to_d() : converted Rectangle to Rectangle<double>
 *
 */
template <typename Type> Rectangle<double> Rectangle<Type>::to_d() const {
  return Rectangle<double>(c_min.to_d(), c_max.to_d());
}

/**
 * method Rectangle<long double> Rectangle<Type>::to_l()
 * @brief returns converted rectangle to Rectangle<long double> (128-bit floating point)
 * @return Rectangle<long double> to_l() : converted Rectangle to Rectangle<long double>
 *
 */
template <typename Type> Rectangle<long double> Rectangle<Type>::to_l() const {
  return Rectangle<long double>(c_min.to_l(), c_max.to_l());
}

/**
 * method void Rectangle<Type>::print()
 * @brief prints rectangle
 *
 */
template <typename Type> void Rectangle<Type>::print() const {
  std::cout << "min = ";
  c_min.print();
  std::cout << "max = ";
  c_max.print();
}

/**
 * method V2d<Type> Square<Type>::cent()
 * @brief Returns center of the square
 * @return V2d<Type> cent() : center of the square
 *
 */
template <typename Type> V2d<Type> Square<Type>::cent() const { return c_cent; }

/**
 * method Type Square<Type>::d()
 * @brief Returns the half length of the segments
 * @return Type d() : the half length of the segments
 *
 */
template <typename Type> Type Square<Type>::d() const { return c_d; }

/**
 * method Square<Type> Square<Type>::enlarge(Type fact)
 * @brief Returns the expanded square by a factor of fact
 * @return Square<Type> enlarge(Type fact) : the expanded square by a factor of fact
 *
 */
template <typename Type> Square<Type> Square<Type>::enlarge(Type fact) const {
  return Square<Type>(c_cent, c_d * fact);
}

/**
 * method std::array<Triangle<Type>, 2> Square<Type>::to_tri()
 * @brief Returns the pair of the splitted triangles
 * @return std::array<Triangle<Type>, 2> to_tri() : the pair of the splitted triangles
 *
 */
template <typename Type>
std::array<Triangle<Type>, 2> Square<Type>::to_tri() const {
  return {
      Triangle<Type>({this->min(), V2d<Type>(this->max().x(), this->min().y()),
                      this->max()}),
      Triangle<Type>({this->max(), V2d<Type>(this->min().x(), this->max().y()),
                      this->min()})};
}

/**
 * Square<Type> include_all(std::vector<Square<Type>> sqs)
 * @brief Return a square that include all of the squares sqs
 * @return Square<Type> include_all() : a square that include all of the squares sqs
 * @param[in] sqs : a set of squares
 *
 */
template <typename Type>
Square<Type> include_all(std::vector<Square<Type>> sqs) {
  Type xmax =
      std::max_element(sqs.begin(), sqs.end(),
                       [](auto a, auto b) { return a.max().x() < b.max().x(); })
          ->max()
          .x();
  Type ymax =
      std::max_element(sqs.begin(), sqs.end(),
                       [](auto a, auto b) { return a.max().y() < b.max().y(); })
          ->max()
          .y();
  Type xmin =
      std::min_element(sqs.begin(), sqs.end(),
                       [](auto a, auto b) { return a.min().x() < b.min().x(); })
          ->min()
          .x();
  Type ymin =
      std::min_element(sqs.begin(), sqs.end(),
                       [](auto a, auto b) { return a.min().y() < b.min().y(); })
          ->min()
          .y();
  Type d = std::max(xmax - xmin, ymax - ymin) / 2.0;
  V2d<Type> cent((xmax + xmin) / 2.0, (ymax + ymin) / 2.0);
  return Square<Type>(cent, d);
}

/**
 * method V2d<Type> Square<Type>::cent()
 * @brief Returns center of the square
 * @return V2d<Type> cent() : center of the square
 *
 */
template <class Type> V2d<Type> Circle<Type>::cent() const { return c_cent; }

/**
 * method Type Circle<Type>::r()
 * @brief Returns radious
 * @return Type r() : radious of the circle
 *
 */
template <class Type> Type Circle<Type>::r() const { return c_r; }

/**
 * method Type Circle<Type>::rho(V2d<Type> p)
 * @brief Returns distance from the center of the circle to the point p
 * @return Type rho() : distance from the center
 *
 */
template <class Type> Type Circle<Type>::rho(V2d<Type> p) const {
  return (p - c_cent).norm() / c_r;
}

/**
 * method Type Circle<Type>::area()
 * @brief Returns area of the circle
 * @return Type area() : Area of the circle
 *
 */
template <class Type> Type Circle<Type>::area() const {
  return M_PI * c_r * c_r;
}

/**
 * method Type Circle<Type>::cross_a(Seg2d<Type> seg)
 * @brief Returns area of the circle
 * @return Type area() : Area of the circle
 *
 */
template <class Type> bool Circle<Type>::in(V2d<Type> p) const {
  return (p - c_cent).norm() < c_r;
}

/**
 * method Type Circle<Type>::cross_a(Seg2d<Type> seg)
 * @brief Returns the ratio a to get the cross point
 * @return Type cross_a() : Returns the ratio a to get the cross point
 *
 */
template <typename Type> Type Circle<Type>::cross_a(Seg2d<Type> seg) const {
  auto a0 = seg.ends()[0] - c_cent;
  auto a1 = seg.ends()[1] - c_cent;
  Type a = (a1 - a0).n2();
  Type b = 2.0 * (a0 * (a1 - a0));
  Type c = a0 * a0 - c_r * c_r;
  Type ca0 = (-b + std::sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
  Type ca1 = (-b - std::sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
  // std::cout << "ca0 = " << ca0 << ", ca1 = " << ca1 << std::endl;
  if (ca0 >= 0.0 && ca0 <= 1.0)
    return ca0;
  else if (ca1 >= 0.0 && ca1 <= 1.0)
    return ca1;
  return 0.0;
}

/**
 * method TriInOut<Type> Circle<Type>::inout(Triangle<Type> tri)
 * @brief Returns "in" and "out" apexes
 * @return TriInOut<Type> Circle<Type>::inout() : "in" and "out" apexes
 *
 */
template <typename Type>
TriInOut<Type> Circle<Type>::inout(Triangle<Type> tri) const {
  auto aps = tri.apex();
  TriInOut<Type> io;
  for (auto ap : aps) {
    if (this->in(ap))
      io.in.push_back(ap);
    else
      io.out.push_back(ap);
  }
  return io;
}

/**
 * method V2d<Type> Circle<Type>::cross(Seg2d<Type> seg)
 * @brief Returns the cross point of the segment
 * @return V2d<Type> Circle<Type>::cross() : Cross point
 *
 */
template <typename Type> V2d<Type> Circle<Type>::cross(Seg2d<Type> seg) const {
  return seg.ends()[0] + (seg.ends()[1] - seg.ends()[0]) * this->cross_a(seg);
}

/**
 * method Type Circle<Type>::area_th(Type th)
 * @brief Returns the area of the crescent defined by th
 * @return Type Circle<Type>::area_th() : Area of the crescent
 *
 */
template <typename Type> Type Circle<Type>::area_th(Type th) const {
  if (th > std::pow(std::numeric_limits<Type>::epsilon(), 1.0 / 3.0))
    return c_r * c_r * (th - std::sin(th) * std::cos(th));
  return 2.0 / 3.0 * th * th * th * c_r * c_r;
}

/**
 * method AreaTri<Type> Circle<Type>::areatri(Triangle<Type> tri)
 * @brief Returns the areas of "in", "out", "crescent" regions
 * @return AreaTri<Type> Circle<Type>::areatri() : the areas of "in", "out", "crescent" regions
 *
 */
template <typename Type>
AreaTri<Type> Circle<Type>::areatri(Triangle<Type> tri) const {
  auto io = this->inout(tri);
  if (io.in.size() < 1)
    return AreaTri<Type>({0.0, 0.0, tri.area()});
  else if (io.in.size() == 1) {
    Seg2d<Type> sg0({io.in[0], io.out[0]});
    Seg2d<Type> sg1({io.in[0], io.out[1]});
    V2d<Type> cr0 = this->cross(sg0);
    V2d<Type> cr1 = this->cross(sg1);
    Type th = std::asin(0.5 * (cr1 - cr0).norm() / c_r);
    Triangle<Type> tri_in({io.in[0], cr0, cr1});
    return AreaTri<Type>({tri_in.area() + this->area_th(th), this->area_th(th),
                          tri.area() - tri_in.area() - this->area_th(th)});
  } else if (io.in.size() == 2) {
    Seg2d<Type> sg0({io.in[0], io.out[0]});
    Seg2d<Type> sg1({io.in[1], io.out[0]});
    V2d<Type> cr0 = this->cross(sg0);
    V2d<Type> cr1 = this->cross(sg1);
    Type th = std::asin(0.5 * (cr1 - cr0).norm() / c_r);
    Triangle<Type> tri_out({io.out[0], cr0, cr1});
    return AreaTri<Type>({tri.area() - tri_out.area() + this->area_th(th),
                          this->area_th(th),
                          tri_out.area() - this->area_th(th)});
  }
  return AreaTri<Type>({tri.area(), 0.0, 0.0});
}

/**
 * method Type Circle<Type>::ratio_ia(Triangle<Type> tri)
 * @brief Returns the ratio of the "in" area and all area of the triangle
 * @return Type Circle<Type>::ratio_ia() : the ratio of the "in" area and all area of the triangle
 *
 */
template <typename Type> Type Circle<Type>::ratio_ia(Triangle<Type> tri) const {
  auto at = this->areatri(tri);
  return at.in / tri.area();
}

/**
 * method Type Circle<Type>::swell(Type th)
 * @brief Returns the thickness of te swell defined by th
 * @return Type Circle<Type>::swell() : the thickness of te swell
 *
 */
template <typename Type> Type Circle<Type>::swell(Type th) const {
  if (th >= 3.0 * std::sqrt(std::numeric_limits<Type>::epsilon()))
    return 1.0 - std::cos(th);
  return 0.5 * th * th * c_r;
}

/**
 * method Type Circle<Type>::swell(Seg2d<Type> seg)
 * @brief Returns the thickness of te swell defined by the segment
 * @return Type Circle<Type>::swell() : the thickness of te swell
 *
 */
template <typename Type> Type Circle<Type>::swell(Seg2d<Type> seg) const {
  return this->swell(seg.angle(c_cent));
}

/**
 * method Type Circle<Type>::swell(Triangle<Type> tri)
 * @brief Returns the thickness of te swell defined by the triangle
 * @return Type Circle<Type>::swell() : the thickness of te swell
 *
 */
template <typename Type> Type Circle<Type>::swell(Triangle<Type> tri) const {
  auto sgs = tri.sides();
  std::vector<Type> swl;
  for (size_t i = 0; i < sgs.size(); i++)
    swl.push_back(this->swell(sgs[i]));
  return *std::max_element(swl.begin(), swl.end());
}

/**
 * method Type Circle<Type>::angle(Triangle<Type> tri)
 * @brief Returns the angle of the triangle 
 * @return Type Circle<Type>::angle() : the angle of the triangle 
 *
 */
template <typename Type> Type Circle<Type>::angle(Triangle<Type> tri) const {
  auto sds = tri.sides();
  return std::max_element(
             sds.begin(), sds.end(),
             [&](auto a, auto b) { return a.angle(c_cent) < b.angle(c_cent); })
      ->angle(c_cent);
}

/**
 * method Type Circle<Type>::gap(Triangle<Type> tri)
 * @brief Returns the gap between the sides of the triangle and the circl
 * @return Type Circle<Type>::gap() : 
 *
 */
template <typename Type> Type Circle<Type>::gap(Triangle<Type> tri) const {
  auto sds = tri.sides();
  return std::min_element(sds.begin(), sds.end(),
                          [&](auto a, auto b) {
                            return a.distance(c_cent) < b.distance(c_cent);
                          })
             ->distance(c_cent) -
         c_r;
}

/**
 * method Square<Type> Circle<Type>::include()
 * @brief Returns the smallest square that include the circle
 * @return Square<Type> Circle<Type>::include() : the smallest square that include the circle
 *
 */
template <typename Type> Square<Type> Circle<Type>::include() const {
  return Square<Type>(c_cent, c_r);
}

/**
 * method Circle<float> Circle<Type>::to_f()
 * @brief returns the converted circle to Circle<float> (32-bit floating point)
 * @return Circle<float> Circle<Type>::to_f() : the converted circle to Circle<float> (32-bit floating point)
 *
 */
template <typename Type> Circle<float> Circle<Type>::to_f() const {
  return Circle<float>(c_cent.to_f(), static_cast<float>(c_r));
}

/**
 * method Circle<double> Circle<Type>::to_d()
 * @brief returns the converted circle to Circle<double> (64-bit floating point)
 * @return Circle<double> Circle<Type>::to_d() : the converted circle to Circle<double> (64-bit floating point)
 *
 */
template <typename Type> Circle<double> Circle<Type>::to_d() const {
  return Circle<double>(c_cent.to_d(), static_cast<double>(c_r));
}

/**
 * method Circle<long double> Circle<Type>::to_l()
 * @brief returns the converted circle to Circle<long double> (128-bit floating point)
 * @return Circle<long double> Circle<Type>::to_l() : the converted circle to Circle<double> (128-bit floating point)
 *
 */
template <typename Type> Circle<long double> Circle<Type>::to_l() const {
  return Circle<long double>(c_cent.to_l(), static_cast<long double>(c_r));
}

/**
 * method void Circle<Type>::print()
 * @brief prints circle
 *
 */
template <typename Type> void Circle<Type>::print() const {
  std::cout << "cent = ";
  c_cent.print();
  std::cout << "r = " << c_r << std::endl;
}

/**
 * method V2d<Type> Ellipse<Type>::cent()
 * @brief Returns the center of the ellipse
 * @return V2d<Type> Ellipse<Type>::cent() : Center
 *
 */
template <class Type> V2d<Type> Ellipse<Type>::cent() const { return c_cent; }

/**
 * method V2d<Type> Ellipse<Type>::ab()
 * @brief Returns the Semi-major and semi-minor axes
 * @return V2d<Type> Ellipse<Type>::ab() : Semi-major and semi-minor axes
 *
 */
template <class Type> V2d<Type> Ellipse<Type>::ab() const { return c_ab; }

/**
 * method Type Ellipse<Type>::ph()
 * @brief Returns the inclination
 * @return Type Ellipse<Type>::ph() : Inclination
 *
 */
template <typename Type> Type Ellipse<Type>::ph() const { return c_ph; }

/**
 * method Type Ellipse<Type>::r(V2d<Type> p)
 * @brief Returns the projected radious to the direction of the point p
 * @return Type Ellipse<Type>::r() : the projected radious to the direction of the point p
 *
 */
template <typename Type> Type Ellipse<Type>::r(V2d<Type> p) const {
  V2d<Type> dp = (p - c_cent) ^ (-c_ph);
  Type cs = dp.x() / dp.norm();
  Type sn = dp.y() / dp.norm();
  return std::sqrt(c_ab.x() * c_ab.x() * cs * cs +
                   c_ab.y() * c_ab.y() * sn * sn);
}

/**
 * method Type Ellipse<Type>::rho(V2d<Type> p)
 * @brief Returns the normalize elliptical distance to the point p
 * @return Type Ellipse<Type>::rho() : the normalized elliptical distance of the point p
 *
 */
template <typename Type> Type Ellipse<Type>::rho(V2d<Type> p) const {
  V2d<Type> dp = (p - c_cent) ^ (-c_ph);
  return std::sqrt((dp.x() / c_ab.x()) * (dp.x() / c_ab.x()) +
                   (dp.y() / c_ab.y()) * (dp.y() / c_ab.y()));
}

/**
 * method Type Ellipse<Type>::area()
 * @brief Returns the area of the ellipse
 * @return Type Ellipse<Type>::area() : the area
 *
 */
template <class Type> Type Ellipse<Type>::area() const {
  return M_PI * c_ab.x() * c_ab.y();
}

/**
 * method Type Ellipse<Type>::dist(V2d<Type> p)
 * @brief Returns the normalize elliptical distance to the point p
 * @return Type Ellipse<Type>::dist(V2d<Type> p) : the normalized elliptical distance of the point p
 *
 */
template <class Type> Type Ellipse<Type>::dist(V2d<Type> p) const {
  V2d<Type> ep = (p - c_cent) ^ (-c_ph);
  return std::sqrt(ep.x() * ep.x() / (c_ab.x() * c_ab.x()) +
                   ep.y() * ep.y() / (c_ab.y() * c_ab.y()));
}

/**
 * method bool Ellipse<Type>::in(V2d<Type> p)
 * @brief Returns if  the point p is inside or not inside
 * @return bool Ellipse<Type>::in() : the point p is inside or not inside
 *
 */
template <class Type> bool Ellipse<Type>::in(V2d<Type> p) const {
  return this->dist(p) < 1.0;
}

/**
 * method V2d<Type> Ellipse<Type>::eqvCirc01(V2d<Type> p)
 * @brief Returns the equivalent point of p for the unit circle
 * @return V2d<Type> Ellipse<Type>::eqvCirc01() : the equivalent point
 *
 */
template <class Type> V2d<Type> Ellipse<Type>::eqvCirc01(V2d<Type> p) const {
  V2d<Type> p0 = p ^ (-c_ph);
  return V2d<Type>(p0.x() / c_ab.x(), p0.y() / c_ab.y());
}

/**
 * method Seg2d<Type> Ellipse<Type>::eqvCirc01(Seg2d<Type> seg)
 * @brief Returns the equivalent segment for the unit circle
 * @return Seg2d<Type> Ellipse<Type>::eqvCirc01() : the equivalent segment
 *
 */
template <class Type>
Seg2d<Type> Ellipse<Type>::eqvCirc01(Seg2d<Type> seg) const {
  return Seg2d<Type>(
      {this->eqvCirc01(seg.ends()[0]), this->eqvCirc01(seg.ends()[1])});
}

/**
 * method Type Ellipse<Type>::cross_a(Seg2d<Type> seg)
 * @brief Returns the equivalent segment facotr (cross_a) for the unit circle
 * @return Type Ellipse<Type>::cross_a() : the unit circle equivalent segment factor
 *
 */
template <typename Type> Type Ellipse<Type>::cross_a(Seg2d<Type> seg) const {
  Seg2d<Type> eqsg = this->eqvCirc01(seg);
  Circle<Type> circ(V2d<Type>(0.0, 0.0), 1.0);
  return circ.cross_a(eqsg);
}

/**
 * method V2d<Type> Ellipse<Type>::cross(Seg2d<Type> seg)
 * @brief Returns the cross point of the segment
 * @return V2d<Type> Ellipse<Type>::cross() : the cross point of the segment
 *
 */
template <typename Type> V2d<Type> Ellipse<Type>::cross(Seg2d<Type> seg) const {
  Seg2d<Type> eqsg = this->eqvCirc01(seg);
  Circle<Type> circ(V2d<Type>(0.0, 0.0), 1.0);
  V2d<Type> crc = circ.cross(eqvCirc01(eqsg));
  return (V2d<Type>(crc.x() * c_ab.x(), crc.y() * c_ab.y()) ^ c_ph) + c_cent;
}

/**
 * method Type Ellipse<Type>::ratio_ia(Triangle<Type> tri)
 * @brief Returns the ratio of the "in" area and all area of the triangle
 * @return Type Ellipse<Type>::ratio_ia() : the ratio of the "in" area and all area of the triangle
 *
 */
template <typename Type>
Type Ellipse<Type>::ratio_ia(Triangle<Type> tri) const {
  Triangle<Type> eqvtri({this->eqvCirc01(tri.apex()[0]),
                         this->eqvCirc01(tri.apex()[1]),
                         this->eqvCirc01(tri.apex()[2])});
  Circle<Type> circ(V2d<Type>(0.0, 0.0), 1.0);
  return circ.ratio_ia(eqvtri);
}

/**
 * method Square<Type> Ellipse<Type>::include()
 * @brief Returns the smallest square that include the circle
 * @return Square<Type> Ellipse<Type>::include() : the smallest square that include the circle
 *
 */
template <typename Type> Square<Type> Ellipse<Type>::include() const {
  Type d = std::max(c_ab.x(), c_ab.y());
  return Square<Type>(c_cent, d);
}

/**
 * method Ellipse<float> Circle<Type>::to_f()
 * @brief returns the converted ellipse to Ellipse<float> (32-bit floating point)
 * @return Ellipse<float> Ellipse<Type>::to_f() : the converted circle to Ellipse<float> (32-bit floating point)
 *
 */
template <typename Type> Ellipse<float> Ellipse<Type>::to_f() const {
  return Ellipse<float>(c_cent.to_f(), c_ab.to_f().x(), c_ab.to_f().y(),
                        static_cast<float>(c_ph));
}

/**
 * method Ellipse<double> Ellipse<Type>::to_d()
 * @brief returns the converted ellipse to Ellipse<double> (64-bit floating point)
 * @return Ellipse<double> Ellipse<Type>::to_d() : the converted ellipse to Ellipse<double> (64-bit floating point)
 *
 */
template <typename Type> Ellipse<double> Ellipse<Type>::to_d() const {
  return Ellipse<double>(c_cent.to_d(), c_ab.to_d().x(), c_ab.to_d().y(),
                         static_cast<double>(c_ph));
}

/**
 * method Ellipse<long double> Ellipse<Type>::to_l()
 * @brief returns the converted ellipse to Ellipse<long double> (128-bit floating point)
 * @return Ellipse<long double> Ellipse<Type>::to_l() : the converted Ellipse to Ellipse<double> (128-bit floating point)
 *
 */
template <typename Type> Ellipse<long double> Ellipse<Type>::to_l() const {
  return Ellipse<long double>(c_cent.to_l(), c_ab.to_l().x(), c_ab.to_l().y(),
                              static_cast<long double>(c_ph));
}

/**
 * method void Ellipse<Type>::print()
 * @brief prints ellipse
 *
 */
template <typename Type> void Ellipse<Type>::print() const {
  std::cout << "cent = ";
  c_cent.print();
  std::cout << "ab = ";
  c_ab.print();
  std::cout << "ph = " << c_ph << std::endl;
}

/**
 * V2d<To> conv(V2d<From> v)
 * @brief Convert from V2d<From> to V2d<To> ; From and To : float, dlouble, long double
 * @return V2d<To> conv() : Converted V2d to To
 * @param[in] v : V2d to convert
 *
 */
template <typename To, typename From> V2d<To> conv(V2d<From> v) {
  return V2d<To>(static_cast<To>(v.x()), static_cast<To>(v.y()));
}

/**
 * M2d<To> conv(M2d<From> m)
 * @brief Convert from M2d<From> to M2d<To> ; From and To : float, dlouble, long double
 * @return M2d<To> conv() : Converted M2d from From to To
 * @param[in] m : M2d to convert
 *
 */
template <typename To, typename From> M2d<To> conv(M2d<From> m) {
  return M2d<To>(static_cast<To>(m.xx()), static_cast<To>(m.xy()),
                 static_cast<To>(m.yx()), static_cast<To>(m.yy()));
}

/**
 * Rectangle<To> conv(Rectangle<From> rec)
 * @brief Convert from Rectangle<From> to Rectangle<To> ; From and To : float, dlouble, long double
 * @return Rectangle<To> conv() : Converted Rectangle from From to To
 * @param[in] rec : Rectangle to convert
 *
 */
template <typename To, typename From> Rectangle<To> conv(Rectangle<From> rec) {
  return Rectangle<To>(conv<To, From>(rec.min()), conv<To, From>(rec.max()));
}

/**
 * Circle<To> conv(Circle<From> cl)
 * @brief Convert from Circle<From> to Circle<To> ; From and To : float, dlouble, long double
 * @return Circle<To> conv() : Converted Circle from From to To
 * @param[in] cl : Circle to convert
 *
 */
template <typename To, typename From> Circle<To> conv(Circle<From> cl) {
  return Circle<To>(conv<To, From>(cl.cent()), static_cast<To>(cl.r()));
}

/**
 * Ellipse<To> conv(Ellipse<From> el)
 * @brief Convert from Ellipse<From> to Ellipse<To> ; From and To : float, dlouble, long double
 * @return Ellipse<To> conv() : Converted Ellipse from From to To
 * @param[in] el : Ellipse to convert
 *
 */
template <typename To, typename From> Ellipse<To> conv(Ellipse<From> el) {
  return Ellipse<To>(conv<To, From>(el.cent()), static_cast<To>(el.ab().x()),
                     static_cast<To>(el.ab().y()), static_cast<To>(el.ph()));
}

// Explicit instantiations of templates of classes
// For float
template class Angle<float>;
template class RA<float>;
template class Dec<float>;
template class Longitude<float>;
template class Latitude<float>;
template class V2d<float>;
template class Seg2d<float>;
template class M2d<float>;
template class Triangle<float>;
template class Rectangle<float>;
template class Square<float>;
template class Circle<float>;
template class Ellipse<float>;
template class V3d<float>;
template class M3d<float>;
template class Rotate_x<float>;
template class Rotate_y<float>;
template class Rotate_z<float>;
template class Align_x<float>;
template class Align_y<float>;
template class Align_z<float>;
template class DirCos<float>;

//For double
template class Angle<double>;
template class RA<double>;
template class Dec<double>;
template class Longitude<double>;
template class Latitude<double>;
template class V2d<double>;
template class Seg2d<double>;
template class M2d<double>;
template class Triangle<double>;
template class Rectangle<double>;
template class Square<double>;
template class Circle<double>;
template class Ellipse<double>;
template class V3d<double>;
template class M3d<double>;
template class Rotate_x<double>;
template class Rotate_y<double>;
template class Rotate_z<double>;
template class Align_x<double>;
template class Align_y<double>;
template class Align_z<double>;
template class DirCos<double>;

//For long double
template class Angle<long double>;
template class RA<long double>;
template class Dec<long double>;
template class Longitude<long double>;
template class Latitude<long double>;
template class V2d<long double>;
template class Seg2d<long double>;
template class M2d<long double>;
template class Triangle<long double>;
template class Rectangle<long double>;
template class Square<long double>;
template class Circle<long double>;
template class Ellipse<long double>;
template class V3d<long double>;
template class M3d<long double>;
template class Rotate_x<long double>;
template class Rotate_y<long double>;
template class Rotate_z<long double>;
template class Align_x<long double>;
template class Align_y<long double>;
template class Align_z<long double>;
template class DirCos<long double>;

// For conv(), include_all()
template V2d<float> conv<float>(V2d<float>);
template V2d<float> conv<float>(V2d<double>);
template V2d<float> conv<float>(V2d<long double>);
template V2d<double> conv<double>(V2d<float>);
template V2d<double> conv<double>(V2d<double>);
template V2d<double> conv<double>(V2d<long double>);
template V2d<long double> conv<long double>(V2d<float>);
template V2d<long double> conv<long double>(V2d<double>);
template V2d<long double> conv<long double>(V2d<long double>);

template Square<float> include_all(std::vector<Square<float>> sps);
template Square<double> include_all(std::vector<Square<double>> sps);
template Square<long double> include_all(std::vector<Square<long double>> sps);
