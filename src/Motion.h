/**
 * @file Motion.h
 *
 * @brief header file for motions
 * @author Fumio Abe
 *
 */

// define _INC_Motion_h_ to load only once
#ifndef _INC_Motion_h_
#define _INC_Motion_h_

// include headers
#include "Geom.h"
#include "Lens.h"
#include "Source.h"

/**
 * template class Time
 * @brief time in verious units
 *
 * Type : double,  double,  long double
 * Aliases : fTime (Time<double),
 *           dTime (Time<double),
 *           lTime (Time<long double)
 *
 */
template <class Type> class Time {
protected:
  /// @brief The time in day, protected
  Type c_t;

public:
  /// @brief Blank constructor
  Time() : c_t(0.0) {};
  /// @brief Constructor to specify time in day
  /// @param[in] day : time in day
  Time(Type day) : c_t(day) {};
  /// @brief Change time in second
  /// @param[in] sec : time in second
  void from_sec(Type sec) { c_t = sec / (60.0 * 60.0 * 24.0); };
  /// @brief Change time in minute
  /// @param[in] min : time in minute
  void from_min(Type min) { c_t = min / (60.0 * 24.0); };
  /// @brief Change time in hour
  /// @param[in] hour : time in hour
  void from_hour(Type hour) { c_t = hour / 24.0; };
  /// @brief Change time in day
  /// @param[in] day : time in day
  void from_day(Type day) { c_t = day; };
  /// @brief Change time in year
  /// @param[in] year : time in year
  void from_year(Type year) { c_t = year * 365.24218944; };
  /// @brief Return time in second
  /// @return Type sec() : time in second
  Type sec() const { return c_t * (60.0 * 60.0 * 24.0); };
  /// @brief Return time in minute
  /// @return Type min() : time in minute
  Type min() const { return c_t * (60.0 * 24.0); };
  /// @brief Return time in hour
  /// @return Type hour() : time in hour
  Type hour() const { return c_t * 24.0; };
  /// @brief Return time in day
  /// @return Type day() : time in day
  Type day() const { return c_t; };
  /// @brief Return time in year
  /// @return Type year() : time in year
  Type year() const { return c_t / 365.24218944; };
  /// @brief Return a series of Time between this and Time<Type> t
  /// @param[in] n : number of divisions
  /// @param[in] t : maximum end of time
  /// @return std::vector<Time<Type>> range() : a series of Time<Type>
  std::vector<Time<Type>> range(int n, Time<Type> t) const {
    Type dt = (t.day() - c_t) / static_cast<Type>(n);
    std::vector<Time<Type>> rt;
    for (int i = 0; i < n; i++)
      rt.push_back(Time<Type>(c_t + dt * i));
    return rt;
  };
  /// @brief Return a series of Time between this and Type t (in day)
  /// @param[in] n : number of divisions
  /// @param[in] t : maximum end of time in day
  /// @return std::vector<Time<Type>> range() : a series of Time<Type>
  std::vector<Time<Type>> range(int n, Type t) const {
    Type dt = (t - c_t) / static_cast<Type>(n);
    std::vector<Time<Type>> rt;
    for (int i = 0; i < n; i++)
      rt.push_back(Time<Type>(c_t + dt * i));
    return rt;
  };
  /// @brief Return a series of Time between t1 and t2
  /// @param[in] n : number of divisions
  /// @param[in] t1 : minimum end of time
  /// @param[in] t2 : maximum end of time
  /// @return std::vector<Time<Type>> range() : a series of Time<Type>
  std::vector<Time<Type>> range(int n, Time<Type> t1, Time<Type> t2) const {
    Type dt = (t2.day() - t1.day()) / static_cast<Type>(n);
    std::vector<Time<Type>> rt;
    for (int i = 0; i < n; i++)
      rt.push_back(Time<Type>(t1.day() + dt * i));
    return rt;
  };
  /// @brief Return a series of Time between t1 and t2
  /// @param[in] n : number of divisions
  /// @param[in] t1 : minimum end of time in day
  /// @param[in] t2 : maximum end of time in day
  /// @return std::vector<Time<Type>> range() : a series of Time<Type>
  std::vector<Time<Type>> range(int n, Type t1, Type t2) const {
    Type dt = (t2 - t1) / static_cast<Type>(n);
    std::vector<Time<Type>> rt;
    for (int i = 0; i < n; i++)
      rt.push_back(Time<Type>(t1 + dt * i));
    return rt;
  };
  /// @brief Add time
  /// @param[in] tt : time to add
  /// @return Time<Type> operator+(Time tt) : result
  Time<Type> operator+(Time tt) { return Time<Type>(c_t + tt.day()); };
  /// @brief Subtract time
  /// @param[in] tt : time to add
  /// @return Time<Type> operator-(Time tt) : result
  Time<Type> operator-(Time tt) { return Time<Type>(c_t - tt.day()); };
  /// @brief Convert Time to double (32-bit floating point)
  /// @return Time<double> to_f() :  converted Time to Time<float>
  Time<float> to_f() const { return Time<float>(static_cast<float>(c_t)); };
  /// @brief Convert Time to double (64-bit doubleing point)
  /// @return Time<double> to_d() :  converted Time to Time<double>
  Time<double> to_d() const { return Time<double>(static_cast<double>(c_t)); }
  /// @brief Convert Time to long double (128-bit doubleing point)
  /// @return Time<long double> to_l() :  converted Time to Time<long double>
  Time<long double> to_l() const {
    return Time<long double>(static_cast<long double>(c_t));
  }
  /// @brief Print parameters
  void print() const { std::cout << "t = " << c_t << " days" << std::endl; }
};

/**
 * template class JD245
 * @brief time in Jukian day - 2450000.0
 *
 * Type : double,  double,  long double
 * Aliases : fJD245 (JD245<double),
 *           dJD245 (JD245<double),
 *           lJD245 (JD245<long double)
 *
 */
template <class Type> class JD245 : public Time<Type> {
private:
public:
  /// @brief Constructor
  /// @param[in] jd245 : Julian day - 2450000.0
  JD245(Type jd245) : Time<Type>(jd245) {};
  /// @brief Constructor in calendar day
  /// @param[in] year : calendar year
  /// @param[in] month : calendar month
  /// @param[in] day : calendar day
  /// @param[in] hour : calendar hour
  /// @param[in] min : calendar minute
  /// @param[in] sec : calendar second
  JD245(Type year, Type month, Type day, Type hour, Type min, Type sec)
      : Time<Type>(std::floor(365.25 * (year - 1.0)) +
                   std::floor((year - 1.0) / 400.0) -
                   std::floor((year - 1.0) / 100.0) +
                   std::floor(30.59 * ((month + 12.0) - 2.0)) + day +
                   (hour + 12.0) / 24.0 + min / (24.0 * 60.0) +
                   sec / (24.0 * 60.0 * 60.0) - 678912 - 50000.0) {};
  /// @brief Return Julian day
  /// @return Type jd() : Julian day
  Type jd() const { return this->c_t - 2450000.0; };
  /// @brief Return MJD (Julian day - 2400000.0)
  /// @return Type mjd() : MJD
  Type mjd() const { return this->c_t + 50000.0; };
  /// @brief Convert JD245 to double (32-bit floating point)
  /// @return JD245<double> to_f() :  converted JD245 to JD245<float>
  JD245<double> to_f() const {
    return JD245<double>(static_cast<float>(this->c_t));
  };
  /// @brief Convert JD245 to double (64-bit doubleing point)
  /// @return JD245<double> to_d() :  converted JD245 to JD245<double>
  JD245<double> to_d() const {
    return JD245<double>(static_cast<double>(this->c_t));
  };
  /// @brief Convert JD245 to long double (128-bit doubleing point)
  /// @return JD245<long double> to_l() :  converted JD245 to JD245<long double>
  JD245<long double> to_l() const {
    return JD245<long double>(static_cast<long double>(this->c_t));
  };
};

/**
 * template class Linear
 * @brief 2d constant velocity linear motion
 *
 * Type : double,  double,  long double
 * Aliases : fLinear (Linear<double),
 *           dLinear (Linear<double),
 *           lLinear (Linear<long double)
 *
 */
template <class Type> class Linear {
private:
  /// @brief 2d position at t0 (ordinary closest position)
  V2d<Type> c_p0;
  /// @brief Time 0 : ordinaly closest time
  Time<Type> c_t0;
  /// @brief 2d velocity
  V2d<Type> c_v;

public:
  /// @brief Blank constructor : (0.0, 0.1) at t0 = 0.0, v = (0.1, 0.0)
  Linear()
      : c_p0(V2d<Type>(0.0, 0.1)), c_t0(Time<Type>(0.0)),
        c_v(V2d<Type>(0.1, 0.0)) {};
  /// @brief Constructor
  /// @param[in] p0 : position at origin
  /// @param[in] t0 : time of origin
  /// @param[in] v :  velocity
  Linear(V2d<Type> p0, Time<Type> t0, V2d<Type> v)
      : c_p0(p0), c_t0(t0), c_v(v) {};
  /// @brief Constructor
  /// @param[in] p0 : position at origin
  /// @param[in] t0 : time of origin (in day)
  /// @param[in] v :  velocity
  Linear(V2d<Type> p0, Type t0, V2d<Type> v)
      : c_p0(p0), c_t0(Time(t0)), c_v(v) {};
  /// @brief Constructor
  /// @param[in] u0 : closest distance
  /// @param[in] t0 : closest time
  /// @param[in] phi : slope of the orbit in x-y plane
  /// @param[in] te : Einstein radius crossing time
  Linear(Type u0, Time<Type> t0, Angle<Type> phi, Time<Type> te)
      : c_p0(V2d<Type>(u0 * std::cos(phi.rad() + M_PI / 2.0),
                       -u0 * std::sin(phi.rad() + M_PI / 2.0))),
        c_t0(t0), c_v(V2d<Type>(1.0 / te.day() * std::cos(phi.rad()),
                                1.0 / te.day() * std::sin(phi.rad()))) {};
  /// @brief Constructor
  /// @param[in] u0 : closest distance
  /// @param[in] t0 : closest time
  /// @param[in] phi : slope of the orbit in x-y plane in radian
  /// @param[in] te : Einstein radius crossing time
  Linear(Type u0, Time<Type> t0, Type phi, Time<Type> te)
      : c_p0(V2d<Type>(u0 * std::cos(phi + M_PI / 2.0),
                       -u0 * std::sin(phi + M_PI / 2.0))),
        c_t0(t0), c_v(V2d<Type>(1.0 / te.day() * std::cos(phi),
                                1.0 / te.day() * std::sin(phi))) {};
  /// @brief Return p0 : position at origin
  /// @return V2d<Type> p0() : position at y0 (closest position)
  V2d<Type> p0() const;
  /// @brief Return t0 : position at origin
  /// @return Time<Type> t0() : time at p0
  Time<Type> t0() const;
  /// @brief Return v : 2d vector
  /// @return V2d<Type> v() : 2d velocity
  V2d<Type> v() const;
  /// @brief Return u0 : closest distance
  /// @return Type u0() : closest distance
  Type u0() const;
  /// @brief Return phi : slope in x-y coordinate
  /// @return Angle<Type> phi() : slope of the orbit
  Angle<Type> phi() const;
  /// @brief Return te : Einstein radious crossing time
  /// @return Time<Type> te() : slope of the orbit
  Time<Type> te() const;
  /// @brief Return 2d position at a time t
  /// @param[in] t
  /// @return V2d<Type> p(Time<Type> t) : 2d position at a time t
  V2d<Type> p(Time<Type> t) const;
  /// @brief Return distance at a time t
  /// @param[in] t
  /// @return Type u(Time<Type> t) : distance at a time t
  Type u(Time<Type> t) const;
  /// @brief Convert Linear to double (32-bit floating point)
  /// @return Linear<double> to_f() :  converted Linear to Linear<float>
  Linear<float> to_f() const {
    return Linear<float>(c_p0.to_f(), c_t0.to_f(), c_v.to_f());
  };
  /// @brief Convert Linear to double (64-bit doubleing point)
  /// @return Linear<double> to_d() :  converted Linear to Linear<double>
  Linear<double> to_d() const {
    return Linear<double>(c_p0.to_d(), c_t0.to_d(), c_v.to_d());
  };
  /// @brief Convert Linear to long double (128-bit doubleing point)
  /// @return Linear<long double> to_l() :  converted Linear to Linear<long
  /// double>
  Linear<long double> to_l() const {
    return Linear<long double>(c_p0.to_l(), c_t0.to_l(), c_v.to_l());
  };
  /// @brief Print parameters
  void print() const;
};

/**
 * template class Length
 * @brief Length in verious units
 *
 * Type : double,  double,  long double
 * Aliases : fLength (Length<double),
 *           dLength (Length<double),
 *           lLength (Length<long double)
 *
 */
template <class Type> class Length {
private:
  /// @brief Length in astronomical unit
  Type c_au;

public:
  /// @brief Blank constructor : 0.0 AU
  Length() : c_au(0.0) {};
  /// @brief Constructor
  /// @param[in] au : length in au (distance between Sun and Earth)
  Length(Type au) : c_au(au) {};
  /// @brief Init length in meter
  /// @param[in] m :  length in meter
  void from_m(Type m) { c_au = m / 149597870700.0; };
  /// @brief Init length in kilo meter
  /// @param[in] km :  length in kilo meter
  void from_km(Type km) { this->from_m(km * 1000.0); };
  /// @brief Init length in astronomical unit
  /// @param[in] au :  length in astronomical unit
  void from_au(Type au) { c_au = au; };
  /// @brief Init length in persec
  /// @param[in] pc :  length in persec
  void from_pc(Type pc) { this->from_au(pc * 60.0 * 60.0 * 180.0 / M_PI); };
  /// @brief Return length in meter
  /// @return Type m() :  length in meter
  Type m() const { return c_au * 149597870700.0; };
  /// @brief Return length in kilo meter
  /// @return Type km() :  length in kilo meter
  Type km() const { return this->m() / 1000.0; };
  /// @brief Return length in astronomical unit
  /// @return Type au() :  length in astronomical unit
  Type au() const { return c_au; };
  /// @brief Return length in persec
  /// @return Type pc() :  length in persec
  Type pc() const { return c_au / (60.0 * 60.0 * 180.0 / M_PI); };
  /// @brief Convert Length to double (32-bit floating point)
  /// @return Length<double> to_f() :  converted Length to Length<float>
  Length<float> to_f() const { return Length<float>(static_cast<float>(c_au)); }
  /// @brief Convert Length to double (64-bit doubleing point)
  /// @return Length<double> to_d() :  converted Length to Length<double>
  Length<double> to_d() const {
    return Length<double>(static_cast<double>(c_au));
  }
  /// @brief Convert Length to long double (128-bit doubleing point)
  /// @return Length<long double> to_l() :  converted Length to Length<long
  /// double>
  Length<long double> to_l() const {
    return Length<long double>(static_cast<long double>(c_au));
  }
  /// @brief Print parameters
  void print() const { std::cout << "au = " << c_au << std::endl; }
};

/**
 * template class V2dU
 * @brief 2d vector with unit
 *
 * Type : double,  double,  long double
 * Aliases : fV2dU (V2dU<double),
 *           dV2dU (V2dU<double),
 *           lV2dU (V2dU<long double)
 *
 */
template <class Type> class V2dU : public V2d<Type> {
private:
public:
  /// @brief Blank constructor (0.0, 0.0)
  V2dU() : V2d<Type>(0.0, 0.0) {};
  /// @brief Constructor in astronomical unit
  /// @param[in] x : x in AU
  /// @param[in] y : y in AU
  V2dU(Type x, Type y) : V2d<Type>(x, y) {};
  /// @brief Constructor in astronomical unit vector
  /// @param[in] vv : 2d vector in astronomical unit
  V2dU(V2d<Type> vv) : V2d<Type>(vv.x(), vv.y()) {};
  /// @brief Constructor in Length
  /// @param[in] x : x in Length
  /// @param[in] y : y in Length
  V2dU(Length<Type> x, Length<Type> y) : V2d<Type>(x.au(), y.au()) {};
  /// @brief Return 2d vector in astronomical unit
  /// @return V2d<Type> v2d() : 2d vector in AU
  V2d<Type> v2d() const { return *this; };
  /// @brief Return 2d vector in meter
  /// @return V2d<Type> m() : 2d vector in meter
  V2d<Type> m() const {
    Length<Type> x(V2d<Type>::x());
    Length<Type> y(V2d<Type>::y());
    return V2d<Type>(x.m(), y.m());
  };
  /// @brief Return 2d vector in kilo meter
  /// @return V2d<Type> km() : 2d vector in kilo meter
  V2d<Type> km() const {
    Length<Type> x(V2d<Type>::x());
    Length<Type> y(V2d<Type>::y());
    return V2d<Type>(x.km(), y.km());
  };
  /// @brief Return 2d vector in astronomical unit
  /// @return V2d<Type> au() : 2d vector in AU
  V2d<Type> au() const {
    Length<Type> x(V2d<Type>::x());
    Length<Type> y(V2d<Type>::y());
    return V2d<Type>(x.au(), y.au());
  };
  /// @brief Return 2d vector in persec
  /// @return V2d<Type> pc() : 2d vector in persec
  V2d<Type> pc() const {
    Length<Type> x(V2d<Type>::x());
    Length<Type> y(V2d<Type>::y());
    return V2d<Type>(x.pc(), y.pc());
  };
  /// @brief Add a vector
  /// @param[in] vv : the vector to add
  /// @return V2dU<Type> operator+(V2dU<Type> vv) : result
  V2dU<Type> operator+(V2dU<Type> vv) const {
    return V2dU<Type>(this->v2d() + vv.v2d());
  };
  /// @brief Subtract a vector
  /// @param[in] vv : the vector to subtract
  /// @return V2dU<Type> operator-(V2dU<Type> vv) : result
  V2dU<Type> operator-(V2dU<Type> vv) const {
    return V2dU<Type>(this->v2d() - vv.v2d());
  };
  /// @brief Inner product
  /// @param[in] vv : the vector to multiply
  /// @return Length<Type> operator*(V2dU<Type> vv) : result
  Length<Type> operator*(V2dU<Type> vv) const {
    return Length<Type>(this->v2d() * vv.v2d());
  };
  /// @brief Multiply by a number
  /// @param[in] dd : the number to multiply
  /// @return V2dU<Type> operator*(Type dd) : result
  V2dU<Type> operator*(Type dd) const { return V2dU<Type>(this->v2d() * dd); };
  /// @brief Divide by a number
  /// @param[in] dd : the number to divide
  /// @return V2dU<Type> operator/(Type dd) : result
  V2dU<Type> operator/(Type dd) const { return V2dU<Type>(this->v2d() / dd); };
  /// @brief Multiply by a matrix
  /// @param[in] mm : the matrix to multiply
  /// @return V2dU<Type> operator*(M2d<Type> mm) : result
  V2dU<Type> operator*(M2d<Type> mm) const {
    return V2dU<Type>(this->v2d() * mm);
  };
  /// @brief Return if two vectors are equal or not
  /// @param[in] vv : the vector to compare
  /// @return bool operator==(V2dU<Type> vv) : they are equal (true) or not
  /// (false)
  bool operator==(V2dU<Type> vv) const { return this->v2d() == vv.v2d(); };
  /// @brief Convert V2dU to double (32-bit floating point)
  /// @return V2dU<double> to_f() :  converted V2dU to V2dU<float>
  V2dU<float> to_f() const {
    return V2dU<float>(V2d<Type>::to_f().x(), V2d<Type>::to_f().y());
  };
  /// @brief Convert V2dU to double (64-bit doubleing point)
  /// @return V2dU<double> to_d() :  converted V2dU to V2dU<double>
  V2dU<double> to_d() const {
    return V2dU<double>(V2d<Type>::to_d().x(), V2d<Type>::to_d().y());
  };
  /// @brief Convert V2dU to long double (128-bit doubleing point)
  /// @return V2dU<long double> to_l() :  converted V2dU to V2dU<long double>
  V2dU<long double> to_l() const {
    return V2dU<long double>(V2d<Type>::to_l().x(), V2d<Type>::to_l().y());
  };
};

/**
 * template class V3dU
 * @brief 3d vector with unit
 *
 * Type : double,  double,  long double
 * Aliases : fV3dU (V3dU<double),
 *           dV3dU (V3dU<double),
 *           lV3dU (V3dU<long double)
 *
 */
template <class Type> class V3dU : public V3d<Type> {
private:
public:
  /// @brief Blank constructor (0.0, 0.0, 0.0)
  V3dU() : V3d<Type>(0.0, 0.0, 0.0) {};
  /// @brief Constructor in astronomical unit
  /// @param[in] x : x in AU
  /// @param[in] y : y in AU
  /// @param[in] z : z in AU
  V3dU(Type x, Type y, Type z) : V3d<Type>(x, y, z) {};
  /// @brief Constructor in astronomical unit vector
  /// @param[in] vv : 3d vector in astronomical unit
  V3dU(V3d<Type> vv) : V3d<Type>(vv.x(), vv.y(), vv.z()) {};
  /// @brief Constructor in Length
  /// @param[in] x : x in Length
  /// @param[in] y : y in Length
  /// @param[in] z : z in Length
  V3dU(Length<Type> x, Length<Type> y, Length<Type> z)
      : V3d<Type>(x.au(), y.au(), z.au()) {};
  /// @brief Return 3d vector in astronomical unit
  /// @return V3d<Type> v2d() : 3d vector in AU
  V3d<Type> v3d() const { return *this; };
  /// @brief Return 3d vector in meter
  /// @return V3d<Type> m() : 3d vector in meter
  V3d<Type> m() const {
    Length<Type> x(V3d<Type>::x());
    Length<Type> y(V3d<Type>::y());
    Length<Type> z(V3d<Type>::z());
    return V3d<Type>(x.m(), y.m(), z.m());
  };
  /// @brief Return 3d vector in kilo meter
  /// @return V3d<Type> km() : 3d vector in kilo meter
  V3d<Type> km() const {
    Length<Type> x(V3d<Type>::x());
    Length<Type> y(V3d<Type>::y());
    Length<Type> z(V3d<Type>::z());
    return V3d<Type>(x.km(), y.km(), z.km());
  };
  /// @brief Return 3d vector in astronomical unit
  /// @return V3d<Type> au() : 3d vector in AU
  V3d<Type> au() const {
    Length<Type> x(V3d<Type>::x());
    Length<Type> y(V3d<Type>::y());
    Length<Type> z(V3d<Type>::z());
    return V3d<Type>(x.au(), y.au(), z.au());
  };
  /// @brief Return 3d vector in persec
  /// @return V3d<Type> pc() : 3d vector in persec
  V3d<Type> pc() const {
    Length<Type> x(V3d<Type>::x());
    Length<Type> y(V3d<Type>::y());
    Length<Type> z(V3d<Type>::z());
    return V3d<Type>(x.pc(), y.pc(), z.pc());
  };
  /// @brief Add a vector
  /// @param[in] vv : the vector to add
  /// @return V3dU<Type> operator+(V3dU<Type> vv) : result
  V3dU<Type> operator+(V3dU<Type> vv) const {
    return V3dU<Type>(this->v3d() + vv.v3d());
  };
  /// @brief Subtract a vector
  /// @param[in] vv : the vector to subtract
  /// @return V3dU<Type> operator-(V3dU<Type> vv) : result
  V3dU<Type> operator-(V3dU<Type> vv) const {
    return V3dU<Type>(this->v3d() - vv.v3d());
  };
  /// @brief Inner product
  /// @param[in] vv : the vector to multiply
  /// @return Length<Type> operator*(V3dU<Type> vv) : result
  Length<Type> operator*(V3dU<Type> vv) const {
    return Length<Type>(this->v3d() * vv.v3d());
  };
  /// @brief Multiply by a number
  /// @param[in] dd : the number to multiply
  /// @return V3dU<Type> operator*(Type dd) : result
  V3dU<Type> operator*(Type dd) const { return V3dU<Type>(this->v3d() * dd); };
  /// @brief Divide by a number
  /// @param[in] dd : the number to divide
  /// @return V2dU<Type> operator/(Type dd) : result
  V3dU<Type> operator/(Type dd) const { return V3dU<Type>(this->v3d() / dd); };
  /// @brief Multiply by a matrix
  /// @param[in] mm : the matrix to multiply
  /// @return V3dU<Type> operator*(M3d<Type> mm) : result
  V3dU<Type> operator*(M3d<Type> mm) const {
    return V3dU<Type>(mm * this->v3d());
  };
  /// @brief Return if two vectors are equal or not
  /// @param[in] vv : the vector to compare
  /// @return bool operator==(V3dU<Type> vv) : they are equal (true) or not
  /// (false)
  bool operator==(V3dU<Type> vv) const { return this->v3d() == vv.v3d(); };
  /// @brief Convert V3dU to double (32-bit floating point)
  /// @return V3dU<double> to_f() :  converted V3dU to V3dU<float>
  V3dU<double> to_f() const {
    return V3dU<double>(V3d<Type>::to_f().x(), V3d<Type>::to_f().y(),
                        V3d<Type>::to_f().z());
  };
  /// @brief Convert V3dU to double (64-bit doubleing point)
  /// @return V3dU<double> to_d() :  converted V3dU to V2dU<double>
  V3dU<double> to_d() const {
    return V3dU<double>(V3d<Type>::to_d().x(), V3d<Type>::to_d().y(),
                        V3d<Type>::z());
  };
  /// @brief Convert V3dU to long double (128-bit doubleing point)
  /// @return V3dU<long double> to_l() :  converted V3dU to V2dU<long double>
  V3dU<long double> to_l() const {
    return V3dU<long double>(V3d<Type>::to_l().x(), V3d<Type>::to_l().y(),
                             V3d<Type>::to_l().z());
  };
};

/**
 * template class EllipseU
 * @brief 2d ellipse with unit
 *
 * Type : double,  double,  long double
 * Aliases : fEllipseU (EllipseU<double),
 *           dEllipseU (EllipseU<double),
 *           lEllipseU (EllipseU<long double)
 *
 */
template <class Type> class EllipseU : public Ellipse<Type> {
private:
public:
  /// @brief Blank constructor. All null.
  EllipseU() : Ellipse<Type>(V2d<Type>(0.0, 0.0), 0.0, 0.0) {};
  /// @brief Constructor with cent, a, b, and ph
  /// @param[in] cent : center of the ellipse
  /// @param[in] a : Semi-mjor axis
  /// @param[in] b : semi-minor axis
  /// @param[in] ph : inclination
  EllipseU(V2dU<Type> cent, Length<Type> a, Length<Type> b, Angle<Type> ph)
      : Ellipse<Type>(cent.au(), a.au(), b.au(), ph.rad()) {};
  /// @brief Constructor with ell
  /// @param[in] ell : Ellipse
  EllipseU(Ellipse<Type> ell) : Ellipse<Type>(ell) {};
  /// @brief Return Ellipse without unit (AU)
  /// @return Ellipse<Type> v2d() : Ellipse with no dimension (AU)
  Ellipse<Type> v2d() const { return *this; };
  /// @brief Return Ellipse in meter
  /// @return Ellipse<Type> m() : Ellipse in meter
  Ellipse<Type> m() const {
    V2dU<Type> cent = V2dU<Type>(Ellipse<Type>::cent());
    V2dU<Type> ab = V2dU<Type>(Ellipse<Type>::ab());
    Type ph = Ellipse<Type>::ph();
    return Ellipse(cent.m(), ab.m().x(), ab.m().y(), ph);
  };
  /// @brief Return Ellipse in Astronomical Unit
  /// @return Ellipse<Type> au() : Ellipse in Astronomical Unit
  Ellipse<Type> au() const {
    V2dU<Type> cent = V2dU<Type>(Ellipse<Type>::cent());
    V2dU<Type> ab = V2dU<Type>(Ellipse<Type>::ab());
    Type ph = Ellipse<Type>::ph();
    return Ellipse(cent.au(), ab.au().x(), ab.au().y(), ph);
  };
  /// @brief Return Ellipse in persec
  /// @return Ellipse<Type> pc() : Ellipse in persec
  Ellipse<Type> pc() const {
    V2dU<Type> cent = V2dU<Type>(Ellipse<Type>::cent());
    V2dU<Type> ab = V2dU<Type>(Ellipse<Type>::ab());
    Type ph = Ellipse<Type>::ph();
    return Ellipse(cent.pc(), ab.pc().x(), ab.pc().y(), ph);
  };
  /// @brief Convert to double
  /// @return EllipseU<double> to_f() : converted EllipseU to float (32-bit
  /// floating point)
  EllipseU<float> to_f() const {
    return EllipseU<float>(Ellipse<Type>::to_f());
  };
  /// @brief Convert to double
  /// @return EllipseU<double> to_d() : converted EllipseU to double (64-bit
  /// doubleing point)
  EllipseU<double> to_d() const {
    return EllipseU<double>(Ellipse<Type>::to_d());
  };
  /// @brief Convert to long double
  /// @return EllipseU<long double> to_l() : converted EllipseU to long double
  /// (128-bit doubleing point)
  EllipseU<long double> to_l() const {
    return EllipseU<long double>(Ellipse<Type>::to_l());
  };
};

/// @brief struct to express Kepler orbit
template <typename Type> struct Orbit {
  /// @brief Center of mass
  V2dU<Type> cm = V2dU<Type>(0.0, 0.0);
  /// @brief Semi-mjor axis
  Length<Type> a = Length<Type>(0.0);
  /// @brief Eccetricity
  Type e = 0.0;
  /// @brief Period
  Time<Type> T = Time<Type>(0.0);
  /// @brief Epoc of perihelion
  Time<Type> tp = Time<Type>(0.0);
  /// @brief Argument of perihelion
  Angle<Type> omega = Angle<Type>(0.0);
  /// @brief Inclination
  Angle<Type> incl = Angle<Type>(0.0);
  /// Longitude of ascending
  Angle<Type> Omega = Angle<Type>(0.0);
};

/**
 * template class Kepler
 * @brief class for Kepler orbit
 *
 * Type : double,  double,  long double
 * Aliases : fKepler (Kepler<double),
 *           dKepler (Kepler<double),
 *           lKepler (Kepler<long double)
 *
 */
template <class Type> class Kepler {
protected:
  /// @brief Orbit parameters
  Orbit<Type> c_obt;

public:
  /// @brief Blank constructor
  Kepler() {};
  Kepler(V2dU<Type> cm, Type a, Type e, Time<Type> T, Time<Type> tp) {
    /// brief Constructor with cm, a, e, T, and tp
    /// @param[in] cm : center of mass
    /// @param[in] a : semi-minor axis
    /// @param[in] e : Eccetricity
    /// @param[in] T : Period
    /// @param[in] tp : Epoc of perihelion
    c_obt = {
        cm, a, e, T, tp, Angle<Type>(0.0), Angle<Type>(0.0), Angle<Type>(0.0)};
  };
  /// @brief Constructor with a, e, T, and tp
  /// @param[in] a : semi-minor axis
  /// @param[in] e : Eccetricity
  /// @param[in] T : Period
  /// @param[in] tp : Epoc of perihelion
  Kepler(Type a, Type e, Time<Type> T, Time<Type> tp) {
    c_obt = {
        V2dU<Type>(0.0, 0.0), a, e, T, tp, Angle<Type>(0.0), Angle<Type>(0.0),
        Angle<Type>(0.0)};
  };
  /// @brief Constructor with a, e, T, tp, omega, incl, and Omega
  /// @param[in] a : semi-minor axis
  /// @param[in] e : Eccetricity
  /// @param[in] T : Period
  /// @param[in] tp : Epoc of perihelion
  /// @param[in] omega : Argument of perihelion
  /// @param[in] incl : Inclination
  /// @param[in] Omega : Longitude of ascending
  Kepler(Length<Type> a, Type e, Time<Type> T, Time<Type> tp, Angle<Type> omega,
         Angle<Type> incl, Angle<Type> Omega) {
    c_obt = {V2dU<Type>(0.0, 0.0), a, e, T, tp, omega, incl, Omega};
  };
  /// @brief Constructor with obt
  /// @param[in] obt : All orbit parameters
  Kepler(Orbit<Type> obt) : c_obt(obt) {};
  /// @brief Constructor with name
  /// @param obj : name of the object (ex. "Moon")
  Kepler(std::string obj);
  /// @brief Return orbit parameters
  /// @return Orbit<Type> orbit() : All orbit parameters
  Orbit<Type> orbit() const { return c_obt; };
  /// @brief Return center of mass
  /// @return Center of mass
  V2dU<Type> cent() const {
    return V2dU<Type>(c_obt.cm.x() + c_obt.a.au() * c_obt.e, c_obt.cm.y());
  };
  /// @brief Return Semi-mjor axis
  /// @return Length<Type> b() : semi-minor axis
  Length<Type> b() const {
    return Length<Type>(c_obt.a.au() * std::sqrt(1.0 - c_obt.e * c_obt.e));
  };
  /// @brief Return EllipseU parameters
  /// @return EllipseU<Type> ell() : EllipsU parameters
  EllipseU<Type> ell() const {
    return EllipseU<Type>(c_obt.cm, c_obt.a, this->b(), c_obt.omega);
  };
  /// @brief Return mean motion
  /// @return  : Mean motion (Angle/day)
  Angle<Type> mm() const { return Angle<Type>(2.0 * M_PI / c_obt.T.day()); };
  /// @brief Return mean anomaly
  /// @param[in] t : the time
  /// @return Angle<Type> ma(Time<Type> t) : Mean anomaly of the time
  Angle<Type> ma(Time<Type> t) const {
    return Angle<Type>(this->mm().rad() * (t.day() - c_obt.tp.day()));
  };
  /// @brief Return initial value of average periasis angle
  /// @param[in] t : the time
  /// @return Angle<Type> apa_i(Time<Type> t) : Initial value of Average
  /// periasis angle
  Angle<Type> apa_i(Time<Type> t) const;
  /// @brief Return difference between mean anomaly and eccentric anomaly
  /// @param[in] apa : eccentric anomaly
  /// @param[in] t : the time
  /// @return Angle<Type> dma(Angle<Type> apa, Time<Type> t) : difference
  /// between mean anomaly and eccentric anomaly
  Angle<Type> dma(Angle<Type> apa, Time<Type> t) const {
    return Angle<Type>(this->ma(t).rad() -
                       (apa.rad() - c_obt.e * std::sin(apa.rad())));
  };
  /// @brief Return average periasis angle
  /// @param[in] t : the time
  /// @param[in] eps : precision
  /// @return Angle<Type> apa(Time<Type> t, Type eps) : average periasis angle
  Angle<Type> apa(Time<Type> t, Type eps) const {
    Angle<Type> apa = this->apa_i(t);
    if (this->dma(apa, t) < eps)
      return apa;
    while (this->dma(apa, t) / apa >= eps) {
      apa += this->dma(apa, t) / (1.0 - c_obt.e * std::cos(apa.rad()));
    }
    return apa;
  }
  /// @brief Return average periasis angle at the time
  /// @param[in] t : the time
  /// @return Angle<Type> apa(Time<Type> t, Type eps) : average periasis angle
  Angle<Type> apa(Time<Type> t) const {
    return this->apa(t, 4.0 * std::numeric_limits<Type>::epsilon());
  };
  /// @brief Return 2d position at time t
  /// @param[in] t : a time
  /// @return V2dU<Type> p2d(Time<Type> t) : 2d position at time t
  V2dU<Type> p2d(Time<Type> t) const {
    return V2dU<Type>(c_obt.a.au() * (std::cos(this->apa(t).rad()) - c_obt.e),
                      c_obt.a.au() * std::sqrt(1.0 - c_obt.e * c_obt.e) *
                          std::sin(this->apa(t).rad()));
  };
  /// @brief Return two dimensional position at the time
  /// @param[in] t : the time
  /// @return V2dU<Type> cm2dFromp2d(Time<Type> t) : two dimensional position
  V2dU<Type> cm2dFromp2d(Time<Type> t) const {
    return c_obt.cm - this->p2d(t);
  };
  /// @brief Return three dimensional position at the time
  /// @param[in] t : the time
  /// @return V3d<Type> p3d(Time<Type> t) : theree dimensional position
  V3d<Type> p3d(Time<Type> t) const {
    V3d<Type> cm(c_obt.cm.x(), c_obt.cm.y(), 0.0);
    return V3d<Type>(this->p2d(t).x(), this->p2d(t).y(), 0.0) - cm;
  };
  /// @brief Return three dimensional center of mass position relative to the p3
  /// at the time
  /// @param[in] t : the time
  /// @return VV3d<Type> cm3dFromp3d(Time<Type> t) : theree dimensional cm
  /// position
  V3d<Type> cm3dFromp3d(Time<Type> t) const {
    V3d<Type> cm(c_obt.cm.x(), c_obt.cm.y(), 0.0);
    return cm - this->p3d(t);
  };
  /// @brief Return three dimensional position in the ecliptic space at the time
  /// @param[in] t : the time
  /// @return V3d<Type> EC(Time<Type> t) : three dimensional position in the
  /// ecliptic space
  V3d<Type> EC(Time<Type> t) const {
    return this->p3d(t)
        .rotate_z(c_obt.omega)
        .rotate_x(c_obt.incl)
        .rotate_z(c_obt.Omega);
  };
  /// @brief Return three dimensional cm position in the ecliptic space at the
  /// time
  /// @param[in] t : the time
  /// @return V3d<Type> ECcm(Time<Type> t) : three dimensional cm position in
  /// the ecliptic space
  V3d<Type> ECcm(Time<Type> t) const {
    return this->cm3dFromp3d(t)
        .rotate_z(c_obt.omega)
        .rotate_x(c_obt.incl)
        .rotate_z(c_obt.Omega);
  };
  /// @brief Return three dimensional position in the equatorial space at the
  /// time
  /// @param[in] t : the time
  /// @return V3d<Type> EQ(Time<Type> t) : three dimensional position in the
  /// equatorial space
  V3d<Type> EQ(Time<Type> t) const { return this->EC(t).ECtoEQ(); };
  /// @brief Return three dimensional center of mass position in the equatorial
  /// space
  /// @param[in] t : the time
  /// @return V3d<Type> EQ(Time<Type> t) : three dimensional cm position in the
  /// equatorial space
  V3d<Type> EQcm(Time<Type> t) const { return this->ECcm(t).ECtoEQ(); };

  /// @brief Convert to double (32-bit floating point)
  /// @return Kepler<double> to_f() : converted class to float
  Kepler<float> to_f() const {
    Orbit<float> obt = {c_obt.cm.to_f(),
                        static_cast<float>(c_obt.a.au()),
                        static_cast<float>(c_obt.e),
                        c_obt.T.to_f(),
                        c_obt.tp.to_f(),
                        c_obt.omega.to_f(),
                        c_obt.incl.to_f(),
                        c_obt.Omega.to_f()};
    return Kepler<float>(obt);
  };
  /// @brief Convert to double (64-bit doubleing point)
  /// @return Kepler<double> to_d() : converted class to double
  Kepler<double> to_d() const {
    Orbit<double> obt = {c_obt.cm.to_d(),
                         static_cast<double>(c_obt.a.au()),
                         static_cast<double>(c_obt.e),
                         c_obt.T.to_d(),
                         c_obt.tp.to_d(),
                         c_obt.omega.to_d(),
                         c_obt.incl.to_d(),
                         c_obt.Omega.to_d()};
    return Kepler<double>(obt);
  };
  /// @brief Convert to long double (128-bit doubleing point)
  /// @return Kepler<long double> to_l() : converted class to long double
  Kepler<long double> to_l() const {
    Orbit<long double> obt = {c_obt.cm.to_l(),
                              static_cast<long double>(c_obt.a.au()),
                              static_cast<long double>(c_obt.e),
                              c_obt.T.to_l(),
                              c_obt.tp.to_l(),
                              c_obt.omega.to_l(),
                              c_obt.incl.to_l(),
                              c_obt.Omega.to_l()};
    return Kepler<long double>(obt);
  };
};

// template <class Type> class ObtEarth : public Kepler<Type> {
// private:
// public:
//   ObtEarth()
//       : Kepler<Type>({V2d_0<Type>, Length<Type>(1.00), 0.01670,
//       Time<Type>(0.0),
//                       Time<Type>(0.0), Angle<Type>(0.0), Angle<Type>(0.0),
//                       Angle<Type>(0.0)}) {
//     // Time<Type> T;
//     Kepler<Type>::c_obt.tp = JD245<Type>(2024, 01, 03, 00, 39, 0.0);
//     Kepler<Type>::c_obt.T.from_year(1.000017);
//   }
//   ObtEarth<double> to_f() const {
//     return ObtEarth<double>();
//   }
//   ObtEarth<double> to_d() const {
//     return ObtEarth<double>();
//   }
//   ObtEarth<long double> to_l() const {
//     return ObtEarth<long double>();
//   }
// };

/**
 * template class PiE (under construction)
 * @brief class for orbital parallax
 *
 * Type : double,  double,  long double
 * Aliases : fPiE (PiE<double),
 *           dPiE (PiE<double),
 *           lPiE (PiE<long double)
 *
 */
template <class Type> class PiE {
private:
public:
};

/**
 * template class Binary
 * @brief class for Binary orbit
 *
 * Type : double,  double,  long double
 * Aliases : fBinary (Binary<double),
 *           dBinary (Binary<double),
 *           lBinary (Binary<long double)
 *
 */
template <class Type> class Binary {
private:
  /// @brief Maaseses
  std::array<Type, 2> c_m;
  /// @brief Kepler parameters
  std::array<Kepler<Type>, 2> c_kep;

public:
  /// @brief Constructor with kep1, m1, and m2
  /// @param[in] kep1 : Kepler parameter for obj1
  /// @param[in] m1 : Mass of the obi1
  /// @param[in] m2 : Mass of the oby2
  Binary(Kepler<Type> kep1, Type m1, Type m2) : c_m({m1, m2}) {
    c_kep = {kep1,
             Kepler<Type>({kep1.orbit().cm,
                           Length<Type>(kep1.orbit().a.au() * m2 / (m1 + m2)),
                           kep1.orbit().e, kep1.orbit().T, kep1.orbit().tp,
                           kep1.orbit().omega, kep1.orbit().incl,
                           Angle<Type>(kep1.orbit().Omega.rad() + M_PI)})};
  };
  /// @brief Constructor with obt, m1, and m2
  /// @param[in] obt : Kepler parameters
  /// @param[in] m1 : Mass of the obi1
  /// @param[in] m2 : Mass of the oby2
  Binary(Orbit<Type> obt, Type m1, Type m2) : c_m({m1, m2}) {
    c_kep = {Kepler<Type>(Length<Type>(obt.a.au() * m1 / (m1 + m2)), obt.e,
                          obt.T, obt.tp, obt.omega, obt.incl, obt.Omega),
             Kepler<Type>(Length<Type>(obt.a.au() * m2 / (m1 + m2)), obt.e,
                          obt.T, obt.tp, obt.omega, obt.incl,
                          Angle<Type>(obt.Omega.rad() + M_PI))};
  };
  /// @brief Return Kepler parameters
  /// @return std::array<Kepler<Type>, 2> kepler() : Kepler parameters
  std::array<Kepler<Type>, 2> kepler() const { return c_kep; };
  /// @brief Return equatorial 3d positions
  /// @param[in] t : the time
  /// @return std::array<V3d<Type>, 2> EQ(Time<Type> t) : equatorial 3d
  /// positions
  std::array<V3d<Type>, 2> EQ(Time<Type> t) const {
    return {c_kep[0].EQ(t), c_kep[1].EQ(t)};
  };
  /// @brief Return equatorial 2d position
  /// @param[in] t : the time
  /// @return V2d<Type> p(Time<Type> t) : equatorial 2d position
  V2d<Type> p(Time<Type> t) const {
    return {V2d<Type>(this->EQ(t)[0].x(), this->EQ(t)[0].y())};
  };
  /// @brief Return 2d position of the center of mass
  /// @return V2dU<Type> cm() : the center of mass
  V2dU<Type> cm() const { return c_kep[0].orbit().cm; };
  /// @brief Convert to double (32-bit floating point)
  /// @return Binary<double> to_f() : converted class to float
  Binary<float> to_f() const {
    return Binary<float>(c_kep[0].to_f(), static_cast<double>(c_m[0]),
                         static_cast<float>(c_m[1]));
  };
  /// @brief Convert to double (64-bit doubleing point)
  /// @return Binary<double> to_d() : converted class to double
  Binary<double> to_d() const {
    return Binary<double>(c_kep[0].to_d(), static_cast<double>(c_m[0]),
                          static_cast<double>(c_m[1]));
  };
  /// @brief Convert to long double (128-bit doubleing point)
  /// @return Binary<long double> to_l() : converted class to long double
  Binary<long double> to_l() const {
    return Binary<long double>(c_kep[0].to_l(),
                               static_cast<long double>(c_m[0]),
                               static_cast<long double>(c_m[1]));
  };
};

/**
 * template class SrcMotion
 * @brief class for motion of the source star
 *
 * Type : double,  double,  long double
 * Aliases : fSrcMotion (SrcMotion<double),
 *           dSrcMotion (SrcMotion<double),
 *           lSrcMotion (SrcMotion<long double)
 *
 */
template <class Type> class SrcMotion {
private:
  /// @brief Source parameters
  Source<Type> c_src;
  /// @brief Function to express motion
  std::function<V2d<Type>(Time<Type>)> c_f;

public:
  /// @brief Constructor with src
  /// @param[in] src : source parameters
  SrcMotion(Source<Type> src) : c_src(src) {};
  /// @brief Constructor with src and l (linear motion)
  /// @param[in] src : source parameters
  /// @param[in] l : parameters of linear motion
  SrcMotion(Source<Type> src, Linear<Type> l)
      : c_src(src), c_f([l](Time<Type> t) { return l.p(t); }) {};
  /// @brief Constructor with src and f (general motion)
  /// @param[in] src : source parameters
  /// @param[in] f : funtion to express motion
  SrcMotion(Source<Type> src, std::function<V2d<Type>(Time<Type>)> f)
      : c_src(src) {};
  /// @brief Return source parameters
  /// @return Source<Type> src() : source parameters
  Source<Type> src() const { return c_src; };
  /// @brief Return source parameters at the time t
  /// @param[in] t : the time
  /// @return Source<Type> src(Time<Type> t) : source parameters at the time t
  Source<Type> src(Time<Type> t) const {
    Source<Type> src = this->src();
    src.move_to(c_f(t));
    return src;
  }
  std::function<V2d<Type>(Time<Type>)> f() const { return c_f; }
  // SrcMotion<float> to_f() const {
  //   return SrcMotion<float>(c_src, [this](Time<float> t) {
  //     auto tt = Time<Type>(static_cast<Type>(t.day()));
  //     return this->p(t).to_f();
  //   });
  /// @brief Set source
  /// @param[in] src : Source
  void SetSrc(Source<Type> src) { c_src = src; };
  /// @brief Return the 2d position of the source at the time t
  /// @param[in] t : the time
  /// @return V2d<Type> p(Time<Type> t) : 2d position of the source
  V2d<Type> p(Time<Type> t) const;
  /// @brief Return distance from p0 at time t
  /// @param[in] t : the time
  /// @param[in] p0 : a position
  /// @return Type u(Time<Type> t, V2d<Type> p0) : distance from p0
  Type u(Time<Type> t, V2d<Type> p0) const;
  /// @brief Return distance from origin
  /// @param[in] t : the time
  /// @return Type u(Time<Type> t) : distance from the origin
  Type u(Time<Type> t) const;
  /// @brief Return the series of 2d positions of the source
  /// @param[in] t : a series of the time
  /// @return std::vector<V2d<Type>> p(std::vector<Time<Type>> t) : a series of
  /// the 2d positions of the source
  std::vector<V2d<Type>> p(std::vector<Time<Type>> t) const;
  /// @brief Return the series of distances from the position p0
  /// @param[in] t : a series of the time
  /// @param[in] p0 : a position
  /// @return std::vector<Type> u(std::vector<Time<Type>> t, V2d<Type> p0) : the
  /// series of distances
  std::vector<Type> u(std::vector<Time<Type>> t, V2d<Type> p0) const;
  /// @brief Return the series of distances from the origin
  /// @param[in] t : a series of the time
  /// @return std::vector<Type> u(std::vector<Time<Type>> t) : the series of
  /// distances
  std::vector<Type> u(std::vector<Time<Type>> t) const;
  /// @brief Convert to double (32-bit floating point)
  /// @return SrcMotion<double> to_f() : converted class to float
  SrcMotion<float> to_f() const;
  /// @brief Convert to double (64-bit doubleing point)
  /// @return SrcMotion<double> to_d() : converted class to double
  SrcMotion<double> to_d() const;
  /// @brief Convert to long double (128-bit doubleing point)
  /// @return SrcMotion<long double> to_l() : converted class to long double
  SrcMotion<long double> to_l() const;
  /// @brief Print parameters
  void print() const;
};

/**
 * template class LensMotion
 * @brief class for motion of the source star
 *
 * Type : double,  double,  long double
 * Aliases : fLensMotion (LensMotion<double),
 *           dLensMotion (LensMotion<double),
 *           lLensMotion (LensMotion<long double)
 *
 */
template <class Type> class LensMotion {
private:
  /// @brief Lens parameters
  Lens<Type> c_l;
  /// @brief Function to express motion
  std::function<V2d<Type>(Time<Type>)> c_f;

public:
  /// @brief Blank constructor
  LensMotion() {};
  /// @brief Construcyor with a lens
  LensMotion(Lens<Type> l)
      : c_l(l), c_f([l](Time<Type> t) { return l.p(); }) {};
  /// @brief Constructor with l and f
  LensMotion(Lens<Type> l, std::function<V2d<Type>(Time<Type>)> f)
      : c_l(l), c_f(f) {};
  // std::function<V2d<Type>(Time<Type>)> f() const {return c_f;};
  /// @brief Return the lens position at time t
  /// @param[in] t : the time
  /// @return V2d<Type> p(Time<Type> t) : 2d position of the lens at t
  V2d<Type> p(Time<Type> t) const { return c_f(t); };
  /// @brief Return lens parameters
  /// @return Lens<Type> l() : lens parameters
  Lens<Type> l() const { return c_l; };
  /// @brief Return the lens parameters at time t
  /// @param[in] t : the time
  /// @return Lens<Type> l(Time<Type> t) : lens parameters at t
  Lens<Type> l(Time<Type> t) const {
    Lens<Type> l = c_l;
    l.move_to(c_f(t));
    return l;
  };
  /// @brief Return the series of lens parameters at a series of time t
  /// @param[in] t : a series of time
  /// @return std::vector<Lens<Type>> l(std::vector<Time<Type>> t) : a series of
  /// lens parameters at t
  std::vector<Lens<Type>> l(std::vector<Time<Type>> t) const;
  /// @brief Convert to double (32-bit floating point)
  /// @return LensMotion<double> to_f() : converted class to float
  LensMotion<float> to_f() const;
  /// @brief Convert to double (64-bit doubleing point)
  /// @return LensMotion<double> to_d() : converted class to double
  LensMotion<double> to_d() const;
  /// @brief Convert to long double (128-bit doubleing point)
  /// @return LensMotion<long double> to_l() : converted class to long double
  LensMotion<long double> to_l() const;
  /// @brief Print parameters
  void print() const;
};

/**
 * template class MlMotion
 * @brief class for motion of the multiple lens system
 *
 * Type : double,  double,  long double
 * Aliases : fMlMotion (MlMotion<double),
 *           dMlMotion (MlMotion<double),
 *           lLensMotion (MlMotion<long double)
 *
 */
template <class Type> class MlMotion : public std::vector<LensMotion<Type>> {
  using size_type = typename std::vector<LensMotion<Type>>::size_type;

public:
  using std::vector<LensMotion<Type>>::vector;
  /// @brief Square bracket access to each element
  LensMotion<Type> &operator[](size_type i);
  /// @brief Square bracket access to each element
  const LensMotion<Type> &operator[](size_type i) const;
  // MlMotion() : std::vector<LensMotion<Type>>(LensMotion<Type>()) {};
  /// @brief Constructor
  /// @param[in] ml : multiple lens parameters
  MlMotion(Mlens<Type> ml) {
    for (size_t i = 0; i < ml.size(); i++)
      this->push_back(LensMotion(ml[i]));
  };
  /// @brief Constructor with ml and f
  /// @param[in] ml : multiple lens parameters
  /// @param[in] f : array functions to express motions of the lenses
  MlMotion(Mlens<Type> ml,
           std::vector<std::function<V2d<Type>(Time<Type>)>> f) {
    for (size_t i = 0; i < ml.size(); i++)
      this->push_back(LensMotion(ml[i], f[i]));
  };
  /// @return the series of mass ratios
  /// @return std::vector<Type> q() : the series of mass ratio
  std::vector<Type> q() const {
    std::vector<Type> q;
    for (size_t i = 0; i < this->size(); i++)
      q.push_back(this->at(i).l().q());
    return q;
  };
  /// @return the series of positions of the lenses
  /// @return std::vector<V2d<Type>> p(Time<Type> t) : the series of positions
  /// at time t
  std::vector<V2d<Type>> p(Time<Type> t) const {
    std::vector<V2d<Type>> p;
    for (size_t i = 0; i < this->size(); i++)
      p.push_back(this->at(i).p(t));
    return p;
  }
  /// @return the multiple lens parameters at time t
  /// @param[in] t : the time
  /// @return Mlens<Type> ml(Time<Type> t) : multiple lens parameters at time t
  Mlens<Type> ml(Time<Type> t) const {
    Mlens<Type> ml;
    for (size_t i = 0; i < this->size(); i++)
      ml.push_back(this->at(i).l(t));
    return ml;
  };

  Mlens<Type> ml() const {
    Mlens<Type> ml;
    for (size_t i = 0; i < this->size(); i++)
      ml.push_back(this->at(i).l());
    return ml;
  };
  MlMotion<float> to_f() const {
    std::vector<std::function<V2d<float>(Time<float>)>> vf;
    for (size_t i = 0; i < this->size(); i++)
      vf.push_back([this, i](Time<float> tf) {
        auto tt = Time<Type>(static_cast<Type>(tf.day()));
        return this->p(tt)[i].to_f();
      });
    return MlMotion<float>(this->ml().to_f(), vf);
  };
  MlMotion<double> to_d() const {
    std::vector<std::function<V2d<double>(Time<double>)>> vf;
    for (size_t i = 0; i < this->size(); i++)
      vf.push_back([this, i](Time<double> tf) {
        auto tt = Time<Type>(static_cast<Type>(tf.day()));
        return this->p(tt)[i].to_d();
      });
    return MlMotion<double>(this->ml().to_d(), vf);
  };
  MlMotion<long double> to_l() const {
    std::vector<std::function<V2d<long double>(Time<long double>)>> vf;
    for (size_t i = 0; i < this->size(); i++)
      vf.push_back([this, i](Time<long double> tf) {
        auto tt = Time<Type>(static_cast<Type>(tf.day()));
        return this->p(tt)[i].to_l();
      });
    return MlMotion<long double>(this->ml().to_l(), vf);
  };
};

// template <class Type> class Parallax {
// private:
// public:
// };
//
// template <class Type> class Xallarap {
// private:
// public:
// };
//
// template <class Type> class Kepler {
// private:
// public:
// };
//
// template <class Type> class Motion {
// private:
//   Linear<Type> c_in;
//   Parallax<Type> c_par;
//   Xallarap<Type> c_xal;
//   Kepler<Type> c_kep;
// public:
//
// };

//
// Suppress instantations for double
//
extern template class Linear<float>;
extern template class Length<float>;
extern template class V2dU<float>;
extern template class V3dU<float>;
extern template class Time<float>;
extern template class JD245<float>;
extern template class EllipseU<float>;
extern template struct Orbit<float>;
extern template class Kepler<float>;
// extern template class ObtEarth<double>;
extern template class PiE<float>;
extern template class Binary<float>;
extern template class SrcMotion<float>;
extern template class LensMotion<float>;
extern template class MlMotion<float>;

/// @brief Alias of Linear<float>
using fLinear = Linear<float>;
/// @brief Alias of Length<float>
using fLength = Length<float>;
/// @brief Alias of V2dU<float>
using fV2dU = V2dU<float>;
/// @brief Alias of V3dU<double>
using fV3dU = V3dU<float>;
/// @brief Alias of Time<double>
using fTime = Time<float>;
/// @brief Alias of JD245<float>
using fJD245 = JD245<float>;
/// @brief Alias of EllipseU<float>
using fEllipseU = EllipseU<float>;
/// @brief Alias of Orbit<float>
using fOrbit = Orbit<float>;
/// @brief Alias of Kepler<float>
using fKepler = Kepler<float>;
// using fObtEarth = ObtEarth<float>;
/// @brief Alias of PiE<float>
using fPiE = PiE<float>;
/// @brief Alias of Binary<float>
using fBinary = Binary<float>;
/// @brief Alias of SrcMotion<float>
using fSrcMotion = SrcMotion<float>;
/// @brief Alias of LensMotion<float>
using fLensMotion = LensMotion<float>;
/// @brief Alias of MlMotion<float>
using fMlMotion = MlMotion<float>;

//
// Suppress instantations for double
//
extern template class Linear<double>;
extern template class Length<double>;
extern template class V2dU<double>;
extern template class V3dU<double>;
extern template class Time<double>;
extern template class JD245<double>;
extern template class EllipseU<double>;
extern template struct Orbit<double>;
extern template class Kepler<double>;
// extern template class ObtEarth<double>;
extern template class PiE<double>;
extern template class Binary<double>;
extern template class SrcMotion<double>;
extern template class LensMotion<double>;
extern template class MlMotion<double>;

/// @brief Alias of Linear<double>
using dLinear = Linear<double>;
/// @brief Alias of Length<double>
using dLength = Length<double>;
/// @brief Alias of V2dU<double>
using dV2dU = V2dU<double>;
/// @brief Alias of V3dU<double>
using dV3dU = V3dU<double>;
/// @brief Alias of Time<double>
using dTime = Time<double>;
/// @brief Alias of JD245<double>
using dJD245 = JD245<double>;
/// @brief Alias of EllipseU<double>
using dEllipseU = EllipseU<double>;
/// @brief Alias of EllipseU<double>
using dOrbit = Orbit<double>;
/// @brief Alias of Kepler<double>
using dKepler = Kepler<double>;
// using dObtEarth = ObtEarth<double>;
/// @brief Alias of PiE<double>
using dPiE = PiE<double>;
/// @brief Alias of Binary<double>
using dBinary = Binary<double>;
/// @brief Alias of SrcMotion<double>
using dSrcMotion = SrcMotion<double>;
/// @brief Alias of LensMotion<double>
using dLensMotion = LensMotion<double>;
/// @brief Alias of MlMotion<double>
using dMlMotion = MlMotion<double>;

//
// Suppress instantations for long double
//
extern template class Linear<long double>;
extern template class Length<long double>;
extern template class V2dU<long double>;
extern template class V3dU<long double>;
extern template class Time<long double>;
extern template class JD245<long double>;
extern template class EllipseU<long double>;
extern template struct Orbit<long double>;
extern template class Kepler<long double>;
// extern template class ObtEarth<long double>;
extern template class PiE<long double>;
extern template class Binary<long double>;
extern template class SrcMotion<long double>;
extern template class LensMotion<long double>;
extern template class MlMotion<long double>;

/// @brief Alias of Linear<long double>
using lLinear = Linear<long double>;
/// @brief Alias of length<long double>
using lLength = Length<long double>;
/// @brief Alias of V2dU<long double>
using lV2dU = V2dU<long double>;
/// @brief Alias of V3dU<long double>
using lV3dU = V3dU<long double>;
/// @brief Alias of Time<long double>
using lTime = Time<long double>;
/// @brief Alias of JD245<long double>
using lJD245 = JD245<long double>;
/// @brief Alias of EllipseU<long double>
using lEllipseU = EllipseU<long double>;
/// @brief Alias of Orbit<long double>
using lOrbit = Orbit<long double>;
/// @brief Alias of Kepler<long double>
using lKepler = Kepler<long double>;
// using lObtEarth = ObtEarth<long double>;
/// @brief Alias of PiE<long double>
using lPiE = PiE<long double>;
/// @brief Alias of Binary<long double>
using lBinary = Binary<long double>;
/// @brief Alias of SrcMotion<long double>
using lSrcMotion = SrcMotion<long double>;
/// @brief Alias of LensMotion<long double>
using lLensMotion = LensMotion<long double>;
/// @brief Alias of MlMotion<long double>
using lMlMotion = MlMotion<long double>;

#endif
