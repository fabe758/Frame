/**
 * @file Geom.h
 *
 * @brief header file for functions, classes to handle geomeries
 * @author Fumio Abe
 *
 */

// define _INC_Geom_h_ to load only once
#ifndef _INC_Geom_h_
#define _INC_Geom_h_

// include all c++ headers
#include <bits/stdc++.h>
#include <math.h>

//
//  Blank definitions
//
template <class Type> class M2d;
template <class Type> class RA;
template <class Type> class Dec;
template <class Type> class Longitude;
template <class Type> class Latitude;

/**
 * template class Angle
 * @brief angle in verious units
 *
 * Type : float,  double,  long double
 * Aliases : fAngle (Angle<float),
 *           dAngle (Angle<double),
 *           lAngle (Angle<long double)
 *
 */
template <class Type> class Angle {
protected:
  /// @brief The angle in radian, protected
  Type c_rad;

public:
  /// @brief Blank constructor
  Angle() : c_rad(0.0) {};
  /// @brief Constructor to specify angle in radian
  /// @param[in] rad : angle in radian
  Angle(Type rad) : c_rad(rad) {};
  /// @brief Change angle in radian
  /// @param[in] rad : angle in radian
  void from_rad(Type rad) { c_rad = rad; };
  /// @brief Change angle in degree
  /// @param[in] deg : angle in degree
  void from_deg(Type deg) { c_rad = M_PI / 180.0 * deg; };
  /// @brief Change angle in minutes
  /// @param[in] min : angle in minute
  void from_min(Type min) { c_rad = M_PI / 180.0 / 60.0 * min; };
  /// @brief Change angle in seconds
  /// @param[in] sec : angle in second
  void from_sec(Type sec) { c_rad = M_PI / 180.0 / 60.0 / 60.0 * sec; };
  /// @brief Change angle in Right Ascension
  /// @param[in] hh : hour in Right Ascension
  /// @param[in] mm : minute in Right Ascension
  /// @param[in] ss : second in Right Ascension
  void from_ra(Type hh, Type mm, Type ss) {
    c_rad = M_PI / 24.0 * (hh + mm / 60.0 + ss / 60.0 / 60.0);
  };
  /// @brief Change angle in Declination
  /// @param[in] dd : degree in declination
  /// @param[in] mm : minute in declination
  /// @param[in] ss : second in declination
  void from_dec(Type dd, Type mm, Type ss) {
    if (dd >= 0.0)
      c_rad = M_PI / 180.0 * (dd + mm / 60.0 + ss / 60.0 / 60.0);
    else
      c_rad = M_PI / 180.0 * (dd - mm / 60.0 - ss / 60.0 / 60.0);
  };
  /// @brief Change angle in Longitude
  /// @param[in] dd : degree in Longitude
  /// @param[in] mm : minute in Longitude
  /// @param[in] ss : second in Longitude
  /// @param[in] ew : East (E) or West (W)
  void from_lon(Type dd, Type mm, Type ss, char ew) {
    if (ew == 'E')
      c_rad = M_PI / 180.0 * (dd + mm / 60.0 + ss / 60.0 / 60.0);
    else if (ew == 'W')
      c_rad = -M_PI / 180.0 * (dd + mm / 60.0 + ss / 60.0 / 60.0);
    else {
      if (dd >= 0)
        c_rad = M_PI / 180.0 * (dd + mm / 60.0 + ss / 60.0 / 60.0);
      else
        c_rad = M_PI / 180.0 * (dd - mm / 60.0 - ss / 60.0 / 60.0);
    }
  };
  /// @brief Return in radian
  /// @return Type rad() : radian
  Type rad() const { return c_rad; };
  /// @brief Return in degree
  /// @return Type deg() : degree
  Type deg() const { return c_rad * 180.0 / M_PI; };
  /// @brief Return in minute
  /// @return Type min() : minute
  Type min() const { return c_rad * 180.0 * 60.0 / M_PI; };
  /// @brief Return in second
  /// @return Type sec() : second
  Type sec() const { return c_rad * 180.0 * 60.0 * 60.0 / M_PI; };
  /// @brief Return in right ascension
  /// @return Type ra() : right ascension
  RA<Type> ra() const { return RA<Type>(c_rad); };
  /// @brief Return in declination
  /// @return Type dec() : declination
  Dec<Type> dec() const { return Dec<Type>(c_rad); };
  /// @brief Return in longitude
  /// @return Type longitude() : longitude
  Longitude<Type> longitude() const { return Longitude<Type>(c_rad); };
  /// @brief Return in latitude
  /// @return Type latitude() : latitude
  Latitude<Type> latitude() const { return Latitude<Type>(c_rad); };
  /// @brief Add th
  /// @param[in] th : angle to add
  /// @return Angle<Type> operator+(Angle<Type> th) : result
  Angle<Type> operator+(Angle<Type> th) const {
    return Angle<Type>(c_rad + th.rad());
  };
  /// @brief Subtract by th
  /// @param[in] th : angle to subtract
  /// @return Angle<Type> operator-(Angle<Type> th) : result
  Angle<Type> operator-(Angle<Type> th) const {
    return Angle<Type>(c_rad - th.rad());
  };
  /// @brief Accumulate by th
  /// @param[in] th : angle to accumulate
  void operator+=(Angle th) { c_rad += th.rad(); };
  /// @brief Deaccumulate by th
  /// @param[in] th : angle to deaccumulate
  void operator-=(Angle<Type> th) { c_rad -= th.rad(); };
  /// @brief Multiply by a number
  /// @param[in] aa : the number to multiply
  /// @return Angle<Type> operator*(Type aa) : result
  Angle<Type> operator*(Type aa) const { return Angle<Type>(c_rad * aa); };
  /// @brief Divide by a number
  /// @param[in] aa : the number to divide
  /// @return Angle<Type> operator/(Type aa) : result
  Angle<Type> operator/(Type aa) const { return Angle<Type>(c_rad / aa); };
  /// @brief Divide by an angle
  /// @param[in] th : the angle to divide
  /// @return Type operator/(Angle<Type> th) : result
  Type operator/(Angle<Type> th) const { return c_rad / th.rad(); };
  /// @brief Return if this angle is less than the angle th
  /// @param[in] th : the angle to compare
  /// @return bool operator<(Angle<Type> th) : result (true or false)
  bool operator<(Angle<Type> th) const { return c_rad < th.rad(); };
  /// @brief Return if this angle is less than or equal to the angle th
  /// @param[in] th : the angle to compare
  /// @return bool operator<=(Angle<Type> th) : result (true or false)
  bool operator<=(Angle<Type> th) const { return c_rad <= th.rad(); };
  /// @brief Return if this angle is greater than the angle th
  /// @param[in] th : the angle to compare
  /// @return bool operator>(Angle<Type> th) : result (true or false)
  bool operator>(Angle<Type> th) const { return c_rad > th.rad(); };
  /// @brief Return if this angle is greater than or equal to the angle th
  /// @param[in] th : the angle to compare
  /// @return bool operator>=(Angle<Type> th) : result (true or false)
  bool operator>=(Angle<Type> th) const { return c_rad >= th.rad(); };
  /// @brief Convert to float (32-bit floating point)
  /// @return Angle<float> to_f() : converted result to float
  Angle<float> to_f() const { return Angle<float>(static_cast<float>(c_rad)); }
  /// @brief Convert to double (64-bit floating point)
  /// @return Angle<double> to_d() : converted result to double
  Angle<double> to_d() const {
    return Angle<double>(static_cast<double>(c_rad));
  }
  /// @brief Convert to long double (128-bit floating point)
  /// @return Angle<long double> to_l() : converted result to long double
  Angle<long double> to_l() const {
    return Angle<long double>(static_cast<long double>(c_rad));
  }
};

/**
 * template class RA
 * @brief angle in Right Ascension
 *
 * Type : float,  double,  long double
 * Aliases : fRA (RA<float),
 *           dRA (RA<double),
 *           lRA (RA<long double)
 *
 */
template <class Type> class RA : public Angle<Type> {
private:
public:
  /// @brief Blank constructor, RA = 00:00:00.0
  RA() { this->c_rad = 0.0; };
  /// @brief Constructor in radian
  /// @param[in] theta : RA in Angle
  RA(Angle<Type> theta) : Angle<Type>(theta.rad()) {};
  /// @brief Constructor in Right Ascension
  /// @param[in] hh : hour in Right Ascension
  /// @param[in] mm : minute in Right Ascension
  /// @param[in] ss : second in Right Ascension
  RA(Type hh, Type mm, Type ss)
      : Angle<Type>(M_PI / 12.0 * (hh + mm / 60.0 + ss / 60.0 / 60.0)) {};
  /// @brief Returns hh
  /// @return Type hh() : hour
  Type hh() const { return std::floor(this->c_rad / M_PI * 12.0); };
  /// @brief Returns mm
  /// @return Type mm() : minutes
  Type mm() const {
    return static_cast<Type>(
        static_cast<int>(this->c_rad / M_PI * 12.0 * 60.0) % 60);
  };
  /// @brief Returns ss
  /// @return Type ss() : second
  Type ss() const {
    return 12.0 * 60.0 * 60.0 / M_PI *
           (this->c_rad - M_PI / 12.0 * this->hh() -
            M_PI / 12.0 / 60.0 * this->mm());
  };
  /// @brief Returns RA by string
  /// @return std::string ra() : right ascension
  std::string ra() const {
    std::ostringstream ss;
    ss << this->hh();
    ss << " " << this->mm();
    ss << " " << this->ss();
    return ss.str();
  };
  /// @brief Returns RA<float>
  /// @brief RA<float> to_f() : RA<float>
  RA<float> to_f() const {
    return RA<float>(Angle<float>(static_cast<float>(this->c_rad)));
  };
  /// @brief Returns RA<double>
  /// @return RA<double> to_d() : RA(double)
  RA<double> to_d() const {
    return RA<double>(Angle<double>(static_cast<double>(this->c_rad)));
  };
  /// @brief Returns RA<long double>
  /// @return RA<long double> : RA<long double>
  RA<long double> to_l() const {
    return RA<long double>(
        Angle<long double>(static_cast<long double>(this->c_rad)));
  };
};

/**
 * template class Dec
 * @brief angle in Declination
 *
 * Type : float,  double,  long double
 * Aliases : fDec (Dec<float),
 *           dDec (Dec<double),
 *           lDec (Dec<long double)
 *
 */
template <class Type> class Dec : public Angle<Type> {
private:
public:
  /// @brief Blank constructor, Dec = 00:00:00.0
  Dec() { this->c_rad = 0.0; };
  /// @brief Constructor in Angle
  /// @param[in] theta : Dec in Angle
  Dec(Angle<Type> theta) : Angle<Type>(theta.rad()) {};
  /// @brief Constructor in Declination
  /// @param[in] dd : degree in Dec
  /// @param[in] mm : minutes in Dec
  /// @param[in] ss : second in Dec
  Dec(Type dd, Type mm, Type ss) : Angle<Type>() {
    if (dd >= 0.0)
      this->c_rad = M_PI / 180.0 * (dd + mm / 60.0 + ss / 60.0 / 60.0);
    else
      this->c_rad = M_PI / 180.0 * (dd - mm / 60.0 - ss / 60.0 / 60.0);
  };
  /// @brief Returns dd
  /// @return Type dd() : degree
  Type dd() const {
    if (this->c_rad >= 0.0)
      return std::floor(this->c_rad / M_PI * 180.0);
    return std::ceil(this->c_rad / M_PI * 180.0);
  };
  /// @brief Returns mm
  /// @return Type mm() : minutes
  Type mm() const {
    if (this->c_rad >= 0.0)
      return std::floor((Angle<Type>::deg() - this->dd()) * 60.0);
    return std::abs(std::ceil((Angle<Type>::deg() - this->dd()) * 60.0));
  };
  /// @brief Returns ss
  /// @return Type ss() : seconds
  Type ss() const {
    if (this->c_rad >= 0.0)
      return ((Angle<Type>::deg() - this->dd()) * 60.0 - this->mm()) * 60.0;
    return std::abs(((Angle<Type>::deg() - this->dd()) * 60.0 + this->mm()) *
                    60.0);
  };
  /// @brief Returns declination by string
  /// @return std::string dec() : declination
  std::string dec() const {
    std::ostringstream ss;
    ss << this->dd();
    ss << " " << this->mm();
    ss << " " << this->ss();
    return ss.str();
  };
  /// @brief Returns Dec<float>
  /// @return Dec<float> to_f() : Dec<float>
  Dec<float> to_f() const {
    return Dec<float>(Angle<float>(static_cast<float>(this->c_rad)));
  };
  /// @brief Returns Dec<double>
  /// @return Dec<double> to_d() : Dec<double>
  Dec<double> to_d() const {
    return Dec<double>(Angle<double>(static_cast<double>(this->c_rad)));
  }
  /// @brief Returns Dec<long double>
  /// @return Dec<long double> to_l() : Dec<long double>
  Dec<long double> to_l() const {
    return Dec<long double>(
        Angle<long double>(static_cast<long double>(this->c_rad)));
  }
};

/**
 * template class Longitude
 * @brief angle in Longitude
 *
 * Type : float,  double,  long double
 * Aliases : fLongitude (Longitude<float),
 *           dLongitude (Longitude<double),
 *           lLongitude (Longitude<long double)
 *
 */
template <class Type> class Longitude : public Angle<Type> {
private:
public:
  /// @brief Blank constructor, Longitude = 00:00:00.0 E
  Longitude() { this->c_rad = 0.0; };
  /// @brief Constructor in Angle
  /// @param[in] lon : Longitude in Angle
  Longitude(Angle<Type> lon) : Angle<Type>() { this->c_rad = lon.rad(); };
  /// @brief Constructor in Longitude : E if dd >= 0.0, W if dd< 0.0
  /// @param[in] dd : degree in Longitude
  /// @param[in] mm : minutes in Longitude
  /// @param[in] ss : second in Longitude
  Longitude(Type dd, Type mm, Type ss) : Angle<Type>() {
    if (dd >= 0.0)
      Angle<Type>::from_deg(dd + mm / 60.0 + ss / 60.0 / 60.0);
    else
      Angle<Type>::from_deg(dd - mm / 60.0 - ss / 60.0 / 60.0);
  };
  /// @brief Constructor in Longitude
  /// @param[in] dd : degree in Longitude
  /// @param[in] mm : minutes in Longitude
  /// @param[in] ss : second in Longitude
  /// @param[in] ew : East (E) or West (W) in Longitude
  Longitude(Type dd, Type mm, Type ss, char ew) {
    Angle<Type>::from_lon(dd, mm, ss, ew);
    if (ew == 'E')
      this->from_deg(dd + mm / 60.0 + ss / 60.0);
    else if (ew == 'W')
      this->from_deg(-dd - mm / 60.0 - ss / 60.0 / 60.0);
    else {
      if (dd >= 0.0)
        this->from_deg(dd + mm / 60.0 + ss / 60.0 / 60.0);
      else
        this->from_deg(dd - mm / 60.0 - ss / 60.0 / 60.0);
    }
  };
  /// @brief Returns dd
  /// @return Type dd() : degrees
  Type dd() const {
    if (Angle<Type>::c_rad >= 0.0)
      return std::floor(this->deg());
    return std::ceil(this->deg());
  };
  /// @brief Returns mm
  /// @return Type mm() : minutes
  Type mm() const {
    return std::floor(std::abs(this->deg() - this->dd()) * 60.0);
  };
  /// @brief Returns ss
  /// @return Type ss() : seconds
  Type ss() const {
    return (std::abs(this->deg() - this->dd()) * 60.0 - this->mm()) * 60.0;
  };
  /// @brief Returns ew
  /// @return char ew() : 'E' for East, 'W' for West
  char ew() const {
    if (Angle<Type>::c_rad >= 0.0)
      return 'E';
    return 'W';
  };
  /// @brief Returns longitude by string
  /// @return std::string longitude() : longitude
  std::string longitude() const {
    std::ostringstream ss;
    if (this->dd() >= 0.0)
      ss << this->dd();
    else
      ss << -this->dd();
    ss << " " << this->mm();
    ss << " " << this->ss();
    ss << " " << this->ew();
    return ss.str();
  };
  /// @brief Returns Longitude<float> (32-bit floating point)
  /// @return Longitude<float> to_f() : converted Longitude to float
  Longitude<float> to_f() const {
    return Longitude<float>(
        Angle<float>(static_cast<float>(Angle<Type>::c_rad)));
  };
  /// @brief Returns Longitude<double>
  /// @return Longitude<double> to_d() : converted Longitude to double
  Longitude<double> to_d() const {
    return Longitude<double>(
        Angle<double>(static_cast<double>(Angle<Type>::c_rad)));
  };
  /// @brief Returns Longitude<long double>
  /// @return Longitude<long double> to_l() : converted Longitude to long double
  Longitude<long double> to_l() const {
    return Longitude<long double>(
        Angle<long double>(static_cast<long double>(Angle<Type>::c_rad)));
  };
};

/**
 * template class Latitude
 * @brief angle in Latitude
 *
 * Type : float,  double,  long double
 * Aliases : fLatitude (Latitude<float),
 *           dLatitude (Latitude<double),
 *           lLatitude (Latitude<long double)
 *
 */
template <class Type> class Latitude : public Angle<Type> {
private:
public:
  /// @brief Blank constructor, Latitude = 00:00:00.0 N
  Latitude() : Angle<Type>() { Angle<Type>::c_rad = 0.0; };
  /// @brief Constructor in Angle
  /// @param[in] lat : Latitude in Angle
  Latitude(Angle<Type> lat) : Angle<Type>() { Angle<Type>::c_rad = lat.rad(); };
  /// @brief Constructor in Longitude : N if dd >= 0.0, S if dd < 0.0
  /// @param[in] dd : degree in Latitude
  /// @param[in] mm : minutes in Latitude
  /// @param[in] ss : second in Latitude
  Latitude(Type dd, Type mm, Type ss) : Angle<Type>() {
    if (dd >= 0.0)
      Angle<Type>::from_deg(dd + mm / 60.0 + ss / 60.0 / 60.0);
    else
      Angle<Type>::from_deg(dd - mm / 60.0 - ss / 60.0 / 60.0);
  };
  /// @brief Constructor in Latitude
  /// @param[in] dd : degree in Latitude
  /// @param[in] mm : minutes in Latitude
  /// @param[in] ss : second in Latitude
  /// @param[in] ns : North (N) or South (S) in Latitude
  Latitude(Type dd, Type mm, Type ss, char ns) : Angle<Type>() {
    if (ns == 'N')
      Angle<Type>::from_deg(dd + mm / 60.0 + ss / 60.0 / 60.0);
    else if (ns == 'S')
      Angle<Type>::from_deg(-(dd + mm / 60.0 + ss / 60.0 / 60.0));
    else if (dd >= 0.0)
      Angle<Type>::from_deg(dd + mm / 60.0 + ss / 60.0 / 60.0);
    else
      Angle<Type>::from_deg(dd - mm / 60.0 - ss / 60.0 / 60.0);
  }
  /// @brief Return : Type dd
  /// @return Type dd() : degrees
  Type dd() const {
    if (Angle<Type>::c_rad >= 0.0)
      return std::floor(Angle<Type>::deg());
    return -std::ceil(Angle<Type>::deg());
  };
  /// @brief Return : Type mm
  /// @return Type mm() : minutes
  Type mm() const {
    if (Angle<Type>::c_rad >= 0.0)
      return std::floor(std::abs(Angle<Type>::deg() - this->dd()) * 60.0);
    return std::floor(std::abs(Angle<Type>::deg() + this->dd()) * 60.0);
  };
  /// @brief Return : Type ss
  /// @return ss() : seconds
  Type ss() const {
    if (Angle<Type>::c_rad >= 0.0)
      return std::abs(Angle<Type>::deg() - this->dd()) * 60.0 * 60.0 -
             this->mm() * 60.0;
    return std::abs(Angle<Type>::deg() + this->dd()) * 60.0 * 60.0 -
           this->mm() * 60.0;
  };
  /// @brief Return : char ns
  char ns() const {
    if (Angle<Type>::c_rad >= 0.0)
      return 'N';
    return 'S';
  };
  /// @brief Return : string latitude
  /// @return std::string latitude() : latitude
  std::string latitude() const {
    std::ostringstream ss;
    if (this->dd() >= 0.0)
      ss << this->dd();
    else
      ss << -this->dd();
    ss << " " << this->mm();
    ss << " " << this->ss();
    ss << " " << this->ns();
    return ss.str();
  };
  /// @brief Return : Latitude<float>
  /// @return Latitude<float> to_f() : Latitude<float>
  Latitude<float> to_f() const {
    return Latitude<float>(
        Angle<float>(static_cast<float>(Angle<Type>::c_rad)));
  };
  /// @brief Return : Latitude<double>
  /// @return Latitude<double> to_d() : Latitude<double>
  Latitude<double> to_d() const {
    return Latitude<double>(
        Angle<double>(static_cast<double>(Angle<Type>::c_rad)));
  };
  /// @brief Return : Latitude<long double>
  /// @return Latitude<long double> to_l() : Latitude<long double>
  Latitude<long double> to_l() const {
    return Latitude<long double>(
        Angle<long double>(static_cast<long double>(Angle<Type>::c_rad)));
  };
};

/**
 * template class V2d
 * @brief two dimensional vector
 *
 * Type : float,  double,  long double
 * Aliases : fV2d (V2d<float),
 *           dV2d (V2d<double),
 *           lV2d (V2d<long double)
 *
 */
template <class Type> class V2d {
private:
  /// @brief x component
  Type c_x;
  /// @brief y component
  Type c_y;

public:
  /// @brief Blank constructor, x = 0.0, y = 0.0
  V2d() : c_x(0.0), c_y(0.0) {};
  /// @brief Constructor specifying x and y
  /// @param[in] x : x component
  /// @param[in] y : y component
  V2d(Type x, Type y) : c_x(x), c_y(y) {};
  /// brief Returns x component
  /// Type x() : x component
  Type x() const;
  /// @brief Returns y component
  /// @return Type y() : y component
  Type y() const;
  /// Returns norm
  /// @return Type norm() : norm
  Type norm() const;
  /// @brief Returns square of norm
  /// @return Type n2() : square of norm
  Type n2() const;
  /// @brief Returns unit vector with the same direction
  /// @return V2d<Type> unit() : unit vector
  V2d<Type> unit() const;
  /// @brief Returns this point is on the same line of v1 and v2
  /// @param[in] v1 : a vector
  /// @raram[in] v2 : another vector
  /// @return bool online() : true if on line, false if not
  bool online(V2d<Type> v1, V2d<Type> v2);
  /// @brief Returns whether this point is between v1 and v2 or not
  /// @param[in] v1 : a vector
  /// @raram[in] v2 : another vector
  /// @return bool between() : true if this point is between v1 and v2, false if
  /// not
  bool between(V2d<Type> v1, V2d<Type> v2);
  /// @brief Returns added vector of this vector and vv
  /// @param[in] vv : a vector
  /// @return V2d<Type> : added vector
  V2d<Type> operator+(V2d<Type> vv) const;
  /// @brief Returns subtracted vector by vv
  /// @param[in] vv : a vector
  /// @return V2d<Type> : subtracted vector by vv
  V2d<Type> operator-(V2d<Type> vv) const;
  /// @brief Returns inner product with vv
  /// @param[in] vv : a vector
  /// @return Type : inner product with vv
  Type operator*(V2d<Type> vv) const;
  /// @brief Returns the vector multiplied by a number dd
  /// @param[in] dd : a number
  /// @return V2d<Type> : the vector multiplied by dd
  V2d<Type> operator*(Type dd) const;
  /// @brief Returns the vector multiplied by a 2d matrix mm
  /// @param mm : a 2d matrix
  /// @return V2d<Type> : the vector multiplied by 2d matrix mm
  V2d<Type> operator*(M2d<Type> mm) const;
  /// @brief Returns a vector divided by a number dd
  /// @param[in] dd : a number
  /// @return V2d<Type> : the vector divided by a number dd
  V2d<Type> operator/(Type dd) const;
  /// @brief Accumulates a vector vv
  /// @param[in] vv : a vector
  void operator+=(V2d<Type> vv);
  /// @brief deaccumulates a vector vv
  /// @param[in] vv : a vector
  void operator-=(V2d<Type> vv);
  /// @brief returns this vector is equal to the vector vv or not
  /// @param[in] vv : a vector
  /// @return bool : true if vv and this vector is equal, false if not
  bool operator==(V2d<Type> vv) const;
  /// @brief returns rotated vector by the angle ph
  /// @param[in] ph : rotation angle in radian
  /// @return V2d<Type> : rotated vector
  V2d<Type> operator^(Type ph) const;
  /// @brief returns rotated vector by the angle ph
  /// @param[in] ph : rotation angle in Angle<Type>
  /// @return V2d<Type> : rotated vector
  V2d<Type> operator^(Angle<Type> ph) const { return *this ^ ph.rad(); };
  /// @brief Returns z component of the vector product with vv
  /// @return Type : z component of the vector product
  Type operator%(V2d<Type> vv) const;
  /// @brief Returns converted vector to float
  /// @return V2d<float> to_f() : V2d<float>
  V2d<float> to_f() const;
  /// @brief Returns converted vector to double
  /// @return V2d<double> to_d() : V2d<double>
  V2d<double> to_d() const;
  /// @brief Returns converted vector to long double
  /// @return V2d<long double> to_l() : V2d<long double>
  V2d<long double> to_l() const;
  /// @brief Print parameters
  void print() const;
};

/// @brief Global const V2d_0<Type> : all parameters are zero
template <typename Type> const V2d<Type> V2d_0 = V2d<Type>(0.0, 0.0);

/**
 * template class Seg2d
 * @brief two dimensional segment
 *
 * Type : float,  double,  long double
 * Aliases : fSeg2d (Seg2d<float),
 *           dSeg2d (Seg2d<double),
 *           lSeg2d (Seg2d<long double)
 *
 */
template <class Type> class Seg2d {
protected:
  /// @brief Two end points (protected)
  std::array<V2d<Type>, 2> c_ends;

public:
  /// @brief Blank constructor, (0.0, 0.0) and (0.0, 0.0)
  Seg2d() : c_ends({V2d_0<Type>, V2d_0<Type>}) {};
  /// @brief Comstructor to initialize
  /// @param[in] ends  : teo ends
  Seg2d(std::array<V2d<Type>, 2> ends) : c_ends(ends) {};
  /// @brief Returns ends
  /// @return std::array<V2d<Type>, 2> ends() : ends
  std::array<V2d<Type>, 2> ends() const;
  /// @brief returns wheter point p is on the segment or not
  /// @param[in] p : a point
  /// @return bool on() : true if p is on the segment, false if not
  bool on(V2d<Type> p) const;
  /// @brief Returns the length of the segment
  /// @return Type length() : length
  Type length() const;
  /// method Type Seg2d<Type>::alpha(V2d<Type> p)
  /// @brief returns the length between end[0] and foot point of point p
  /// @param[in] p : a point on the plane
  /// @return Type Seg2d<Type>::alpha(V2d<Type> p) : the length between the
  /// end[0] and foot point of point p
  Type alpha(V2d<Type> p) const;
  /// @brief returns the foot point of the point p
  /// @param[in] p : a point on the plane
  /// @return V2d<Type> Seg2d<Type>::foot() : the foot point
  V2d<Type> foot(V2d<Type> p) const;
  /// @brief returns the closest point on the segment to the point p
  /// @param[in] p : a point on the plane
  /// @return V2d<Type> Seg2d<Type>::foot() : the foot point
  V2d<Type> closest(V2d<Type> p) const;
  /// @brief returns the closest point on the segment to the point p
  /// @param[in] p : a point on the plane
  /// @return Type Seg2d<Type>::distance() : the distance of the point p
  Type distance(V2d<Type> p) const;
  /// @brief returns the opening angle of the segment seen from the point p
  /// @param[in] p : a point on the plane
  /// @return Type Seg2d<Type>::angle() : the opening angle
  Type angle(V2d<Type> p) const;
  /// @brief returns converted segment to Seg2d<float> (32-bit floating point)
  /// @return Seg2d<float> : converted segment to Seg2d<float>
  Seg2d<float> to_f() const;
  /// @brief returns converted segment to Seg2d<double> (64-bit floating point)
  /// @return Seg2d<double> : converted segment to Seg2d<double>
  Seg2d<double> to_d() const;
  /// @brief returns converted segment to Seg2d<long double> (128-bit floating
  /// point)
  /// @return Seg2d<long double> : converted segment to Seg2d<long double>
  Seg2d<long double> to_l() const;
  ///  @brief prints end points
  void print();
};

/**
 * template class M2d
 * Type : float,  double,  long double
 * @brief 2x2 matrix class
 *
 */
template <class Type> class M2d {
protected:
  /// @brief xx component
  Type c_xx = 0.0; // xx component
  /// @brief xy component
  Type c_xy = 0.0; // xy component
  /// @brief yx component
  Type c_yx = 0.0; // yx component
  /// @brief yy component
  Type c_yy = 0.0; // yy component
public:
  /// @brief Blank constructor, xx = xy = yx = yy = 0.0
  M2d() : c_xx(0.0), c_xy(0.0), c_yx(0.0), c_yy(0.0) {};
  /// @brief Constructor to initialize components
  /// @param[in] xx : xx component
  /// @param[in] xy : xy component
  /// @param[in] yx : yx component
  /// @param[in] yy : yy component
  M2d(Type xx, Type xy, Type yx, Type yy)
      : c_xx(xx), c_xy(xy), c_yx(yx), c_yy(yy) {};
  /// @brief Returns xx component
  /// @return Type : xx component
  Type xx() const;
  /// @brief Returns xy component
  /// @return Type : xy component
  Type xy() const;
  /// @brief Returns yx component
  /// @return Type : yx component
  Type yx() const;
  /// @brief Returns yy component
  /// @return Type : yy component
  Type yy() const;
  /// @brief Returns the matrix added by mm
  /// @param[in] mm : a matrix to add
  /// @return M2d<Type> : matrix added by mm
  M2d<Type> operator+(M2d<Type> mm) const;
  /// @brief Returns the matrix subtracted by mm
  /// @param[in] mm : a matrix to subtract
  /// @return M2d<Type> : matrix subtracted by mm
  M2d<Type> operator-(M2d<Type> mm) const;
  /// @brief Returns the matrix mutiplied by mm
  /// @param[in] mm : a matrix to multiply
  /// @return M2d<Type> : matrix multiplied by mm
  M2d<Type> operator*(M2d<Type> mm) const;
  /// @brief Returns the mutiplied matrix by a vector vv
  /// @param[in] vv : a vector to multiply
  /// @return V2d<Type> : multiplied by vv
  V2d<Type> operator*(V2d<Type> vv) const;
  /// @brief Returns the matrix mutiplied by a number ss
  /// param[in] ss : a number to multiply
  /// @return V2d<Type> : the matrix multiplied by ss
  M2d<Type> operator*(Type ss) const;
  /// @brief Returns the transposed matrix
  /// @return M2d<Type> : the transposed matrix
  M2d<Type> operator!() const;
  /// @brief Accumulates by the matrix mm
  /// @param[in] mm : a matrix to add
  void operator+=(M2d<Type> mm);
  /// @brief Returns the determinant of the matrix
  /// @return Type M2d<Type>::det() : the determinant of the matrix
  Type det() const;
  /// @brief Returns the inverse matrix
  /// @return M2d<Type> M2d<Type>::inv() : the inverse matrix
  M2d<Type> inv() const;
  /// @brief Returns eigen values
  /// @return std::array<Type, 2> lambda() : eigenvalues
  std::array<Type, 2> lambda() const;
  /// @brief Returns converted matrix to M2d<float> (32-bit floating point)
  /// @return M2d<float> to_f() : converted matrix to M2d<float>
  M2d<float> to_f() const;
  /// @brief returns converted matrix to M2d<double> (64-bit floating point)
  /// @return M2d<double> to_d() : converted matrix to M2d<double>
  M2d<double> to_d() const;
  /// @brief returns converted matrix to M2d<long double> (128-bit floating
  /// point)
  /// @return M2d<double> to_l() : converted matrix to M2d<long double>
  M2d<long double> to_l() const;
  /// @brief prints matrix
  void print() const; // Print
};

/// @brief Global const M2d_0<Type> : all parameters are zero
template <typename Type> const M2d<Type> M2d_0 = M2d<Type>(0.0, 0.0, 0.0, 0.0);

template <class Type> class Circle;  // Blank definition
template <class Type> class Ellipse; // Blank definition

/**
 * template class Triangle
 * Type : float,  double,  long double
 * @brief 2d Triangle class
 *
 */
template <class Type> class Triangle {
protected:
  /// @brief three apexes (protected)
  std::array<V2d<Type>, 3> c_apex;

public:
  /// @brief Blank constructor, (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)
  Triangle() : c_apex({V2d_0<Type>, V2d_0<Type>, V2d_0<Type>}) {};
  /// @brief Constructor specifying three apexes
  /// @param[in] apex : three apexes
  Triangle(std::array<V2d<Type>, 3> apex) : c_apex(apex) {};
  /// @brief Returns center
  /// @return center
  V2d<Type> cent() const;
  /// @brief Returns apexes
  /// @return Apeexes
  std::array<V2d<Type>, 3> apex() const;
  /// @brief Returns smallest x and smallest y
  /// @return Smallest x and smallest y
  V2d<Type> min() const;
  /// @brief Returns greatest x and greatest y
  /// @return Greatest x and greatest y
  V2d<Type> max() const;
  /// @Returns sides
  /// @return sides
  std::array<Seg2d<Type>, 3> sides() const;
  /// @brief Returns area
  /// @return area
  Type area() const;
  /// @brief Returns closest distance from the point p
  /// @param[in] p : a point on the plane
  /// @return Closest distance
  Type closestd(V2d<Type> p) const;
  /// @brief Returns if the point p is inside or outside
  /// @param[in] p : A point
  /// @return true if p is inside, false if p is not inside
  bool in(V2d<Type> p) const;
  /// @brief Returns if the circle cl is inside or not
  /// @param[in] cl : a circle
  /// @return true if cl is inside, false if cl is not inside
  bool in(Circle<Type> cl) const;
  /// @brief Returns if the circle cl is inside or not with margin
  /// @param[in] cl : a circle
  /// @param[in] fact : margin factor
  /// @return true if cl is inside, false if cl is not inside
  bool in(Circle<Type> cl, Type fact) const;
  /// @brief Returns if the Ellipse el is inside or not
  /// @param[in] el : a ellipse
  /// @return true if el is inside, false if el is not inside
  bool in(Ellipse<Type> el) const;
  /// @brief Returns if the ellipse el is inside or not with margin
  /// @param[in] el : a ellipse
  /// @param[in] fact : margin factor
  /// @return true if el is inside, false if el is not inside
  bool in(Ellipse<Type> el, Type fact) const;
  /// @brief Returns if the triangle is inside of the circle
  /// @param[in] cl : a circle on the plane
  /// @return bool insideof(Circle<Type> cl) : the triangle is inside of
  /// the circle (true : inside, false : not inside)
  bool insideof(Circle<Type> cl) const;
  /// @brief Returns converted triangle to Triangle<float> (32-bit floating
  /// point)
  /// @return Triangle<float> to_f() : converted Triangle to Triangle<float>
  Triangle<float> to_f() const;
  /// @brief returns converted triangle to Triangle<double> (64-bit floating
  /// point)
  /// @return Triangle<double> to_d() : converted Triangle to Triangle<double>
  Triangle<double> to_d() const;
  /// @brief returns converted triangle to Triangle<long double> (128-bit
  /// floating point)
  /// @return Triangle<long double> to_l() : converted Triangle to Triangle<long
  /// double>
  Triangle<long double> to_l() const;
  /// @brief prints triangle
  void print() const;
};

/**
 * template class Rectangle
 * Type : float,  double,  long double
 * @brief Class to express rectangle
 *
 */
template <class Type> class Rectangle {
protected:
  /// @brief lower left corner
  V2d<Type> c_min;
  /// @brief upper right corner
  V2d<Type> c_max;

public:
  /// @brief Blank constructor
  Rectangle() : c_min(V2d_0<Type>), c_max(V2d_0<Type>) {};
  /// @brief Constructor to initialize
  /// @param[in] min : lower left corner
  /// @param[in] max : upper right corner
  Rectangle(V2d<Type> min, V2d<Type> max) : c_min(min), c_max(max) {};
  /// @brief Returns lower left corner
  /// @return V2d<Type> min() : lower left corner
  V2d<Type> min() const;
  /// @brief Returns upper right corner
  /// @return V2d<Type> max() : upper right corner
  V2d<Type> max() const;
  /// @brief Returns center of the rectangle
  /// @return V2d<Type> cent() : center
  V2d<Type> cent() const;
  /// @brief Returns width
  /// @return Type d() : width
  Type d() const;
  /// @brief Returns area
  /// @return Type area() : area
  Type area() const;
  /// brief Returns if point p is in or not
  /// @param[in] p : a point on the plane
  /// @return bool in() : p is in (true : yes, false : no)
  bool in(V2d<Type> p) const;
  /// @brief Returns corners
  /// @return std::array<V2d<Type>, 4> corners() : four corners
  std::array<V2d<Type>, 4> corners() const;
  /// @brief Returns segments
  /// @return std::array<Seg2d<Type>, 4> segs() : four segments
  std::array<Seg2d<Type>, 4> segs() const;
  /// @brief returns converted rectangle to Rectangle<float> (32-bit floating
  /// point)
  /// @return Rectangle<float> to_f() : converted Rectangle to Rectangle<float>
  Rectangle<float> to_f() const;
  /// @brief returns converted rectangle to Rectangle<double> (64-bit floating
  /// point)
  /// @return Rectangle<double> to_d() : converted Rectangle to
  /// Rectangle<double>
  Rectangle<double> to_d() const;
  /// @brief returns converted rectangle to Rectangle<long double> (128-bit
  /// floating point)
  /// @return Rectangle<long double> to_l() : converted Rectangle to
  /// Rectangle<long double>
  Rectangle<long double> to_l() const;
  /// @brief prints rectangle
  void print() const;
};

/// @brief Global const Rectangle_0<Type> : all parameters are zero
template <typename Type>
Rectangle<Type> const Rectangle_0 = Rectangle<Type>(V2d_0<Type>, V2d_0<Type>);

/**
 * template class Square
 * Type : float,  double,  long double
 * @brief Class to express Square
 *
 */
template <class Type> class Square : public Rectangle<Type> {
private:
  /// @brief Center of the suare
  V2d<Type> c_cent;
  /// @brief half width
  Type c_d;

public:
  /// @brief Blank constructor (0.0, 0.0), 0.0
  Square() : c_cent(V2d_0<Type>), c_d(0.0) {};
  /// @brief Constructor for initialize
  /// @param[in] cent : Center
  /// @param[in] d : half width
  Square(V2d<Type> cent, Type d)
      : Rectangle<Type>(cent - V2d<Type>(d, d), cent + V2d<Type>(d, d)),
        c_cent(cent), c_d(d) {};
  /// @brief Returns center of the square
  /// @return V2d<Type> cent() : center of the square
  V2d<Type> cent() const;
  /// @brief Returns the half length of the segments
  /// @return Type d() : the half length of the segments
  Type d() const;
  /// @brief Returns the expanded square by a factor of fact
  /// @return Square<Type> enlarge(Type fact) : the expanded square by a factor
  /// of fact
  Square<Type> enlarge(Type fact) const;
  /// @brief Returns the pair of the splitted triangles
  /// @return std::array<Triangle<Type>, 2> to_tri() : the pair of the splitted
  /// triangles
  std::array<Triangle<Type>, 2> to_tri() const;
};

/**
 * template struct TriInOut
 * @brief struct to express in apexes and out apexes of the triangle
 */
template <typename Type> struct TriInOut {
  /// @brief apexes in the Circle
  std::vector<V2d<Type>> in;
  /// @brief apexes out of the circle
  std::vector<V2d<Type>> out;
};

/**
 * template struct AreaTri
 * @brief Struct to express "in" area, "out" area and "crescent" area to the
 * circle
 */
template <typename Type> struct AreaTri {
  /// @brief area in the circle of the triangle
  Type in;
  /// @brief crescent area between the circle and the triangle
  Type cres;
  /// @brief area out of the circle
  Type out;
};

/**
 * template class Circle
 * Type : float,  double,  long double
 * @brief Class to express Circle
 *
 */
template <class Type> class Circle {
protected:
  /// @brief Center of the circle
  V2d<Type> c_cent = V2d<Type>(0.0, 0.0);
  /// @brief Radious
  Type c_r = 0.0;

public:
  /// @brief Blank constructor (0.0, 0.0), 0.0
  Circle() : c_cent(V2d<Type>(0.0, 0.0)), c_r(0.0) {};
  Circle(Type r) : c_cent(V2d<Type>(0.0, 0.0)), c_r(r) {};
  /// @brief Constructor to define
  /// @param[in] cent : Center
  /// @param[in] r : Radius
  Circle(V2d<Type> cent, Type r) : c_cent(cent), c_r(r) {};
  /// @brief Constructor with cl, and cent
  /// @param[in] cl : circle parameters
  /// @param[in] cent : center of the circle
  Circle(Circle<Type> cl, V2d<Type> cent) : c_cent(cent) { c_r = cl.r(); };
  /// @brief Returns center of the square
  /// @return V2d<Type> cent() : center of the square
  V2d<Type> cent() const;
  /// @brief Returns radious
  /// @return Type r() : radious of the circle
  Type r() const;
  /// @brief Returns distance from the center of the circle to the point p
  /// @param[in] p : a point
  /// @return Type rho() : distance from the center
  Type rho(V2d<Type> p) const;
  /// @brief Returns area of the circle
  /// @return Type area() : Area of the circle
  Type area() const;
  /// @brief Returns whether point p is inside of the circle or not
  /// @param[in] p : a point
  /// @return bool in() : true if p is inside, false if not
  bool in(V2d<Type> p) const;
  /// @brief Returns the ratio a to get the cross point
  /// @param[in] seg : a segment
  /// @return Type cross_a() : Returns the ratio a to get the cross point
  Type cross_a(Seg2d<Type> seg) const;
  /// @brief Returns the cross point of the segment
  /// @param[in] seg : a segment
  /// @return V2d<Type> Circle<Type>::cross() : Cross point
  V2d<Type> cross(Seg2d<Type> seg) const;
  /// @brief Returns "in" and "out" apexes
  /// @param[in] tri : a triangle
  /// @return TriInOut<Type> Circle<Type>::inout() : "in" and "out" apexes
  TriInOut<Type> inout(Triangle<Type> tri) const;
  /// @brief Returns the area of the crescent defined by th
  /// @param[in] th : angle of the sector
  /// @return Type Circle<Type>::area_th() : Area of the crescent
  Type area_th(Type th) const;
  /// @brief Returns the areas of "in", "out", "crescent" regions
  /// @param[in] tri : a triangle
  /// @return AreaTri<Type> Circle<Type>::areatri() : the areas of "in", "out",
  /// "crescent" regions
  AreaTri<Type> areatri(Triangle<Type> tri) const;
  /// @brief Returns the ratio of the "in" area and all area of the triangle
  /// @param[in] tri : a triangle
  /// @return Type Circle<Type>::ratio_ia() : the ratio of the "in" area and all
  /// area of the triangle
  Type ratio_ia(Triangle<Type> tri) const;
  /// @brief Returns the thickness of te swell defined by th
  /// @param[in] th : opening angle of the swell
  /// @return Type Circle<Type>::swell() : the thickness of te swell
  Type swell(Type th) const;
  /// @brief Returns the thickness of te swell defined by the segment
  /// @param[in] seg : a segment
  /// @return Type Circle<Type>::swell() : the thickness of te swell
  Type swell(Seg2d<Type> seg) const;
  /// @brief Returns the thickness of te swell defined by the triangle
  /// @param[in] tri : a triangle
  /// @return Type Circle<Type>::swell() : the thickness of te swell
  Type swell(Triangle<Type> tri) const;
  /// @brief Returns the opening angle of the triangle
  /// @param[in] tri : a triangle
  /// @return Type Circle<Type>::angle() : the opening angle of the triangle
  Type angle(Triangle<Type> tri) const;
  /// @brief Returns the gap between the sides of the triangle and the circl
  /// @param[in] tri : a triangle
  /// @return Type Circle<Type>::gap() :
  Type gap(Triangle<Type> tri) const;
  /// @brief Returns the smallest square that include the circle
  /// @return Square<Type> Circle<Type>::include() : the smallest square that
  /// include the circle
  Square<Type> include() const;
  /// @brief returns the converted circle to Circle<float> (32-bit floating
  /// point)
  /// @return Circle<float> Circle<Type>::to_f() : the converted circle to
  /// Circle<float> (32-bit floating point)
  Circle<float> to_f() const;
  /// @brief returns the converted circle to Circle<double> (64-bit floating
  /// point)
  /// @return Circle<double> Circle<Type>::to_d() : the converted circle to
  /// Circle<double> (64-bit floating point)
  Circle<double> to_d() const;
  /// @brief returns the converted circle to Circle<long double> (128-bit
  /// floating point)
  /// @return Circle<long double> Circle<Type>::to_l() : the converted circle to
  /// Circle<double> (128-bit floating point)
  Circle<long double> to_l() const;
  /// @brief prints circle
  void print() const; // Print
};

/// @brief Global const Circle_0<Type> : all parameters are zero
template <typename Type>
const Circle<Type> Circle_0 = Circle<Type>(V2d_0<Type>, 0.0);

/**
 * template class Ellipse
 * Type : float,  double,  long double
 * @brief Class to express Ellipse
 *
 */
template <class Type> class Ellipse {
protected:
  /// @brief Center of the ellipse
  V2d<Type> c_cent;
  /// @brief Semi-major and semi-minor axes
  V2d<Type> c_ab;
  /// @brief Inclination
  Type c_ph;

public:
  /// @brief Blank constructor
  Ellipse() : c_cent(V2d_0<Type>), c_ab(V2d_0<Type>), c_ph(0.0) {}
  /// @brief Constructor for no inclination
  Ellipse(V2d<Type> cent, Type a, Type b)
      : c_cent(cent), c_ab(V2d<Type>(a, b)) {};
  /// @brief Constructor from ell
  Ellipse(const Ellipse<Type> &ell)
      : c_cent(ell.cent()), c_ab(ell.ab()), c_ph(ell.ph()) {};
  /// @brief Constructor with inclination
  Ellipse(V2d<Type> cent, Type a, Type b, Type ph)
      : c_cent(cent), c_ab(V2d<Type>(a, b)), c_ph(ph) {};
  /// @brief Constructor with shift center
  Ellipse(Ellipse<Type> el, V2d<Type> cent)
      : c_cent(cent), c_ab(el.ab()), c_ph(el.ph()) {};
  /// @brief Constructor with shift center and rotate
  Ellipse(Ellipse<Type> el, V2d<Type> cent, Type ph)
      : c_cent(cent), c_ab(el.ab()), c_ph(ph) {};
  /// @brief Returns the center of the ellipse
  /// @return V2d<Type> Ellipse<Type>::cent() : Center
  V2d<Type> cent() const;
  /// @brief Returns the Semi-major and semi-minor axes
  /// @return V2d<Type> Ellipse<Type>::ab() : Semi-major and semi-minor axes
  V2d<Type> ab() const;
  /// @brief Returns the inclination
  /// @return Type Ellipse<Type>::ph() : Inclination
  Type ph() const;
  /// @brief Returns the projected radious to the direction of the point p
  /// @param[in] p : a point
  /// @return Type Ellipse<Type>::r() : the projected radious to the direction
  /// of the point p
  Type r(V2d<Type> p) const;
  /// @brief Returns the normalize elliptical distance to the point p
  /// @param[in] p : a point
  /// @return Type Ellipse<Type>::rho() : the normalized elliptical distance of
  /// the point p
  Type rho(V2d<Type> p) const;
  /// @brief Returns the area of the ellipse
  // @return Type Ellipse<Type>::area() : the area
  Type area() const;
  /// @brief Returns the normalize elliptical distance to the point p
  /// @param[in] p : a point
  /// @return Type Ellipse<Type>::dist(V2d<Type> p) : the normalized elliptical
  /// distance of the point p
  Type dist(V2d<Type> p) const;
  /// @brief Returns if  the point p is inside or not inside
  /// @param[in] p : a point
  /// @return bool Ellipse<Type>::in() : the point p is inside or not inside
  bool in(V2d<Type> p) const;
  /// @brief Returns the equivalent point of p for the unit circle
  /// @param[in] p : a point
  /// @return V2d<Type> Ellipse<Type>::eqvCirc01() : the equivalent point
  V2d<Type> eqvCirc01(V2d<Type> p) const;
  /// @brief Returns the equivalent segment for the unit circle
  /// @param[in] seg : a segment
  /// @return Seg2d<Type> Ellipse<Type>::eqvCirc01() : the equivalent segment
  Seg2d<Type> eqvCirc01(Seg2d<Type> seg) const;
  /// @brief Returns the equivalent segment facotr (cross_a) for the unit circle
  /// @param[in] seg : a segment
  /// @return Type Ellipse<Type>::cross_a() : the unit circle equivalent segment
  /// factor
  Type cross_a(Seg2d<Type> seg) const;
  /// @brief Returns the cross point of the segment
  /// @param[in] seg : a segment
  /// @return V2d<Type> Ellipse<Type>::cross() : the cross point of the segment
  V2d<Type> cross(Seg2d<Type> seg) const;
  /// @brief Returns the ratio of the "in" area and all area of the triangle
  /// @param[in] tri : a triangle
  /// @return Type Ellipse<Type>::ratio_ia() : the ratio of the "in" area and
  /// all area of the triangle
  Type ratio_ia(Triangle<Type> tri) const;
  /// @brief Returns the smallest square that include the circle
  /// @return Square<Type> Ellipse<Type>::include() : the smallest square that
  /// include the circle
  Square<Type> include() const;
  /// @brief returns the converted ellipse to Ellipse<float> (32-bit floating
  /// point)
  /// @return Ellipse<float> Ellipse<Type>::to_f() : the converted circle to
  /// Ellipse<float> (32-bit floating point)
  Ellipse<float> to_f() const;
  /// @brief returns the converted ellipse to Ellipse<double> (64-bit floating
  /// point)
  /// @return Ellipse<double> Ellipse<Type>::to_d() : the converted ellipse to
  /// Ellipse<double> (64-bit floating point)
  Ellipse<double> to_d() const;
  /// @brief returns the converted ellipse to Ellipse<long double> (128-bit
  /// floating point)
  /// @return Ellipse<long double> Ellipse<Type>::to_l() : the converted Ellipse
  /// to Ellipse<double> (128-bit floating point)
  Ellipse<long double> to_l() const;
  /// @brief prints ellipse
  void print() const; // Print
};

/// @brief Global const Ellipse_0<Type> : all parameters are zero
template <typename Type>
Ellipse<Type> const Ellipse_0 = Ellipse<Type>(V2d_0<Type>, 0.0, 0.0);

/// @brief Global const IncEC<Type> : inclination of Ecliptic
template <typename Type> const Type IncEC = 23.43928 / 180.0 * M_PI;

// Blank class definitions
template <class Type> class M3d;
template <class Type> class Rotate_x;
template <class Type> class Rotate_y;
template <class Type> class Rotate_z;
template <class Type> class DirCos;

/**
 * template class V3d
 * Type : float,  double,  long double
 * @brief Class to express 3d vector
 *
 */
template <class Type> class V3d {
protected:
  /// @brief x component
  Type c_x;
  /// @brief y component
  Type c_y;
  /// @brief z component
  Type c_z;

public:
  /// @brief Blank constructor : all zero
  V3d() : c_x(0.0), c_y(0.0), c_z(0.0) {};
  /// @brief Constructor to define elements
  /// @param[in] x : x-component
  /// @param[in] y : x-component
  /// @param[in] z : x-component
  V3d(Type x, Type y, Type z) : c_x(x), c_y(y), c_z(z) {};
  /// @brief Constructor to define direction cosine
  /// @param[in] phi : longitude
  /// @param[in] theta : latitude
  V3d(Angle<Type> phi, Angle<Type> theta)
      : c_x(std::cos(theta.rad()) * std::cos(phi.rad())),
        c_y(std::cos(theta.rad()) * std::sin(phi.rad())),
        c_z(std::sin(theta.rad())) {};
  /// @brief Return x component
  Type x() const { return c_x; };
  /// @brief Return y component
  Type y() const { return c_y; };
  /// @brief Return z component
  Type z() const { return c_z; };
  /// @brief Return norm of the vector
  Type norm() const { return std::sqrt(c_x * c_x + c_y * c_y + c_z * c_z); };
  /// @brief Return declination of the vector
  Dec<Type> dec() const {
    return Dec<Type>(Angle<Type>(std::asin(c_z / this->norm())));
  };
  /// @brief Return right ascension of the vector
  RA<Type> ra() const {
    Type ac = std::acos(c_x / std::cos(this->dec().rad()) / this->norm());
    if (c_y >= 0.0)
      return RA(Angle<Type>(Angle<Type>(ac)));
    return RA(Angle<Type>(Angle<Type>(2.0 * M_PI - ac)));
  };
  /// @brief Return latitude of the vector
  Latitude<Type> latitude() const {
    return Latitude(Angle<Type>(std::asin(c_z / this->norm())));
  };
  /// @brief Return longitude of the vector
  Longitude<Type> longitude() const {
    Type ac = std::acos(c_x / std::cos(this->latitude().rad()) / this->norm());
    if (c_y >= 0.0)
      return Longitude<Type>(ac);
    return Longitude<Type>(-ac);
  };
  /// @brief Returns unit vector
  V3d<Type> unit() const { return *this / this->norm(); };
  /// @brief Adds a vector vv
  /// @param[in] vv : 3d vector to add
  V3d<Type> operator+(V3d<Type> vv) const {
    return V3d<Type>(c_x + vv.x(), c_y + vv.y(), c_z + vv.z());
  };
  /// @brief Accumulates by a vector vv
  /// @param[in] vv : 3d vector to accumulate
  void operator+=(V3d<Type> vv) {
    c_x += vv.x();
    c_y += vv.y();
    c_z += vv.z();
  };
  /// @brief Subtracts by a vector vv
  /// @param[in] vv : 3d vector to subtract
  V3d<Type> operator-(V3d<Type> vv) const {
    return V3d<Type>(c_x - vv.x(), c_y - vv.y(), c_z - vv.z());
  };
  /// @brief Returns inner product with a vector vv
  /// @param[in] vv : 3d vector to multiply
  Type operator*(V3d<Type> vv) const {
    return c_x * vv.x() + c_y * vv.y() + c_z * vv.z();
  };
  /// @brief Returns the product by the number dd
  /// @param[in] dd : the number to multiply
  V3d<Type> operator*(Type dd) const {
    return V3d<Type>(c_x * dd, c_y * dd, c_z * dd);
  };
  /// @brief Returns the vector product by a vector vv
  /// @param[in] vv : 3d vector to multiply
  V3d<Type> operator%(V3d<Type> vv) const {
    return V3d<Type>(c_y * vv.z() - vv.y() * c_z, -c_x * vv.z() + vv.x() * c_z,
                     c_x * vv.y() - vv.x() * c_y);
  };
  /// @brief Returns if the vector is equal to the vector vv
  /// @param[in] vv : 3d vector to compare
  bool operator==(V3d<Type> vv) const {
    return std::abs(vv.x() - c_x) <= std::numeric_limits<Type>::epsilon() &&
           std::abs(vv.y() - c_y) <= std::numeric_limits<Type>::epsilon() &&
           std::abs(vv.z() - c_z) <= std::numeric_limits<Type>::epsilon();
  }
  /// @brief Returns the vector devided by the number aa
  /// @param[in] aa : the number to divide
  V3d<Type> operator/(Type aa) const {
    return V3d<Type>(c_x / aa, c_y / aa, c_z / aa);
  };
  /// @brief Rotates around x axis
  /// @param[in] phi : the angle to rotate
  V3d<Type> rotate_x(Angle<Type> phi) const {
    return Rotate_x<Type>(phi).m3d() * (*this);
  };
  /// @brief Rotates around y axis
  /// @param[in] phi : the angle to rotate
  V3d<Type> rotate_y(Angle<Type> phi) const {
    return Rotate_y<Type>(phi).m3d() * (*this);
  };
  /// @brief Rotates around z axis
  V3d<Type> rotate_z(Angle<Type> phi) const {
    return Rotate_z<Type>(phi).m3d() * (*this);
  };
  /// @brief Cnverts from Equatorial to Ecliptic
  V3d<Type> EQtoEC() const { return this->rotate_x(Angle<Type>(IncEC<Type>)); };
  /// @brief Cnverts from Ecliptic to Equatorial
  V3d<Type> ECtoEQ() const {
    return this->rotate_x(Angle<Type>(-IncEC<Type>));
  };
  /// @brief Returns V3d<float>
  V3d<float> to_f() const {
    return V3d<float>(static_cast<float>(c_x), static_cast<float>(c_y),
                      static_cast<float>(c_z));
  };
  /// @brief Returns V3d<double>
  V3d<double> to_d() const {
    return V3d<double>(static_cast<double>(c_x), static_cast<double>(c_y),
                       static_cast<double>(c_z));
  };
  /// @brief Returns V3d<long double>
  V3d<long double> to_l() const {
    return V3d<long double>(static_cast<long double>(c_x),
                            static_cast<long double>(c_y),
                            static_cast<long double>(c_z));
  };
  /// @brief Prints components
  void print() const {
    std::cout << "x = " << c_x << ", y = " << c_y << ", z = " << c_z
              << std::endl;
  };
};

/**
 * template class M3d
 * Type : float,  double,  long double
 * @brief Class to express 3d matrix
 *
 */
template <class Type> class M3d {
protected:
  /// @brief x-component (3d vector)
  V3d<Type> c_x;
  /// @brief y-component (3d vector)
  V3d<Type> c_y;
  /// @brief z-component (3d vector)
  V3d<Type> c_z;

public:
  /// @brief Blank constructor
  M3d() : c_x(V3d<Type>()), c_y(V3d<Type>()), c_z(V3d<Type>()) {};
  /// @brief Constructor initialized by M3d
  /// @param[in] &mm : 3d matrix
  M3d(M3d<Type> &mm) : c_x(mm.x()), c_y(mm.y()), c_z(mm.z()) {};
  /// @brief Constructor initialized by V3d's
  /// @param[in] x : 3d vector for x component
  /// @param[in] y : 3d vector for y component
  /// @param[in] z : 3d vector for z component
  M3d(V3d<Type> x, V3d<Type> y, V3d<Type> z) : c_x(x), c_y(y), c_z(z) {};
  /// @brief Constructor initialized by 3x3 numbers
  /// @param[in] arr : array for all component
  M3d<Type>(std::array<std::array<Type, 3>, 3> arr)
      : c_x(V3d<Type>(arr[0][0], arr[0][1], arr[0][2])),
        c_y(V3d<Type>(arr[1][0], arr[1][1], arr[1][2])),
        c_z(V3d<Type>(arr[2][0], arr[2][1], arr[2][2])){};
  /// @brief Return x componet vector
  V3d<Type> x() const { return c_x; };
  /// @brief Return y componet vector
  V3d<Type> y() const { return c_y; };
  /// @brief Return z componet vector
  V3d<Type> z() const { return c_z; };
  /// @brief Adds by a matrix
  /// @param[in] mm : the matrix to add
  /// @return M3d<Type> : added matrix
  M3d<Type> operator+(M3d<Type> mm) const {
    return M3d<Type>(c_x + mm.x(), c_y + mm.y(), c_z + mm.z());
  };
  /// @brief Subtract by a matrix
  /// @param[in] mm : the matrix to subtract
  /// @return M3d<Type> : subtracted matrix
  M3d<Type> operator-(M3d<Type> mm) const {
    return M3d<Type>(c_x - mm.x(), c_y - mm.y(), c_z - mm.z());
  };
  /// @brief Multiply by a matrix
  /// @param[in] vv : the vector to multiply
  /// @return V3d<Type> : multiplied vector
  V3d<Type> operator*(V3d<Type> vv) const {
    return V3d<Type>(c_x * vv, c_y * vv, c_z * vv);
  };
  /// @brief Multiply by a matrix
  /// @param[in] mm : the matrix to multiply
  /// @return M3d<Type> : multiplied matrix
  M3d<Type> operator*(M3d<Type> mm) const {
    return M3d<Type>(mm.x() * c_x.x() + mm.y() * c_y.x() + mm.z() * c_z.x(),
                     mm.x() * c_x.y() + mm.y() * c_y.y() + mm.z() * c_z.y(),
                     mm.x() * c_x.z() + mm.y() * c_y.z() + mm.z() * c_z.z());
  };
  /// @brief Multiply by a number
  /// @param[in] ss : the number to multiply
  /// @return M3d<Type> : multiplied matrix
  M3d<Type> operator*(Type ss) const {
    return M3d<Type>(c_x * ss, c_y * ss, c_z * ss);
  };
  /// @brief Make transpose
  /// @return M3d<Type> operator!() : transposed matrix
  M3d<Type> operator!() const {
    return M3d<Type>(V3d<Type>(c_x.x(), c_y.x(), c_z.x()),
                     V3d<Type>(c_x.y(), c_y.y(), c_z.y()),
                     V3d<Type>(c_x.z(), c_y.z(), c_z.z()));
  };
  /// @brief Accumulate by a matrix
  /// @param[in] mm : the matrix to accumulate
  void operator+=(M3d<Type> mm) {
    c_x += mm.x();
    c_y += mm.y();
    c_z += mm.z();
  };
  /// @brief Return the determinant
  /// @return Type : determinant
  Type det() const {
    return c_x.x() * c_y.y() * c_z.z() + c_y.x() * c_z.y() * c_x.z() +
           c_z.x() * c_x.y() * c_y.z() - c_z.x() * c_y.y() * c_x.z() -
           c_z.y() * c_y.z() * c_x.x() - c_z.z() * c_x.y() * c_y.x();
  };
  /// @brief Return the converted matrix to M3d<float>
  /// @return M3d<float> : determinant in float (32-bit floating points)
  M3d<float> to_f() const {
    return M3d<float>(c_x.to_f(), c_y.to_f(), c_z.to_f());
  };
  /// @brief Return the converted matrix to M3d<double>
  /// @return M3d<double> : determinant in double (64-bit floating points)
  M3d<double> to_d() const {
    return M3d<double>(c_x.to_d(), c_y.to_d(), c_z.to_d());
  };
  /// @brief Return the converted matrix to M3d<long double>
  /// @return M3d<long double> : determinant in long double (128-bit floating
  /// points)
  M3d<long double> to_l() const {
    return M3d<long double>(c_x.to_l(), c_y.to_l(), c_z.to_l());
  };
  /// @brief Prints the matrix
  void print() const {
    std::cout << std::fixed << std::setprecision(8) << "       "
              << "c_x"
              << "         "
              << "c_y"
              << "          "
              << "c_z" << std::endl;
    std::cout << "x  " << c_x.x() << "  " << c_y.x() << "   " << c_z.x()
              << std::endl;
    std::cout << "y  " << c_x.y() << "  " << c_y.y() << "   " << c_z.y()
              << std::endl;
    std::cout << "z  " << c_x.z() << "  " << c_y.z() << "   " << c_z.z()
              << std::endl;
  };
};

/**
 * template class Rotate_x
 * Type : float,  double,  long double
 * @brief Class to express 3d rotation around x-axis
 *
 */
template <class Type> class Rotate_x : public M3d<Type> {
private:
public:
  /// @brief Constructor to specify angle by radian
  /// @param[in] phi :  rotation angle in radian
  Rotate_x(Type phi)
      : M3d<Type>(V3d<Type>(1.0, 0.0, 0.0),
                  V3d<Type>(0.0, std::cos(phi), std::sin(phi)),
                  V3d<Type>(0.0, -std::sin(phi), std::cos(phi))) {};
  /// @brief Constructor to specify angle by Angle
  /// @param[in] phi :  rotation angle by Angle<Type>
  Rotate_x(Angle<Type> phi)
      : M3d<Type>(V3d<Type>(1.0, 0.0, 0.0),
                  V3d<Type>(0.0, std::cos(phi.rad()), std::sin(phi.rad())),
                  V3d<Type>(0.0, -std::sin(phi.rad()), std::cos(phi.rad()))) {};
  /// @brief Return the rotation matrix
  /// @return M3d<Type> : 3d matrix to rotate
  M3d<Type> m3d() const { return M3d<Type>(this->c_x, this->c_y, this->c_z); };
};

/**
 * template class Rotate_y
 * Type : float,  double,  long double
 * @brief Class to express 3d rotation around y-axis
 *
 */
template <class Type> class Rotate_y : public M3d<Type> {
private:
public:
  /// @brief Constructor to specify angle by radian
  /// @param[in] phi :  rotation angle in radian
  Rotate_y(Type phi)
      : M3d<Type>(V3d<Type>(std::cos(phi), 0.0, -std::sin(phi)),
                  V3d<Type>(0.0, 1.0, 0.0),
                  V3d<Type>(std::sin(phi), 0.0, std::cos(phi))) {};
  /// @brief Constructor to specify angle by Angle
  /// @param[in] phi :  rotation angle by Angle<Type>
  Rotate_y(Angle<Type> phi)
      : M3d<Type>(V3d<Type>(std::cos(phi.rad()), 0.0, -std::sin(phi.rad())),
                  V3d<Type>(0.0, 1.0, 0.0),
                  V3d<Type>(std::sin(phi.rad()), 0.0, std::cos(phi.rad()))) {};
  /// @brief Return the rotation matrix
  /// @return M3d<Type> : 3d matrix to rotate
  M3d<Type> m3d() const { return M3d<Type>(this->c_x, this->c_y, this->c_z); };
};

/**
 * template class Rotate_z
 * Type : float,  double,  long double
 * @brief Class to express 3d rotation around z-axis
 *
 */
template <class Type> class Rotate_z : public M3d<Type> {
private:
public:
  /// @brief Constructor to specify angle by radian
  /// @param[in] phi :  rotation angle in radian
  Rotate_z(Type phi)
      : M3d<Type>(V3d<Type>(std::cos(phi), std::sin(phi), 0.0),
                  V3d<Type>(-std::sin(phi), std::cos(phi), 0.0),
                  V3d<Type>(0.0, 0.0, 1.0)) {};
  /// @brief Constructor to specify angle by Angle
  /// @param[in] phi :  rotation angle by Angle<Type>
  Rotate_z(Angle<Type> phi)
      : M3d<Type>(V3d<Type>(std::cos(phi.rad()), std::sin(phi.rad()), 0.0),
                  V3d<Type>(-std::sin(phi.rad()), std::cos(phi.rad()), 0.0),
                  V3d<Type>(0.0, 0.0, 1.0)) {};
  /// @brief Return the rotation matrix
  /// @return M3d<Type> : 3d matrix to rotate
  M3d<Type> m3d() const { return M3d<Type>(this->c_x, this->c_y, this->c_z); };
};

/**
 * template class Align_x
 * Type : float,  double,  long double
 * @brief Class to express 3d rotation to align x-axis to the direction of the
 * vector v
 *
 */
template <class Type> class Align_x {
private:
  /// @brief Rotation around z
  Rotate_z<Type> c_rtz;
  /// @brief Rotation around y
  Rotate_y<Type> c_rty;

public:
  /// @brief Constructor to specify 3d vector v
  /// @param[in] v : the 3d vector
  Align_x(V3d<Type> v)
      : c_rtz(Rotate_z<Type>(v.ra())), c_rty(Rotate_y<Type>(v.dec())) {};
  /// @brief Return 3d matrix to rotate around z axis
  /// @return M3d<Type> rtz() : 3d matrix to rotate around z axis
  M3d<Type> rtz() const { return c_rtz.m3d(); };
  /// @brief Return 3d matrix to rotate around y axis
  /// @return M3d<Type> rty() : 3d matrix to rotate around y axis
  M3d<Type> rty() const { return c_rty.m3d(); };
  /// @brief Return 3d matrix to align direction to v
  /// @return M3d<Type> m3d() : 3d matrix to align direction to v
  M3d<Type> m3d() const { return c_rty.m3d() * c_rtz.m3d(); };
  /// @brief Return converted 3d vector to the new coordinats
  /// @param[in] vv : a 3d vector
  /// @return vv in the new coordinate
  V3d<Type> align(V3d<Type> vv) const { return this->m3d() * vv; };
};

/**
 * template class Align_y
 * Type : float,  double,  long double
 * @brief Class to express 3d rotation to align y-axis to the direction of the
 * vector v
 *
 */
template <class Type> class Align_y {
private:
  /// @brief Rotation around z
  Rotate_z<Type> c_rtz;
  /// @brief Rotation around x
  Rotate_x<Type> c_rtx;

public:
  /// @brief Constructor to specify 3d vector v
  /// @param[in] v : the 3d vector
  Align_y(V3d<Type> v)
      : c_rtz(Rotate_z<Type>(v.ra().rad() - M_PI / 2.0)),
        c_rtx(Rotate_x<Type>(-v.dec().rad())) {};
  /// @brief Return 3d matrix to rotate around z axis
  /// @return M3d<Type> rtz() : 3d matrix to rotate around z axis
  M3d<Type> rtz() const { return c_rtz.m3d(); };
  /// @brief Return 3d matrix to rotate around x axis
  /// @return M3d<Type> rtz() : 3d matrix to rotate around x axis
  M3d<Type> rtx() const { return c_rtx.m3d(); };
  /// @brief Return 3d matrix to align direction to v
  /// @return M3d<Type> m3d() : 3d matrix to align direction to v
  M3d<Type> m3d() const { return c_rtx.m3d() * c_rtz.m3d(); }
  /// @brief Return converted 3d vector to the new coordinats
  /// @param[in] vv : a 3d vector
  /// @return vv in the new coordinate
  V3d<Type> align(V3d<Type> vv) const { return this->m3d() * vv; };
};

/**
 * template class Align_z
 * Type : float,  double,  long double
 * @brief Class to express 3d rotation to align z-axis to the direction of the
 * vector v
 *
 */
template <class Type> class Align_z {
private:
  /// @brief Rotation around z
  Rotate_z<Type> c_rtz;
  /// @brief Rotation around y
  Rotate_y<Type> c_rty;

public:
  /// @brief Constructor to specify 3d vector v
  /// @param[in] v : the 3d vector
  Align_z(V3d<Type> v)
      : c_rtz(Rotate_z<Type>(v.ra())),
        c_rty(Rotate_y<Type>(v.dec().rad() - M_PI / 2.0)) {};
  /// @brief Return 3d matrix to rotate around z axis
  /// @return M3d<Type> rtz() : 3d matrix to rotate around z axis
  M3d<Type> rtz() const { return c_rtz.m3d(); };
  /// @brief Return 3d matrix to rotate around y axis
  /// @return M3d<Type> rtz() : 3d matrix to rotate around y axis
  M3d<Type> rty() const { return c_rty.m3d(); };
  /// @brief Return 3d matrix to align direction to v
  /// @return M3d<Type> m3d() : 3d matrix to align direction to v
  M3d<Type> m3d() const { return c_rty.m3d() * c_rtz.m3d(); };
  /// @brief Return converted 3d vector to the new coordinats
  /// @param[in] vv : a 3d vector
  /// @return vv in the new coordinate
  V3d<Type> align(V3d<Type> vv) const { return this->m3d() * vv; };
};

// template <typename Type> using EQtoEC = Rotate_x<Type>(Angle<Type>(IncEC));

/**
 * template class DirCos
 * Type : float,  double,  long double
 * @brief Direction cosine
 *
 */
template <class Type> class DirCos : public V3d<Type> {
private:
public:
  /// @brief Constructor specifying longitude and latitude
  /// @param[in] lon : longitude
  /// @param[in] lat : latitude
  DirCos(Angle<Type> lon, Angle<Type> lat)
      : V3d<Type>(std::cos(lat.rad()) * std::cos(lon.rad()),
                  std::cos(lat.rad()) * std::sin(lon.rad()),
                  std::sin(lat.rad())) {};
  /// @brief Constructor specifying a vector vv
  /// @param[in] vv : a 3d vector
  DirCos(V3d<Type> vv)
      : V3d<Type>(vv.unit().x(), vv.unit().y(), vv.unit().z()) {};
  /// @brief Return latitude
  /// @return Angle<Type> lat() : Latitude
  Angle<Type> lat() const { return Angle<Type>(std::asin(V3d<Type>::z())); };
  /// @brief Return longitude
  /// @return Angle<Type> lon() : longitude
  Angle<Type> lon() const {
    if (V3d<Type>::y() > 0.0)
      return std::acos(V3d<Type>::x() / std::cos(this->lat().rad()));
    return M_PI - std::acos(V3d<Type>::x() / std::cos(this->lat().rad()));
  };
};

/**
 * template Celestial Coordinates
 * Type : float,  double,  long double
 * @brief Celestial Coordinate
 *
 */
template <class Type> class CelCoord {
private:
  /// @brief Longitude
  Angle<Type> c_lon;
  /// @brief Latitude
  Angle<Type> c_lat;

public:
  /// @brief Blank constructor (all zero)
  CelCoord() : c_lon(Angle<Type>(0.0)), c_lat(Angle<Type>(0.0)) {};
  /// @brief Constructor specifying longitude and latitude
  /// @param[in] lon : Longitude
  /// @param[in] lat : Latitude
  CelCoord(Angle<Type> lon, Angle<Type> lat) : c_lon(lon), c_lat(lat) {};
  /// @brief Constructor specifying right ascension and declination
  /// @param[in] ra : right ascension
  /// @param[in] dc : declination
  CelCoord(RA<Type> ra, Dec<Type> dc) : c_lon(ra), c_lat(dc) {};
  /// @brief Constructor specifying direction cosine
  /// @param[in] dc : direction cosine
  CelCoord(DirCos<Type> dc) : c_lon(dc.lon()), c_lat(dc.lat()) {};
  /// @brief Return longitude
  /// @return Angle<Type> lon() : longitude
  Angle<Type> lon() const { return c_lon; };
  /// @brief Return latitude
  /// @return Angle<Type> lat() : latitude
  Angle<Type> lat() const { return c_lat; };
  /// @brief Return direction cosine
  /// @return DirCos<Type> dircos() : direction cosine
  DirCos<Type> dircos() const { return DirCos<Type>(c_lon, c_lat); };
  /// @brief Return converted Celestial coordinate from Equatorial to Ecliptic
  /// @return CelCoord<Type> EQtoEC() : converted Celestial coordinate from
  /// Equatorial to Ecliptic
  CelCoord<Type> EQtoEC() const { return CelCoord(this->dircos().EQtoEC()); };
  /// @brief Return converted Celestial coordinate from Ecliptic to Equatorial
  /// @return CelCoord<Type> ECtoEQ() : converted Celestial coordinate from
  /// Ecliptic to Equatorial
  CelCoord<Type> ECtoEQ() const { return CelCoord(this->dircos().ECtoEQ()); };
};

// template <typename Type> V3d<Type> EQtoEC(V3d<Type> eq) {
//   return Rotate_x(Angle<Type>(IncEC<Type>)) * eq;
// };
// template <typename Type> CelCoord<Type> EQtoEC(CelCoord<Type> eq) {
//   return CelCoord<Type>(EQtoEC<Type>(eq.dircos()));
// };
// template <typename Type> V3d<Type> ECtoEQ(V3d<Type> eq) {
//   return Rotate_x(Angle<Type>(-IncEC<Type>)) * eq;
// };
// template <typename Type> CelCoord<Type> ECtoEQ(CelCoord<Type> eq) {
//   return CelCoord<Type>(ECtoEQ<Type>(eq.dircos()));
// };

// template <typename Type> using ECtoEQ =
// Rotate_x<Type>(Angle<Type>(-1.0*IncEC));

//
//
//  Numeric exceptions
//
/// @brief Numerical limit
template <typename Type>
const Type Numeric_NaN = std::numeric_limits<Type>::quiet_NaN();

/// @brief Numerical limit for V2d
template <typename Type>
const V2d<Type> V2d_NaN(Numeric_NaN<Type>, Numeric_NaN<Type>);

/// @brief Numerical limit for M2d
template <typename Type>
const M2d<Type> M2d_NaN(Numeric_NaN<Type>, Numeric_NaN<Type>, Numeric_NaN<Type>,
                        Numeric_NaN<Type>);

/// @brief Numerical limit for Rectangle
template <typename Type>
const Rectangle<Type> Rectangle_NaN(V2d_NaN<Type>, V2d_NaN<Type>);

/// @brief Numerical limit for Circle
template <typename Type>
const Circle<Type> Circle_NaN(V2d_NaN<Type>, Numeric_NaN<Type>);

/// @brief Numerical limit
template <typename Type>
const Ellipse<Type> Ellipse_NaN(V2d_NaN<Type>, Numeric_NaN<Type>,
                                Numeric_NaN<Type>);

//
//  Function prototypes
//
template <typename To, typename From> V2d<To> conv(V2d<From> p);
template <typename To, typename From> M2d<To> conv(M2d<From> p);
template <typename To, typename From> Rectangle<To> conv(Rectangle<From> p);
template <typename To, typename From> Circle<To> conv(Circle<From> p);
template <typename To, typename From> Ellipse<To> conv(Ellipse<From> p);

/**
 * Square<Type> include_all(std::vector<Square<Type>> sqs)
 * @brief Return a square that include all of the squares sqs
 * @return Square<Type> include_all() : a square that include all of the squares
 * sqs
 * @param[in] sqs : a set of squares
 *
 */
template <typename Type>
Square<Type> include_all(std::vector<Square<Type>> sqs);

//
//  Suppress instantiations of templates of classes
//
//   For float
//
extern template class Angle<float>;
extern template class RA<float>;
extern template class Dec<float>;
extern template class Longitude<float>;
extern template class Latitude<float>;
extern template class V2d<float>;
extern template class Seg2d<float>;
extern template class M2d<float>;
extern template class Triangle<float>;
extern template class Rectangle<float>;
extern template class Square<float>;
extern template class Circle<float>;
extern template class Ellipse<float>;
extern template class V3d<float>;
extern template class M3d<float>;
extern template class Rotate_x<float>;
extern template class Rotate_y<float>;
extern template class Rotate_z<float>;
extern template class Align_z<float>;
extern template class DirCos<float>;

//  Define simplified class expressions
/// @brief using fAngle = Angle<float>;
using fAngle = Angle<float>;
/// @brief using fRA = RA<float>;
using fRA = RA<float>;
/// @brief using fDec = Dec<float>;
using fDec = Dec<float>;
/// @brief using fLongitude = Longitude<float>;
using fLongitude = Longitude<float>;
/// @brief using fLatitude = Latitude<float>;
using fLatitude = Latitude<float>;
/// @brief using fV2d = V2d<float>;
using fV2d = V2d<float>;
/// @brief using fSeg = Seg2d<float>;
using fSeg = Seg2d<float>;
/// @brief using fM2d = M2d<float>;
using fM2d = M2d<float>;
/// @brief using fTriangle = Triangle<float>;
using fTriangle = Triangle<float>;
/// @brief using fRectangle = Rectangle<float>;
using fRectangle = Rectangle<float>;
/// @brief using fSquare = Square<float>;
using fSquare = Square<float>;
/// @brief using fCircle = Circle<float>;
using fCircle = Circle<float>;
/// @brief using fEllipse = Ellipse<float>;
using fEllipse = Ellipse<float>;
/// @brief using fV3d = V3d<float>;
using fV3d = V3d<float>;
/// @brief using fM3d = M3d<float>;
using fM3d = M3d<float>;
/// @brief using fRotate_x = Rotate_x<float>;
using fRotate_x = Rotate_x<float>;
/// @brief using fRotate_y = Rotate_y<float>;
using fRotate_y = Rotate_y<float>;
/// @brief using fRotate_z = Rotate_z<float>;
using fRotate_z = Rotate_z<float>;
/// @brief using fAlign_x = Align_x<float>;
using fAlign_x = Align_x<float>;
/// @brief using fAlign_y = Align_y<float>;
using fAlign_y = Align_y<float>;
/// @brief using fAlign_z = Align_z<float>;
using fAlign_z = Align_z<float>;
/// @brief using fDirCos = DirCos<float>;
using fDirCos = DirCos<float>;

//
//  For double
//
extern template class Angle<double>;
extern template class RA<double>;
extern template class Longitude<double>;
extern template class Latitude<double>;
extern template class Latitude<double>;
extern template class V2d<double>;
extern template class Seg2d<double>;
extern template class M2d<double>;
extern template class Triangle<double>;
extern template class Rectangle<double>;
extern template class Square<double>;
extern template class Circle<double>;
extern template class Ellipse<double>;
extern template class V3d<double>;
extern template class M3d<double>;
extern template class Rotate_x<double>;
extern template class Rotate_y<double>;
extern template class Rotate_z<double>;
extern template class Align_x<double>;
extern template class Align_y<double>;
extern template class Align_z<double>;
extern template class DirCos<double>;

/// @brief using dAngle = Angle<double>;
using dAngle = Angle<double>;
/// @brief using dRA = RA<double>;
using dRA = RA<double>;
/// @brief using dDec = Dec<double>;
using dDec = Dec<double>;
/// @brief using dLongitude = Longitude<double>;
using dLongitude = Longitude<double>;
/// @brief using dLatitude = Latitude<double>;
using dLatitude = Latitude<double>;
/// @brief using dV2d = V2d<double>;
using dV2d = V2d<double>;
/// @brief using dSeg2d = Seg2d<double>;
using dSeg2d = Seg2d<double>;
/// @brief using dM2d = M2d<double>;
using dM2d = M2d<double>;
/// @brief using dTriangle = Triangle<double>;
using dTriangle = Triangle<double>;
/// @brief using dRectangle = Rectangle<double>;
using dRectangle = Rectangle<double>;
/// @brief using dSquare = Square<double>;
using dSquare = Square<double>;
/// @brief using dCircle = Circle<double>;
using dCircle = Circle<double>;
/// @brief using dEllipse = Ellipse<double>;
using dEllipse = Ellipse<double>;
/// @brief using dV3d = V3d<double>;
using dV3d = V3d<double>;
/// @brief using dM3d = M3d<double>;
using dM3d = M3d<double>;
/// @brief using dRotate_x = Rotate_x<double>;
using dRotate_x = Rotate_x<double>;
/// @brief using dRotate_y = Rotate_y<double>;
using dRotate_y = Rotate_y<double>;
/// @brief using dRotate_z = Rotate_z<double>;
using dRotate_z = Rotate_z<double>;
/// @brief using dAlign_x = Align_x<double>;
using dAlign_x = Align_x<double>;
/// @brief using dAlign_y = Align_y<double>;
using dAlign_y = Align_y<double>;
/// @brief using dAlign_z = Align_z<double>;
using dAlign_z = Align_z<double>;
/// @brief using dDirCos = DirCos<double>;
using dDirCos = DirCos<double>;

//
extern template class Angle<long double>;
extern template class RA<long double>;
extern template class Longitude<long double>;
extern template class Latitude<long double>;
extern template class Latitude<long double>;
extern template class V2d<long double>;
extern template class Seg2d<long double>;
extern template class M2d<long double>;
extern template class Triangle<long double>;
extern template class Rectangle<long double>;
extern template class Square<long double>;
extern template class Circle<long double>;
extern template class Ellipse<long double>;
extern template class V3d<long double>;
extern template class M3d<long double>;
extern template class Rotate_x<long double>;
extern template class Rotate_y<long double>;
extern template class Rotate_z<long double>;
extern template class Align_x<long double>;
extern template class Align_y<long double>;
extern template class Align_z<long double>;
extern template class DirCos<long double>;

/// @brief using lAngle = Angle<long double>;
using lAngle = Angle<long double>;
/// @brief using lRA = RA<long double>;
using lRA = RA<long double>;
/// @brief using lDec = Dec<long double>;
using lDec = Dec<long double>;
/// @brief using lLongitude = Longitude<long double>;
using lLongitude = Longitude<long double>;
/// @brief using lLatitude = Latitude<long double>;
using lLatitude = Latitude<long double>;
/// @brief using lV2d = V2d<long double>;
using lV2d = V2d<long double>;
/// @brief using lSeg2d = Seg2d<long double>;
using lSeg2d = Seg2d<long double>;
/// @brief using lM2d = M2d<long double>;
using lM2d = M2d<long double>;
/// @brief using lTriangle = Triangle<long double>;
using lTriangle = Triangle<long double>;
/// @brief using lRectangle = Rectangle<long double>;
using lRectangle = Rectangle<long double>;
/// @brief using lSquare = Square<long double>;
using lSquare = Square<long double>;
/// @brief using lCircle = Circle<long double>;
using lCircle = Circle<long double>;
/// @brief using lEllipse = Ellipse<long double>;
using lEllipse = Ellipse<long double>;
/// @brief using lV3d = V3d<long double>;
using lV3d = V3d<long double>;
/// @brief using lM3d = M3d<long double>;
using lM3d = M3d<long double>;
/// @brief using lRotate_x = Rotate_x<long double>;
using lRotate_x = Rotate_x<long double>;
/// @brief using lRotate_y = Rotate_y<long double>;
using lRotate_y = Rotate_y<long double>;
/// @brief using lRotate_z = Rotate_z<long double>;
using lRotate_z = Rotate_z<long double>;
/// @brief using lAlign_x = Align_x<long double>;
using lAlign_x = Align_x<long double>;
/// @brief using lAlign_y = Align_y<long double>;
using lAlign_y = Align_y<long double>;
/// @brief using lAlign_z = Align_z<long double>;
using lAlign_z = Align_z<long double>;
/// @brief using lDirCos = DirCos<long double>;
using lDirCos = DirCos<long double>;

#endif
