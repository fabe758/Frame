
/**
 * @file Lcurve.cpp
 +
 * @brief Class library to handle light curve
 * @author Fumio Abe
 * @date 1 February,  2024
 *
 * @details  Theoretical light curve generator
 * @note
 *
 */

/**
 * include headers  for this project
 * all c++ headers are loaded in Geom.h
 *
 */
#include "Lcurve.h" /// for Lcurve

template <typename Type> Source<Type> Lcurve<Type>::src(int it) const {
  return c_srcm[0].src(c_t[it]);
}

template <typename Type> Source<Type> Lcurve<Type>::src(Time<Type> t) const {
  return c_srcm[0].src(t);
}

template <typename Type> Mlens<Type> Lcurve<Type>::ml(int it) const {
  return c_mlm.ml(c_t[it]);
}

template <typename Type> Mlens<Type> Lcurve<Type>::ml(Time<Type> t) const {
  return c_mlm.ml(t);
}

template <typename Type> FrTri<Type> Lcurve<Type>::frtri(int it) const {
  return FrTri<Type>(this->ml(it), this->src(it), c_param);
}

template <typename Type> FrTri<Type> Lcurve<Type>::frtri(Time<Type> t) const {
  return FrTri<Type>(this->ml(t), this->src(t), c_param);
}

template <class Type> void Lcurve<Type>::process(int round) {
  if (c_lcp.bench)
    c_lcp.start = std::clock();
  c_mag = std::vector<std::vector<Type>>(c_srcm.size(),
                                         std::vector<Type>(c_t.size()));
  c_as = std::vector<V2d<Type>>(c_t.size());
  std::vector<Type> br_tot(c_t.size());
  for (size_t is = 0; is < c_srcm.size(); is++) {
    if (c_lcp.record) {
      c_lcp.in_size.push_back(std::vector<int>(c_t.size()));
      c_lcp.unc_size.push_back(std::vector<int>(c_t.size()));
      c_lcp.out_size.push_back(std::vector<int>(c_t.size()));
      c_lcp.in_area.push_back(std::vector<Type>(c_t.size()));
      c_lcp.unc_area.push_back(std::vector<Type>(c_t.size()));
      c_lcp.area_sn.push_back(std::vector<Type>(c_t.size()));
      c_lcp.round.push_back(std::vector<int>(c_t.size()));
    }
    for (size_t it = 0; it < c_t.size(); it++) {
      FrTri<Type> frtri(c_mlm.ml(c_t[it]), c_srcm[is].src(c_t[it]), c_param);
      frtri.process(round);
      c_mag[is][it] = frtri.area_sn() / frtri.src().area();
      br_tot[it] += frtri.area_sn();
      c_as[it] += frtri.sump_area();
      if (c_lcp.record) {
        c_lcp.in_size[is][it] = frtri.in().size();
        c_lcp.unc_size[is][it] = frtri.uncertain().size();
        c_lcp.out_size[is][it] = frtri.out().size();
        c_lcp.in_area[is][it] = frtri.in().area();
        c_lcp.unc_area[is][it] = frtri.uncertain().area();
        c_lcp.area_sn[is][it] = frtri.area_sn();
        c_lcp.round[is][it] = frtri.round();
      }
      //   for (size_t is = 0; is < c_srcm.size(); is++) {
      //     if (is >= c_frtri.size())
      //       c_frtri.push_back({});
      //     for (size_t it = 0; it < c_t.size(); it++) {
      //       c_frtri[is].push_back(
      //           FrTri<Type>(c_mlm.ml(c_t[it]), c_srcm[is].src(c_t[it]),
      //           c_param));
      //     }
      //   }
      //   for (size_t is = 0; is < c_srcm.size(); is++)
      //     for (size_t it = 0; it < c_t.size(); it++)
      //       c_frtri[is][it].process(round);
      for (size_t it = 0; it < c_t.size(); it++)
        c_as[it] = c_as[it] / br_tot[it];
      if (c_lcp.bench)
        c_lcp.stop = std::clock();
    }
  }
}

template <typename Type> void Lcurve<Type>::process_while(int in_size) {
  if (c_lcp.bench)
    c_lcp.start = std::clock();
  c_mag = std::vector<std::vector<Type>>(c_srcm.size(),
                                         std::vector<Type>(c_t.size()));
  c_as = std::vector<V2d<Type>>(c_t.size());
  std::vector<Type> br_tot(c_t.size());
  for (size_t is = 0; is < c_srcm.size(); is++) {
    if (c_lcp.record) {
      c_lcp.in_size.push_back(std::vector<int>(c_t.size()));
      c_lcp.unc_size.push_back(std::vector<int>(c_t.size()));
      c_lcp.out_size.push_back(std::vector<int>(c_t.size()));
      c_lcp.in_area.push_back(std::vector<Type>(c_t.size()));
      c_lcp.unc_area.push_back(std::vector<Type>(c_t.size()));
      c_lcp.area_sn.push_back(std::vector<Type>(c_t.size()));
      c_lcp.round.push_back(std::vector<int>(c_t.size()));
    }
    for (size_t it = 0; it < c_t.size(); it++) {
      FrTri<Type> frtri(c_mlm.ml(c_t[it]), c_srcm[is].src(c_t[it]), c_param);
      frtri.process_while(in_size);
      c_mag[is][it] = frtri.area_sn() / frtri.src().area();
      br_tot[it] += frtri.area_sn();
      c_as[it] += frtri.sump_area();
      if (c_lcp.record) {
        c_lcp.in_size[is][it] = frtri.in().size();
        c_lcp.unc_size[is][it] = frtri.uncertain().size();
        c_lcp.out_size[is][it] = frtri.out().size();
        c_lcp.in_area[is][it] = frtri.in().area();
        c_lcp.unc_area[is][it] = frtri.uncertain().area();
        c_lcp.area_sn[is][it] = frtri.area_sn();
        c_lcp.round[is][it] = frtri.round();
      }
    }
  }
  for (size_t it = 0; it < c_t.size(); it++)
    c_as[it] = c_as[it] / br_tot[it];
  if (c_lcp.bench)
    c_lcp.stop = std::clock();
}

template <typename Type> void Lcurve<Type>::process_while(Type ratio_unc_in) {
  if (c_lcp.bench)
    c_lcp.start = std::clock();
  c_mag = std::vector<std::vector<Type>>(c_srcm.size(),
                                         std::vector<Type>(c_t.size()));
  c_as = std::vector<V2d<Type>>(c_t.size());
  std::vector<Type> br_tot(c_t.size());
  for (size_t is = 0; is < c_srcm.size(); is++) {
    for (size_t it = 0; it < c_t.size(); it++) {
      FrTri<Type> frtri(c_mlm.ml(c_t[it]), c_srcm[is].src(c_t[it]), c_param);
      frtri.process_while(ratio_unc_in);
      c_mag[is][it] = frtri.area_sn() / frtri.src().area();
      br_tot[it] += frtri.area_sn();
      c_as[it] += frtri.sump_area();
    }
  }
  for (size_t it = 0; it < c_t.size(); it++)
    c_as[it] = c_as[it] / br_tot[it];
  if (c_lcp.bench)
    c_lcp.stop = std::clock();
}

template <typename Type>
std::vector<std::vector<Type>> Lcurve<Type>::mags() const {
  return c_mag;
}

template <typename Type> std::vector<Type> Lcurve<Type>::mag() const {
  return c_mag[0];
}

template <typename Type> std::vector<V2d<Type>> Lcurve<Type>::as() const {
  return c_as;
}

// template <typename Type>
// std::vector<std::vector<Type>> Lcurve<Type>::mags() const {
//   std::vector<std::vector<Type>> mag0;
//   for (size_t is = 0; is < c_srcm.size(); is++) {
//     if (is >= mag0.size())
//       mag0.push_back({});
//     for (size_t it = 0; it < c_t.size(); it++) {
//       mag0[is].push_back(c_frtri[is][it].area_sn() /
//                          c_frtri[is][it].src().area());
//     }
//   }
//   return mag0;
// };

// template <typename Type> std::vector<Type> Lcurve<Type>::mag() const {
//   return this->mags()[0];
// };

// template <typename Type> std::vector<V2d<Type>> Lcurve<Type>::as() const
// {
//   std::vector<V2d<Type>> as;
//   std::vector<Type> br_tot;
//   for (size_t it = 0; it < c_t.size(); it++) {
//     as.push_back(V2d_0<Type>);
//     br_tot.push_back(0.0);
//   }
//   for (size_t is = 0; is < c_srcm.size(); is++) {
//     for (size_t it = 0; it < c_t.size(); it++) {
//       br_tot[it] += c_frtri[is][it].area_sn();
//       as[it] += c_frtri[is][it].sump_area();
//     }
//   }
//   for (size_t it = 0; it < c_t.size(); it++)
//     as[it] = as[it] / br_tot[it];
//   return as;
// };

template <typename Type>
std::vector<std::vector<Snap<Type>>> Lcurve<Type>::snap(size_t min,
                                                        size_t max) const {
  std::vector<std::vector<Snap<Type>>> snap;
  for (size_t it = min; it < max; it++)
    for (size_t is = 0; is < c_srcm.size(); is++) {
      // if (is >= snap.size())
      //   snap.push_back({});
      // snap[is].push_back(c_frtri[is][it].snap());
    }
  return snap;
};

template <typename Type>
std::vector<std::vector<Snap<Type>>>
Lcurve<Type>::snap(std::function<bool(Time<Type>)> select) const {
  std::vector<std::vector<Snap<Type>>> snap;
  for (size_t is = 0; is < c_srcm.size(); is++) {
    if (is >= snap.size())
      snap.push_back({});
    // for (size_t it = 0; it < c_t.size(); it++)
    //   if (select(c_t[it]))
    //     snap[is].push_back(c_frtri[is][it].snap());
  }
  return snap;
};

template <typename Type>
std::vector<std::vector<Snap<Type>>>
Lcurve<Type>::snap(std::function<bool(FrTri<Type>)> select) const {
  std::vector<std::vector<Snap<Type>>> snap;
  for (size_t it = 0; it < c_t.size(); it++) {
    bool take = false;
    // for (size_t is = 0; is < c_srcm.size(); is++)
    //   if (select(c_frtri[is][it]))
    //     take = true;
    if (take)
      for (size_t is = 0; is < c_srcm.size(); is++) {
        // if (is >= snap.size())
        //   snap.push_back({});
        // snap[is].push_back(c_frtri[is][it].snap());
      };
  }
  return snap;
};

template <typename Type> Type Lcurve<Type>::ctime() const {
  return static_cast<Type>(c_lcp.stop - c_lcp.start) / CLOCKS_PER_SEC;
}

//
//  Explicit instantiations of templates of classes
//
//   For float
//
template struct LcParam<float>;
template class Lcurve<float>;

//
//  For double
//
template struct LcParam<double>;
template class Lcurve<double>;

//
//  For long double
//
template struct LcParam<long double>;
template class Lcurve<long double>;
