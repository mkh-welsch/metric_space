#ifndef _METRIC_DISTANCE_EMD_HPP
#define _METRIC_DISTANCE_EMD_HPP

namespace metric
{

namespace distance
{
/*** structural similartiy (for images) ***/
template <typename V>
struct EMD{
  using value_type = V;
  using distance_type = value_type;

  std::vector< std::vector< value_type > > C;
  value_type extra_mass_penalty = -1;
  std::vector< std::vector< value_type > > * F;


  explicit EMD(std::vector<std::vector<value_type>> && C_): C(C_) {}
  
  EMD(const std::vector<std::vector<value_type>> & C_,
      const value_type & extra_mass_penalty_ = -1, std::vector<std::vector<value_type>> *F_ = nullptr):
      C(C_), extra_mass_penalty(extra_mass_penalty_), F(F_) {}

  template<typename Container>
  distance_type  operator()(const Container &Pc, const Container &Qc) const;

};

}

}

#include "EMD.cpp"


#endif // Header Guard
