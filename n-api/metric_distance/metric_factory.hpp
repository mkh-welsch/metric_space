#ifndef METRIC_FACTORY_H_GUARD
#define METRIC_FACTORY_H_GUARD
#include <utility>
namespace metric {
namespace distance {

template<typename M>
struct Metric {
  using distance_type = typename M::distance_type;
  using record_type = typename M::record_type;
  M m;
  explicit Metric(M && m_):m(std::move(m_)) {}


  distance_type operator()(const record_type & lhs, const record_type & rhs) {
    return m(lhs,rhs);
  }
};

template<typename M, typename... Args>
M make_metric(Args&&... args) {
  return M(std::forward<Args>(args)...);
}

} // namespace distance
} // namespace metric
#endif
