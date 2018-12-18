#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_cluster
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <vector>
#include "metric_search.hpp"
template<typename T>
struct distance {
  T operator()( const T &lhs,  const T &rhs) const {
    return std::abs(lhs - rhs);
  }
};

BOOST_AUTO_TEST_CASE(ccc) {
  std::vector<int> data = {7,8,9,10,11,12,13};
  metric_search::Tree<int,distance<int>> tree;
  tree.insert(data);
  tree.print();
  std::vector<double> distribution = {0.1, 0.2, 0.3, 0.5};
    std::vector<std::size_t>  IDS = {1,2,3};
  auto result = tree.clastering(distribution, IDS, data);
  std::size_t i = 0;
  for(auto & v : result) {
    std::cout << distribution[i] << " = {";
    for(auto p : v) {
      std::cout << "[" << data[p] << ":" << p << "]"  << ", ";
    }
    std::cout << "}" << std::endl;
    i++;
  }
}
