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
namespace std {
    std::ostream & operator << (std::ostream & ostr, const std::vector<std::size_t> & v) {
        ostr << "[ ";
        for(auto i : v) {
            ostr << i << ", ";
        }
        ostr <<"]";
        return ostr;
    }
}
BOOST_AUTO_TEST_CASE(cluster1) {
    std::vector<int> data = {7,8,9,10,11,12,13};
    metric_search::Tree<int,distance<int>> tree;
    tree.insert(data);
    tree.print();
    std::vector<double> distribution = {0.1, 0.2, 0.3, 0.5};
    std::vector<std::size_t>  IDS = {1,2,3};
    std::vector<int> points = {8,9,10};

    auto result = tree.clustering(distribution, IDS, data);
    auto result2 = tree.clustering(distribution, points);
    std::vector<std::vector<std::size_t>> test_result = {{}, {1}, {0}, {2}};
    BOOST_TEST(result == test_result, boost::test_tools::per_element());
    BOOST_TEST(result2 == test_result, boost::test_tools::per_element());
    // std::size_t i = 0;
    // for(auto & v : result) {
    //   std::cout << distribution[i] << " = {";
    //   for(auto p : v) {
    //     std::cout << "[" << data[p] << ":" << p << "]"  << ", ";
    //   }
    //   std::cout << "}" << std::endl;
    //   i++;
    // }
}

BOOST_AUTO_TEST_CASE(cluster2) {
    std::vector<int> data = {7,8,9,10,11,12,13};
    metric_search::Tree<int,distance<int>> tree;
    tree.insert(data);
    tree.print();
    std::vector<double> distribution = {0.1, 0.2, 0.3, 0.5};
    std::vector<std::size_t>  IDS = {3};
    std::vector<int>  points = {10};
    auto result = tree.clustering(distribution, IDS, data);
    auto result2 = tree.clustering(distribution, points);
    std::vector<std::vector<std::size_t>> test_result = {{}, {3}, {4}, {2}};
    BOOST_TEST(result == test_result, boost::test_tools::per_element());
    BOOST_TEST(result2 == test_result, boost::test_tools::per_element());
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

BOOST_AUTO_TEST_CASE(cluster3) {
    std::vector<int> data = {7,8,9,10,11,12,13};
    metric_search::Tree<int,distance<int>> tree;
    tree.insert(data);
    tree.print();
    std::vector<double> distribution = {0.1, 0.2, 0.5, 0.9};
    std::vector<std::size_t>  IDS = {3};
    std::vector<int>  points = {10};
    auto result = tree.clustering(distribution, IDS, data);
    auto result2 = tree.clustering(distribution, points);
    std::vector<std::vector<std::size_t>> test_result = {{}, {3}, {4,2}, {1,0,5}};
    BOOST_TEST(result == test_result, boost::test_tools::per_element());
    BOOST_TEST(result2 == test_result, boost::test_tools::per_element());
    // std::size_t i = 0;
    // for(auto & v : result) {
    //     std::cout << distribution[i] << " = {";
    //     for(auto p : v) {
    //         std::cout << "[" << data[p] << ":" << p << "]"  << ", ";
    //     }
    //     std::cout << "}" << std::endl;
    //     i++;
    // }
}
