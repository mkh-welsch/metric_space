#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_api
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "metric_search.hpp"
template<typename T>
struct distance {
    int operator()( const T &lhs,  const T &rhs) const {
        return std::abs(lhs - rhs);
    }
};

BOOST_AUTO_TEST_CASE(test_insert) {
    std::vector<int> data = {3,5,-10,50,1,-200,200};
    metric_search::Tree<int,distance<int>> tree;
    for(auto t : data) {
        tree.insert(t);
        BOOST_TEST(tree.check_covering());
    }
}

BOOST_AUTO_TEST_CASE(test_insert_batch) {
    std::vector<int> data = {3,5,-10,50,1,-200,200};
    metric_search::Tree<int,distance<int>> tree;
    tree.insert(data);
    BOOST_TEST(tree.check_covering());
}
BOOST_AUTO_TEST_CASE(test_nn) {
    std::vector<int> data = {3,5,-10,50,1,-200,200};
    metric_search::Tree<int,distance<int>> tree;
    tree.insert(data);
    tree.print();
    BOOST_TEST(tree.nn(200)->data == 200);
    // for(auto d : data) {
    //     auto p = tree.nn(d);
    //     BOOST_TEST(p->data == d);
    // }
}

BOOST_AUTO_TEST_CASE(test_knn) {
    std::vector<int> data = {3,5,-10,50,1,-200,200};
    metric_search::Tree<int,distance<int>> tree;
    tree.insert(data);
    auto k1 = tree.knn(3,15);
    BOOST_TEST(k1.size() == 7);
    BOOST_TEST(k1[0].first->data == 3);
    BOOST_TEST(k1[1].first->data == 1);
    BOOST_TEST(k1[2].first->data == 5);
    BOOST_TEST(k1[3].first->data == -10);
    BOOST_TEST(k1[4].first->data == 50);
    BOOST_TEST(k1[5].first->data == 200);
    BOOST_TEST(k1[6].first->data == -200);
}

BOOST_AUTO_TEST_CASE(test_erase) {
    std::vector<int> data = {3,5,-10,50,1,-200,200};
    metric_search::Tree<int,distance<int>> tree;
    tree.insert(data);
    //    tree.print();
    for(auto d : data) {
        tree.erase(d);
        BOOST_TEST(tree.check_covering());
        //        tree.print();
    }
}

BOOST_AUTO_TEST_CASE(test_erase_root) {
    std::vector<int> data = {3,5,-10,50,1,-200,200};
    metric_search::Tree<int,distance<int>> tree;
    tree.insert(data);
    //    tree.print();
    for(int i = 0; i < 7; i++) {
        auto root = tree.get_root();
        tree.erase(root->data);
        BOOST_TEST(tree.check_covering());
        //        tree.print();
    }
}

BOOST_AUTO_TEST_CASE(test_insert_if) {
    metric_search::Tree<int,distance<int>> tree;
    tree.insert(1);
    BOOST_TEST(!tree.insert_if(2,10));
    BOOST_TEST(tree.insert_if(15,10));
    BOOST_TEST(!tree.insert_if(14,10));
    BOOST_TEST(tree.insert_if(26,10));
}