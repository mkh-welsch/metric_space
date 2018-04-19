/* Signal Epowering Technology   
                       
presents

███╗   ███╗███████╗████████╗██████╗ ██╗ ██████╗    ███████╗███████╗ █████╗ ██████╗  ██████╗██╗  ██╗
████╗ ████║██╔════╝╚══██╔══╝██╔══██╗██║██╔════╝    ██╔════╝██╔════╝██╔══██╗██╔══██╗██╔════╝██║  ██║
██╔████╔██║█████╗     ██║   ██████╔╝██║██║         ███████╗█████╗  ███████║██████╔╝██║     ███████║
██║╚██╔╝██║██╔══╝     ██║   ██╔══██╗██║██║         ╚════██║██╔══╝  ██╔══██║██╔══██╗██║     ██╔══██║
██║ ╚═╝ ██║███████╗   ██║   ██║  ██║██║╚██████╗    ███████║███████╗██║  ██║██║  ██║╚██████╗██║  ██║
╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝ ╚═════╝    ╚══════╝╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝
                                                                            Licensed under MPL 2.0. 
                                                                            Michael Welsch (c) 2018.
                                                                                                   
a library for metric search algorithmen and data containers

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Short summary of the MPL 2.0 intention:
Use it as you like, but if you make changes, you must make the source 
code available under MPL again. So, keep the MPL Software in separated 
files, to get not in conflict with other licences. */

/*** Usage with default L2 metric for abitrary stl containers. ***/

/*
#include "metric_search.hpp"

int main()
{
// here are some data records
std::vector<Distance> v0{0, 1, 1, 1, 1, 1, 2, 3};
std::vector<Distance> v1{1, 1, 1, 1, 1, 2, 3, 4};
std::vector<Distance> v2{2, 2, 2, 1, 1, 2, 0, 0};
std::vector<Distance> v3{3, 3, 2, 2, 1, 1, 0,0};
std::vector<Distance> v4{4, 3, 2, 1, 0, 0, 0, 0};
std::vector<Distance> v5{5, 3, 2, 1, 0, 0, 0, 0};
std::vector<std::vector<Distance>> datafield;
datafield.push_back(v0); datafield.push_back(v1); datafield.push_back(v2); datafield.push_back(v3); datafield.push_back(v4); datafield.push_back(v5);

//initialize the tree
metric_search::Tree<vector<Distance>> cTree;
// metric_search::Tree<vector<Distance>> cTree(v0); // alternative
// metric_search::Tree<vector<Distance>> cTree(datafield); // alternative

// add data records
cTree.push(v1); // pushes the data records into the tree

// add data record under user definied conditions
auto compare_fun1 = [](distance){return (dist(_p) =! distance);}; // record allready exists (default behaviour)
auto compare_fun2 = [](distance){return (dist(_p) < distance);}; // insert only by distance under a threshold.
bool idx cTree.insert(v2, compare_fun); // inserts the data record with lambda function check, gives back the new index if iserted or the colling index of the check

// access data records
auto data_record = cTree[1] // gives back v1
cTree[1] = v3; // change record (avoid this expansiv opteration, because in the tree it can't be overriden und must be earesed and pushed again.)

//erase data record
auto data_erased = cTree.erase(5); // erases the 5. data record and give it back

// find data record
auto data_nn = cTree.nn(v_test); // nearest neigbour
auto vector_of_data_knn = cTree.knn(v_test,3); // k nearest neigbours
auto vector_of_data_ranged = cTree.range(v_test,v1); //gives back the data record in an std::vector.

// tree properties
auto tree_size = cTree.size();

// import data
metric_search::Tree<vector<Distance>> cTree1(container_of_data_records);
metric_search::Tree<vector<Distance>> cTree2(container_of_data_records.begin(),container_of_data_records.end());
auto vector_of_select_data = cTree.toVector(cTree.begin()+1,cTree.end()-2);

// export data
auto vector_of_complete_data = cTree.toVector();
auto vector_of_select_data = cTree.toVector(cTree.begin()+1,cTree.end()-2);

return 0;
} 
*/

/*** Use with custom Metric ***/

/*
#include <eigen3/Eigen/Core> // path to Eigen lib
#include <blaze/Math.h> // path to Blaze Lib
#include "metric_search.hpp"

// custom Metric for STL Containers.
template <typename Container>
struct customMetric
{
    typedef typename Container::value_type T;
    static_assert(
        std::is_floating_point<T>::value, "T must be a float type");

    T operator()(const Container &a, const Container &b) const
    {
        T sum = 0;
        for (auto it1 = a.begin(), it2 = b.begin(); it1 != a.end() || it2 != b.end(); ++it1, ++it2)
        {
            sum += (*it1 - *it2) * (*it1 - *it2);
        }

        return std::sqrt(sum);
    }
};

// L2 Metric for Eigen Vector
template <typename T>
struct recMetric_Eigen
{
    T operator()(const Eigen::VectorXd &p, const Eigen::VectorXd &q) const
    {
        return (p - q).norm();
    }
};
typedef Eigen::VectorXd recType_Eigen;

// L2 Metric for Blaze Vector
template <typename T>
struct recMetric_Blaze
{
    T operator()(const blaze::DynamicVector<T> &p, const blaze::DynamicVector<T> &q) const
    {
        return blaze::l2Norm(p - q);
    }
};




int main()
{
metric_search::Tree<deque<Distance>,customMetric<Distance>> cTree1;

metric_search::Tree<recType_Eigen,recMetric_Eigen> cTree2;

metric_search::Tree<blaze::DynamicVector<Distance>,recMetric_Blaze<Distance>> cTree3;

return 0;
}


*/
#ifndef _METRIC_SEARCH_HPP
#define _METRIC_SEARCH_HPP

#include <atomic>
#include <fstream>
#include <iostream>
#include <stack>
#include <map>
#include <vector>
#include <shared_mutex>
#include <numeric>
#include <cmath>
#include <string>

namespace metric_search
{

/*
  |   |                   |                     
  |   |   _ \   _` |   _` |   _ \   __|         
  ___ |   __/  (   |  (   |   __/  |            
 _|  _| \___| \__,_| \__,_| \___| _|            
                                                                                                                
*/
/*** class for each node of the tree ***/

/*
  _ \         _|             |  |       \  |        |       _)      
  |  |  -_)   _| _` |  |  |  |   _|    |\/ |   -_)   _|   _| |   _| 
 ___/ \___| _| \__,_| \_,_| _| \__|   _|  _| \___| \__| _|  _| \__| 
                                                                    
*/
/*** standard euclidian (L2) Metric ***/
template <typename Container>
struct L2_Metric_STL;

/*
 __ __|              
    |   _| -_)   -_) 
   _| _| \___| \___| 
                     
*/
/*** Cover Tree Implementation ***/
template <class recType, class Metric = L2_Metric_STL<recType>>
class Tree
{
private:
  class Node; // Node Class (see implementation for details)

  /*** Types ***/
  Metric metric_;
  typedef Tree<recType, Metric>::Node NodeType;
  typedef std::shared_ptr<Tree<recType, Metric>::Node> Node_ptr;
  typedef Tree<recType, Metric> TreeType;
  typedef typename std::result_of<Metric(recType, recType)>::type Distance;

  /*** Properties ***/
  Distance base = 2;                  // Base for estemating the covering of the tree
  Node_ptr root;                      // Root of the tree
  std::atomic<int> min_scale;         // Minimum scale
  std::atomic<int> max_scale;         // Minimum scale
  int truncate_level;                 // Relative level below which the tree is truncated
  std::atomic<unsigned> N;            // Number of points in the cover tree
  std::shared_timed_mutex global_mut; // lock for changing the root

  /*** Imlementation Methodes ***/
  template <typename pointOrNodeType>
  std::tuple<std::vector<int>, std::vector<Distance>>
  sortChildrenByDistance(Node_ptr p, pointOrNodeType x) const;

  template <typename pointOrNodeType>
  bool insert_(Node_ptr p, pointOrNodeType x, size_t new_id = -2);

  void nn_(Node_ptr current, Distance dist_current, const recType &p, std::pair<Node_ptr, Distance> &nn) const;
  void knn_(Node_ptr current, Distance dist_current, const recType &p, std::vector<std::pair<Node_ptr, Distance>> &nnList) const;
  void range_(Node_ptr current, Distance dist_current, const recType &p, Distance distance, std::vector<std::pair<Node_ptr, Distance>> &nnList) const;

  void printTree_(NodeType *node_p);

public:
  /*** Coonstructors ***/
  Tree(int truncate = -1, Metric d = Metric());                   // empty tree
  Tree(const recType &p, int truncate = -1, Metric d = Metric()); // cover tree with one data record as root
  ~Tree();                                                        // Destuctor

  /*** Access Operations ***/
  bool insert(const recType &p); // insert data record into the cover tree
  bool erase(const recType &p);  // erase data record into the cover tree
  recType operator[](size_t id); // access a data record by ID

  /*** Nearest Neighbour search ***/
Node_ptr nn(const recType &p) const;                                                // nearest Neighbour
  std::vector<std::pair<Node_ptr, Distance>> knn(const recType &p, unsigned k = 10) const;                 // k-Nearest Neighbours
  std::vector<std::pair<Node_ptr, Distance>> range(const recType &queryPt, Distance distance = 1.0) const; // Range Search

  /** Dev Tools **/
  int levelSize();                        // return the max_level of the tree (= root level)
  std::map<int, unsigned> print_levels(); // print and return level informations
  size_t size();                          // return node size.
  std::vector<recType> toVector();        // return all records in the tree in a std::vector

  bool check_covering() const;
  bool check_covering2() const;
  bool repair_covering();
  void printTree();

  void shift_level(int level);
  void traverse(const std::function<void(Node_ptr)> &f);
  void traverse_child(const std::function<void(Node_ptr)> &f);
};

} // end namespace

#include "metric_search.cpp" // include the implementation

#endif //_METRIC_SEARCH_HPP
