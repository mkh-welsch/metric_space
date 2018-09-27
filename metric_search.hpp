/* Signal Empowering Technology   
                       
presents

███╗   ███╗███████╗████████╗██████╗ ██╗ ██████╗    ███████╗███████╗ █████╗ ██████╗  ██████╗██╗  ██╗
████╗ ████║██╔════╝╚══██╔══╝██╔══██╗██║██╔════╝    ██╔════╝██╔════╝██╔══██╗██╔══██╗██╔════╝██║  ██║
██╔████╔██║█████╗     ██║   ██████╔╝██║██║         ███████╗█████╗  ███████║██████╔╝██║     ███████║
██║╚██╔╝██║██╔══╝     ██║   ██╔══██╗██║██║         ╚════██║██╔══╝  ██╔══██║██╔══██╗██║     ██╔══██║
██║ ╚═╝ ██║███████╗   ██║   ██║  ██║██║╚██████╗    ███████║███████╗██║  ██║██║  ██║╚██████╗██║  ██║
╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝ ╚═════╝    ╚══════╝╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝
                                                                            Licensed under MPL 2.0. 
                                                                            Michael Welsch (c) 2018.
                                                                                                   
a library for metric search algorithms and data containers

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
#include <functional>
namespace metric_search
{
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
    //typedef std::shared_ptr<Tree<recType, Metric>::Node> Node_ptr;
    typedef Tree<recType,Metric>::Node* Node_ptr;
  typedef Tree<recType, Metric> TreeType;
    //  typedef typename std::result_of<Metric(recType, recType)>::type Distance;
    using Distance = typename std::result_of<Metric(recType,recType)>::type;

  /*** Properties ***/
  Distance base = 2;                  // Base for estemating the covering of the tree
  Node_ptr root;                      // Root of the tree
  std::atomic<int> min_scale;         // Minimum scale
  std::atomic<int> max_scale;         // Minimum scale
  int truncate_level;                 // Relative level below which the tree is truncated
  std::atomic<unsigned> N;            // Number of points in the cover tree
  // std::shared_timed_mutex global_mut; // lock for changing the root

  /*** Imlementation Methodes ***/
  template <typename pointOrNodeType>
  std::tuple<std::vector<int>, std::vector<Distance>>
  sortChildrenByDistance(Node_ptr p, pointOrNodeType x) const;

  template <typename pointOrNodeType>
  bool insert_(Node_ptr p, pointOrNodeType x, int new_id = -2);

    void nn_(Node_ptr current, Distance dist_current, const recType &p, std::pair<Node_ptr, Distance> &nn) const;
    std::size_t knn_(Node_ptr current, Distance dist_current, const recType &p, std::vector<std::pair<Node_ptr, Distance>> &nnList, std::size_t nnSize) const;
  void rnn_(Node_ptr current, Distance dist_current, const recType &p, Distance distance, std::vector<std::pair<Node_ptr, Distance>> &nnList) const;

  void print_(NodeType *node_p);

    Node_ptr merge(Node_ptr p, Node_ptr q);
    std::pair<Node_ptr,std::vector<Node_ptr>> mergeHelper(Node_ptr p, Node_ptr q);
    auto findAnyLeaf() -> Node_ptr;
    void extractNode(Node_ptr node);
public:
  /*** Constructors ***/
  Tree(int truncate = -1, Metric d = Metric());                                // empty tree
  Tree(const recType &p, int truncate = -1, Metric d = Metric());              // cover tree with one data record as root
  Tree(const std::vector<recType> &p, int truncate = -1, Metric d = Metric()); // with a vector of data records
  ~Tree();                                                                     // Destuctor

  /*** Access Operations ***/
    bool insert(const recType &p);              // insert data record into the cover tree
    bool insert_if(const recType &p, Distance treshold);              // insert data record into the cover tree only if distance bigger than a treshold
    std::size_t insert_if(const std::vector<recType> &p, Distance treshold); // insert data record into the cover tree
    bool insert(const std::vector<recType> &p); // insert data record into the cover tree
    bool erase(const recType &p);               // erase data record into the cover tree
    recType operator[](size_t id);              // access a data record by ID

  /*** Nearest Neighbour search ***/
  Node_ptr nn(const recType &p) const;                                                                   // nearest Neighbour
  std::vector<std::pair<Node_ptr, Distance>> knn(const recType &p, unsigned k = 10) const;               // k-Nearest Neighbours
  std::vector<std::pair<Node_ptr, Distance>> rnn(const recType &queryPt, Distance distance = 1.0) const; // Range Search

  /*** utilitys ***/
  size_t size(); // return node size.
  void traverse(const std::function<void(Node_ptr)> &f);

  /** Dev Tools **/
  int levelSize();                        // return the max_level of the tree (= root level)
  void print();
  std::map<int, unsigned> print_levels(); // print and return level informations

  std::vector<recType> toVector(); // return all records in the tree in a std::vector

  bool check_covering() const;
  void traverse_child(const std::function<void(Node_ptr)> &f);
    Node_ptr get_root() {
        return root;
    }
    bool empty() {
        return root == nullptr;
    }
    int get_root_level() {
        return root->level;
    }
};

} // end namespace

#include "metric_search.cpp" // include the implementation

#endif //_METRIC_SEARCH_HPP
