/* 
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Signal Empowering Technology ®Michael Welsch
*/

#ifndef _METRIC_GRAPH_HPP
#define _METRIC_GRAPH_HPP


#include "Blaze.h"
#include <stack>
#include <type_traits>


namespace metric {
namespace graph {


/*
// -------
// old Graph code, left for comparative testing only, TODO remove


class Graph_old {
        public:
                explicit Graph_old(size_t nodesNumber);
                ~Graph_old();

                size_t getNodesNumber();
                bool isValid();

        std::vector<std::vector<size_t>> getNeighbours(const size_t nodeIndex, const size_t maxDeep);

        protected:
                size_t nodesNumber;
                bool valid;
                std::vector<std::vector<size_t>> edges;

                size_t modularPow(const size_t base, const size_t exponent, const size_t modulus);

                void buildEdges(const std::vector<std::pair<size_t, size_t>> &edgesPairs);

                void neighboursWalk(const size_t index, std::unordered_map<size_t, size_t> &indexes, size_t deep,
                                                        const size_t maxDeep);
};

class Grid4_old : public Graph_old {
        public:
                explicit Grid4_old(size_t nodesNumber);
                Grid4_old(size_t width, size_t height);
        private:
                void construct(size_t width, size_t height);
};

class Grid6_old : public Graph_old {
public:
        explicit Grid6_old(size_t nodesNumber);
        Grid6_old(size_t width, size_t height);
private:
        void construct(size_t width, size_t height);
};

class Grid8_old : public Graph_old {
        public:
                explicit Grid8_old(size_t nodesNumber);
                Grid8_old(size_t width, size_t height);
        private:
                void construct(size_t width, size_t height);
};

class LPS_old : public Graph_old {
        public:
                explicit LPS_old(size_t nodesNumber);
        private:
                bool miller_rabin_pass(const size_t a, const size_t s,
                                                        const size_t d, const size_t nodesNumber);

                bool millerRabin(const size_t nodesNumber);
};

class Paley_old : public Graph_old {
        public:
                explicit Paley_old(size_t nodesNumber);
};

class Margulis_old : public Graph_old {
        public:
                explicit Margulis_old(size_t nodesNumber);
};



// end of old Graph code
// -------
*/




// Graph based on blaze-lib

//template <typename WeightType = bool, bool isDense = true, bool isSymmetric = true> // incorrect defaults, needed only to test dense
template <typename WeightType = bool, bool isDense = false, bool isSymmetric = true> // correct defaults
class Graph {

private:
    static constexpr bool isWeighted = !std::is_same<WeightType, bool>::value; // used only in getNeighboursOld, TODO remove or update so as in new method

    typedef typename std::conditional<
        isDense,
        blaze::DynamicMatrix<WeightType>,
        blaze::CompressedMatrix<WeightType>
        >::type InnerMatrixType;

    typedef typename std::conditional<
        isSymmetric,
        blaze::SymmetricMatrix<InnerMatrixType>,
        InnerMatrixType
        >::type MatrixType;

public:
    explicit Graph(size_t nodesNumber);
    Graph();
    Graph(const std::vector<std::pair<size_t, size_t>> & edgesPairs);
    ~Graph();

    size_t getNodesNumber();
    bool isValid();

    std::vector<std::vector<size_t>> getNeighboursOld(const size_t nodeIndex, const size_t maxDeep);

    template <typename T = WeightType, bool denseFlag = isDense>
    typename std::enable_if_t<!std::is_same<T, bool>::value /*&& !denseFlag*/, std::vector<std::vector<size_t>>>
    getNeighbours(const size_t nodeIndex, const size_t maxDeep); // not bool

    template <typename T = WeightType, bool denseFlag = isDense>
    typename std::enable_if_t<std::is_same<T, bool>::value && !denseFlag, std::vector<std::vector<size_t>>>
    getNeighbours(const size_t nodeIndex, const size_t maxDeep); // bool, not dense

    template <typename T = WeightType, bool denseFlag = isDense>
    typename std::enable_if_t<std::is_same<T, bool>::value && denseFlag, std::vector<std::vector<size_t>>>
    getNeighbours(const size_t nodeIndex, const size_t maxDeep); // bool, dense

    MatrixType get_matrix();

    void buildEdges(const std::vector<std::pair<size_t, size_t>> &edgesPairs);

protected:
    size_t nodesNumber;
    bool valid;

    MatrixType m;

    size_t modularPow(const size_t base, const size_t exponent, const size_t modulus);

//    void buildEdges(const std::vector<std::pair<size_t, size_t>> &edgesPairs); // moved to public section

};





class Grid4 : public Graph<> {
public:
    explicit Grid4(size_t nodesNumber);
    Grid4(size_t width, size_t height);
private:
    void construct(size_t width, size_t height);
};




class Grid6 : public Graph<> {
public:
    explicit Grid6(size_t nodesNumber);
    Grid6(size_t width, size_t height);
private:
    void construct(size_t width, size_t height);
};




class Grid8 : public Graph<> {
public:
    explicit Grid8(size_t nodesNumber);
    Grid8(size_t width, size_t height);
private:
    void construct(size_t width, size_t height);
};




class Paley : public Graph<> {
public:
    explicit Paley(size_t nodesNumber);
};




class LPS : public Graph<> {
public:
    explicit LPS(size_t nodesNumber);
private:
    bool miller_rabin_pass(const size_t a, const size_t s,
                           const size_t d, const size_t nodesNumber);

    bool millerRabin(const size_t nodesNumber);
};




class Margulis : public Graph<> {
public:
    explicit Margulis(size_t nodesNumber);
};




// random weighted graph for usage as e. g. ESN reservoir, TODO implement

class Random : public Graph<> {
public:
    explicit Random(size_t nodesNumber, double sparsity);
};




}
}


#include "graph.cpp"

#endif // header guards
