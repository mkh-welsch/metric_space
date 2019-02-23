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
namespace mapping {
namespace SOM_details {


// -------
// old Graph code, TODO remove


class Graph {
        public:
                explicit Graph(size_t nodesNumber);
                ~Graph();

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

class Grid4 : public Graph {
        public:
                explicit Grid4(size_t nodesNumber);
                Grid4(size_t width, size_t height);
        private:
                void construct(size_t width, size_t height);
};

class Grid6 : public Graph {
public:
        explicit Grid6(size_t nodesNumber);
        Grid6(size_t width, size_t height);
private:
        void construct(size_t width, size_t height);
};

class Grid8 : public Graph {
        public:
                explicit Grid8(size_t nodesNumber);
                Grid8(size_t width, size_t height);
        private:
                void construct(size_t width, size_t height);
};

class LPS : public Graph {
        public:
                explicit LPS(size_t nodesNumber);
        private:
                bool miller_rabin_pass(const size_t a, const size_t s,
                                                        const size_t d, const size_t nodesNumber);

                bool millerRabin(const size_t nodesNumber);
};

class Paley : public Graph {
        public:
                explicit Paley(size_t nodesNumber);
};

class Margulis : public Graph {
        public:
                explicit Margulis(size_t nodesNumber);
};



// end of old Graph code
// -------





// Graph_blaze


template <typename WeightType = bool, bool isDense = false, bool isSymmetric = true> //  TODO implement isDense & reorder params
class Graph_blaze {

private:
//    typedef typename std::conditional<
//        isWeighted,
//        WeightType,
//        bool
//        >::type ElementType;

    static constexpr bool isWeighted = !std::is_same<WeightType, bool>::value;

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
    explicit Graph_blaze(size_t nodesNumber);
    Graph_blaze();
    ~Graph_blaze();

    size_t getNodesNumber();
    bool isValid();

    std::vector<std::vector<size_t>> getNeighbours(const size_t nodeIndex, const size_t maxDeep);

    std::vector<std::vector<size_t>> getNeighborsNew(const size_t nodeIndex, const size_t maxDeep);

//    typename std::enable_if<std::is_same<WeightType, bool>::value, std::vector<std::vector<size_t>>>::type
//    getNeighborsNew(const size_t nodeIndex, const size_t maxDeep);

//    typename std::enable_if<std::is_same<WeightType, bool>::value, std::vector<std::vector<size_t>>>::type
//    getNeighborsNew(const size_t nodeIndex, const size_t maxDeep);

    MatrixType get_matrix();

    void buildEdges(const std::vector<std::pair<size_t, size_t>> &edgesPairs);

protected:
    size_t nodesNumber;
    bool valid;
//    std::vector<std::vector<size_t>> edges;

    MatrixType m;

    size_t modularPow(const size_t base, const size_t exponent, const size_t modulus);

//    void buildEdges(const std::vector<std::pair<size_t, size_t>> &edgesPairs); // moved to public

//    void neighboursWalk(const size_t index, std::unordered_map<size_t, size_t> &indexes, size_t deep,
//                        const size_t maxDeep);
};





class Grid4_blaze : public Graph_blaze<> {
public:
    explicit Grid4_blaze(size_t nodesNumber);
    Grid4_blaze(size_t width, size_t height);
private:
    void construct(size_t width, size_t height);
};




class Grid6_blaze : public Graph_blaze<> {
public:
    explicit Grid6_blaze(size_t nodesNumber);
    Grid6_blaze(size_t width, size_t height);
private:
    void construct(size_t width, size_t height);
};




class Grid8_blaze : public Graph_blaze<> {
public:
    explicit Grid8_blaze(size_t nodesNumber);
    Grid8_blaze(size_t width, size_t height);
private:
    void construct(size_t width, size_t height);
};




class Paley_blaze : public Graph_blaze<> {
public:
    explicit Paley_blaze(size_t nodesNumber);
};




class LPS_blaze : public Graph_blaze<> {
public:
    explicit LPS_blaze(size_t nodesNumber);
private:
    bool miller_rabin_pass(const size_t a, const size_t s,
                           const size_t d, const size_t nodesNumber);

    bool millerRabin(const size_t nodesNumber);
};




class Margulis_blaze : public Graph_blaze<> {
public:
    explicit Margulis_blaze(size_t nodesNumber);
};




}
}
}


#include "graph.cpp"

#endif // header guards
