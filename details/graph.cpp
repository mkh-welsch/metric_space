#include "graph.hpp"
//#include "graph/connected_components.hpp"

#include <unordered_set> // for old Graph code, TODO remove


namespace metric {
namespace graph {


// -------
// old Graph code, TODO remove



Graph_old::Graph_old(size_t nodesNumber) : nodesNumber(nodesNumber), valid(false) {};

Graph_old::~Graph_old() = default;

size_t Graph_old::getNodesNumber()
{
        return nodesNumber;
}

bool Graph_old::isValid() {
        return valid;
}

std::vector<std::vector<size_t>> Graph_old::getNeighbours(const size_t nodeIndex, const size_t maxDeep)
{
        std::vector<std::vector<size_t>> neighboursList(maxDeep + 1);

        std::unordered_map<size_t, size_t> indexes;
        neighboursWalk(nodeIndex, indexes, 0, maxDeep);

        for (const auto& [index, deep]: indexes) {
                neighboursList[deep].push_back(index);
        }

        return neighboursList;
}

size_t Graph_old::modularPow(const size_t base, const size_t exponent, const size_t modulus)
{
        if (modulus == 1) {
                return 1;
        }
        size_t c = 1;
        for (size_t e = 0; e < exponent; ++e) {
                c = (c * base) % modulus;
        }

        return c;
}

void Graph_old::buildEdges(const std::vector<std::pair<size_t, size_t>> &edgesPairs)
{
        std::vector<std::unordered_set<size_t>> edgesSets(nodesNumber);

        for (const auto& [i, j]: edgesPairs) {
                if (i != j) {
                        edgesSets[i].insert(j);
                        edgesSets[j].insert(i);
                }
        }

        edges.resize(nodesNumber);

for (auto i = 0; i < edgesSets.size(); ++i) {
                for (const auto& j: edgesSets[i]) {
                        edges[i].push_back(j);
                }
        }
}

void Graph_old::neighboursWalk(const size_t index, std::unordered_map<size_t, size_t> &indexes, size_t deep,
                                                  const size_t maxDeep)
{
        if (deep > maxDeep) {
                return;
        }

        auto iterator = indexes.find(index);
        if (iterator != indexes.end()) {
                if (iterator->second <= deep) {
                        return;
                }
        }

        indexes[index] = deep;

        // code for debug, TODO remove
        // std::cout << "added node: " << index << "\n";
        // std::cout << "depth: " << deep << "\n";
        // end of debug code

        for (const auto& i: edges[index]) {
                neighboursWalk(i, indexes, deep + 1, maxDeep);
        }
}

Grid4_old::Grid4_old(size_t nodesNumber) : Graph_old(nodesNumber)
{
        int s = sqrt(nodesNumber);
        if ((s * s) != nodesNumber) {
                valid = false;
        } else {
                construct(s, s);
        }
}

Grid4_old::Grid4_old(size_t width, size_t height) : Graph_old(width * height)
{
        construct(width, height);
}

void Grid4_old::construct(size_t width, size_t height)
{
        edges.resize(width * height);

        for (size_t i = 0; i < height; ++i) {
                for (size_t j = 0; j < width; ++j) {

                        int ii0 = 0, ii1 = 0;
                        int jj0 = 0, jj1 = 0;

                        if (i > 0) {
                                ii0 = -1;
                        }
                        if (j > 0) {
                                jj0 = -1;
                        }

                        if (i < height - 1) {
                                ii1 = 1;
                        }
                        if (j < width - 1) {
                                jj1 = 1;
                        }

                        for (int ii = ii0; ii <= ii1; ++ii) {
                                for (int jj = jj0; jj <= jj1; ++jj) {
                                        if ((ii == 0) or (jj == 0)) {
                                                edges[i * width + j].push_back((i + ii) * width + (j + jj));
                                        }
                                }
                        }
                }
        }

        valid = true;
}

Grid6_old::Grid6_old(size_t nodesNumber) : Graph_old(nodesNumber)
{
        int s = sqrt(nodesNumber);
        if ((s * s) != nodesNumber) {
                valid = false;
        } else {
                construct(s, s);
        }
}

Grid6_old::Grid6_old(size_t width, size_t height) : Graph_old(width * height)
{
        construct(width, height);
}

void Grid6_old::construct(size_t width, size_t height)
{
        edges.resize(width * height);

        for (size_t i = 0; i < height; ++i) {
                for (size_t j = 0; j < width; ++j) {

                        bool odd = true;
                        if (i % 2 == 0) {
    odd = false;
                        }

                        bool up = true;
                        if (i == 0) {
                                up = false;
                        }
                        bool down = true;
                        if (i == height - 1) {
                                down = false;
                        }
                        bool left = true;
                        if (j == 0) {
                                left = false;
                        }

                        bool right = true;
if (j == width - 1) {
        right = false;
}

                        if (up) {
                                edges[i * width + j].push_back((i - 1) * width + (j + 0));
                        }
                        if (down) {
                                edges[i * width + j].push_back((i + 1) * width + (j + 0));
                        }
                        if (left) {
                                edges[i * width + j].push_back((i + 0) * width + (j - 1));
                        }
                        if (right) {
                                edges[i * width + j].push_back((i + 0) * width + (j + 1));
                        }

                        if (!odd and left) {
                                if (up) {
                                        edges[i * width + j].push_back((i - 1) * width + (j - 1));
                                }
                                if (down) {
                                        edges[i * width + j].push_back((i + 1) * width + (j - 1));
                                }
                        }

                        if (odd and right) {
                                if (up) {
                                        edges[i * width + j].push_back((i - 1) * width + (j + 1));
                                }
                                if (down)
                                        edges[i * width + j].push_back((i + 1) * width + (j + 1));
                        }
                }
        }

        valid = true;
}

Grid8_old::Grid8_old(size_t nodesNumber) : Graph_old(nodesNumber)
{
        int s = sqrt(nodesNumber);
        if ((s * s) != nodesNumber) {
                valid = false;
        } else {
                construct(s, s);
        }
}

Grid8_old::Grid8_old(size_t width, size_t height) : Graph_old(width * height)
{
        construct(width, height);
}

void Grid8_old::construct(size_t width, size_t height)
{
        edges.resize(width * height);

        for (size_t i = 0; i < height; ++i) {
                for (size_t j = 0; j < width; ++j) {

                        int ii0 = 0, ii1 = 0;
                        int jj0 = 0, jj1 = 0;

                        if (i > 0) {
                                ii0 = -1;
                        }
                        if (j > 0) {
                                jj0 = -1;
                        }

                        if (i < height - 1) {
                                ii1 = 1;
                        }
                        if (j < width - 1) {
                                jj1 = 1;
                        }

                        for (int ii = ii0; ii <= ii1; ++ii) {
                                for (int jj = jj0; jj <= jj1; ++jj) {
                                        if ((ii != 0) or (jj != 0)) {
                                                edges[i * width + j].push_back((i + ii) * width + (j + jj));
                                        }
                                }
                        }
                }
        }

        valid = true;
}

LPS_old::LPS_old(size_t nodesNumber) : Graph_old(nodesNumber)
{
        if (!millerRabin(nodesNumber)) {
                return;
        }

        std::vector<std::pair<size_t, size_t>> edgesPairs;

        for (size_t i = 0; i < nodesNumber; ++i) {
                if (i == 0) {
                        edgesPairs.emplace_back(0, nodesNumber - 1);
                        edgesPairs.emplace_back(0, 1);
                } else {
                        edgesPairs.emplace_back(i, i - 1);
                        edgesPairs.emplace_back(i, (i + 1) % nodesNumber);
                        edgesPairs.emplace_back(i, modularPow(i, nodesNumber - 2, nodesNumber));
                }
        }

        buildEdges(edgesPairs);

        valid = true;
}

bool LPS_old::millerRabin(const size_t nodesNumber)
{
        srand(time(NULL));

        auto d = nodesNumber - 1;
        auto s = 0;

        while (d % 2 == 0) {
                d >>= 1;
                s += 1;
        }

        for (int repeat = 0; repeat < 20; ++repeat) {
                size_t a = 0;
                while (a == 0) {
                        a = rand() % nodesNumber;
                }

                if (!miller_rabin_pass(a, s, d, nodesNumber)) {
                        return false;
                }
        }
        return true;
}

bool LPS_old::miller_rabin_pass(const size_t a, const size_t s,
                                                        const size_t d, const size_t nodesNumber)
{
        auto p = size_t(pow(a, d)) % nodesNumber;
        if (p == 1) {
                return true;
        }

        for (auto i = 0; i < s - 1; ++i) {
                if (p == nodesNumber - 1) {
                        return true;
                }
                p = (p * p) % nodesNumber;
        }

        return p == nodesNumber - 1;
}

Paley_old::Paley_old(size_t nodesNumber) : Graph_old(nodesNumber)
{
        if (nodesNumber % 4 != 1) {
                return;
        }

        std::vector<std::pair<size_t, size_t>> edgesPairs;

        std::vector<size_t> squareList;

        size_t l = (nodesNumber - 1) / 2;
        squareList.reserve(l);

        for (auto i = 0; i < l; ++i) {
                squareList.push_back(i * i % nodesNumber);
        }

        for (auto i = 0; i < nodesNumber; ++i) {
                for (auto j: squareList) {
                        edgesPairs.emplace_back(i, (i + j) % nodesNumber);
                }
        }

        buildEdges(edgesPairs);

        valid = true;
}

Margulis_old::Margulis_old(size_t nodesNumber) : Graph_old(nodesNumber) {
        int s = sqrt(nodesNumber);
        if ((s * s) != nodesNumber) {
                valid = false;
        } else {
                std::vector<std::pair<size_t, size_t>> edgesPairs;

                for (size_t i = 0; i < s; ++i) {
                        for (size_t j = 0; j < s; ++j) {
                                edgesPairs.emplace_back(i * s + j, ((i + 2 * j) % s) * s + j);
                                edgesPairs.emplace_back(i * s + j, i * s + (2 * i + j) % s);
                                edgesPairs.emplace_back(i * s + j, i * s + (2 * i + j + 1) % s);
                                edgesPairs.emplace_back(i * s + j, ((i + 2 * j + 1) % s) * s + j);
                        }
                }

                buildEdges(edgesPairs);

                valid = true;
        }
}





// end of old Graph code
// -------





// Graph based on blaze-lib



template <typename WeightType, bool isDense, bool isSymmetric>
Graph<WeightType, isDense, isSymmetric>::Graph(size_t nodesNumber) : nodesNumber(nodesNumber), valid(false) {}

template <typename WeightType, bool isDense, bool isSymmetric>
Graph<WeightType, isDense, isSymmetric>::Graph() : nodesNumber(0), valid(false) {}

template <typename WeightType, bool isDense, bool isSymmetric>
Graph<WeightType, isDense, isSymmetric>::Graph(const std::vector<std::pair<size_t, size_t>> & edgesPairs)
{
    buildEdges(edgesPairs);
}

template <typename WeightType, bool isDense, bool isSymmetric>
Graph<WeightType, isDense, isSymmetric>::~Graph() = default;

template <typename WeightType, bool isDense, bool isSymmetric>
size_t Graph<WeightType, isDense, isSymmetric>::getNodesNumber()
{
    return nodesNumber;
}

template <typename WeightType, bool isDense, bool isSymmetric>
bool Graph<WeightType, isDense, isSymmetric>::isValid() {
    return valid;
}

//* // left for comparison
template <typename WeightType, bool isDense, bool isSymmetric>
std::vector<std::vector<size_t>>
Graph<WeightType, isDense, isSymmetric>::getNeighboursOld(const size_t index, const size_t maxDeep)
{

    if (isWeighted)
        return std::vector<std::vector<size_t>>(0); // return empty if weighted, TODO implement weight-based metric if needed

    std::vector<std::vector<size_t>> neighboursList(maxDeep+1);


    std::stack<typename Graph<WeightType, isDense, isSymmetric>::MatrixType::Iterator> iterator_stack;
    std::stack<size_t> row_stack;
    std::unordered_map<size_t, size_t> indices;

    row_stack.push(index);
    iterator_stack.push(m.begin(index)); // stacks are ever processed together (so they are always of the same size)

    indices[row_stack.top()] = 0;

    size_t depth = 1;

    while (true)
    {
        if ( iterator_stack.top() == m.end(row_stack.top()) || depth > maxDeep ) // end of row or max depth reached
        {
            row_stack.pop(); // return to the previous level
            iterator_stack.pop();
            depth--;

            if (depth < 1) // finish
                break;
        }
        else
        {
            row_stack.push(iterator_stack.top()->index()); // enter the next level
            iterator_stack.push(m.begin(row_stack.top()));
            depth++;

            auto search = indices.find(row_stack.top());
            if (search == indices.end() || search->second > depth - 1) // node not exists or its depth is greater than current
                indices[row_stack.top()] = depth - 1; // write node into output
            //std::cout << "checked node: " << row_stack.top() << "\n"; // output for DEBUG, TODO remove
            //std::cout << "depth: " << depth - 1 << "\n";

            continue; // prevent from a step along the level when entered the new level
        }
        iterator_stack.top()++; // step along the level
    }

    for (const auto& [index, deep]: indices) {
        neighboursList[deep].push_back(index);
    }

    return neighboursList;
}

//*/

template <typename WeightType, bool isDense, bool isSymmetric>
size_t Graph<WeightType, isDense, isSymmetric>::modularPow(const size_t base, const size_t exponent, const size_t modulus)
{
    if (modulus == 1) {
        return 1;
    }
    size_t c = 1;
    for (size_t e = 0; e < exponent; ++e) {
        c = (c * base) % modulus;
    }

    return c;
}


template <typename WeightType, bool isDense, bool isSymmetric>
void Graph<WeightType, isDense, isSymmetric>::buildEdges(const std::vector<std::pair<size_t, size_t>> &edgesPairs)
{
    size_t max = 0;
    for (const auto & [i, j]: edgesPairs) { // TODO loop may be removed if max is passed as a parameter
        if (i > max)
            max = i;
        if (j > max)
            max = j;
    }

    max = max + 1;

    if (max>nodesNumber)
        nodesNumber = max;

    m.resize((unsigned long)max, (unsigned long)max);
    m.reset();

    //m.reserve(edgePairs.size()); // TODO optimize via reserve-insert-finalize idiom
    for (const auto & [i, j]: edgesPairs) {
        // std::cout << "m:\n" << m << "\n";
        if (i != j)
        {
            m(i, j) = 1;
            // std::cout << "adding to matrix: i=" << i << ", j=" << j << "\n";
        }

        /* // old code. TODO on the final optimization stage: restore faster insertion via .insert(...) with compile-time check/specialization
        {
            // edgesSets[i].insert(j);
            // edgesSets[j].insert(i);

//            if (directed)
            if (!isSymmetric)
            {
                //if (i < max && j < max)
                    m.insert(i, j, 1); // no need in check whether element exists (??)
            }
            else
            {
                //m.set(i, j, 1);
                //m.set(j, i, 1);
//                if (m.find(i, j) == m.end(i))
//                    m.insert(i, j, 1);
//                if (m.find(j, i) == m.end(j))
//                    m.insert(j, i, 1);

                // code for debugging, TODO remove
//                if (i > 7 || j > 7)
//                    std::cout << "exceeded\n";
                // end of code for debugging

                //if (i < max && j < max)
                {
                    m(i, j) = 1;
//                    m(j, i) = 1;
                }

            }

        }
        */
    }

//    edges.resize(nodesNumber);

//    for (auto i = 0; i < edgesSets.size(); ++i) {
//        for (const auto& j: edgesSets[i]) {
//            edges[i].push_back(j);
//        }
//    }
}




template <typename WeightType, bool isDense, bool isSymmetric>
template <typename T, bool denseFlag>
typename std::enable_if_t<!std::is_same<T, bool>::value, std::vector<std::vector<size_t>>>
Graph<WeightType, isDense, isSymmetric>::getNeighbours(const size_t index, const size_t maxDeep)
{
    //std::cout << "non-default value type specialization called\n";
    return std::vector<std::vector<size_t>>(0); // return empty if weighted, TODO implement weight-based metric if needed
}


template <typename WeightType, bool isDense, bool isSymmetric>
template <typename T, bool denseFlag>
typename std::enable_if_t<std::is_same<T, bool>::value && !denseFlag, std::vector<std::vector<size_t>>>
Graph<WeightType, isDense, isSymmetric>::getNeighbours(const size_t index, const size_t maxDeep)
{
    std::cout << "Sparse & default value type specialization called\n"; // TODO remove

    std::vector<std::vector<size_t>> neighboursList(maxDeep+1);

    if (index >= m.columns())
        return neighboursList;

    std::vector<size_t> parents;
    std::vector<size_t> children;
    std::vector<bool> nodes = std::vector<bool>(m.columns(), false);

    parents.push_back(index);
    neighboursList[0].push_back(index);
    nodes[index] = true;

    size_t depth = 1;

    while (depth<=maxDeep)
    {
        typename MatrixType::Iterator it;

        for (auto el : parents)
        {
            for (it = m.begin(el); it!=m.end(el); it++) // for dense and for sparse m.end(...) has different meaning!..
            {
                //* // tested on sparse only
                if (!nodes[it->index()])
                    neighboursList[depth].push_back(it->index()); // write node into output
                children.push_back(it->index());
                nodes[it->index()] = true;
                //std::cout << "added node: " << it->index() << "\n"; // output for DEBUG, TODO remove
                //std::cout << "depth: " << depth << "\n";
                //*/

                /* // code for dense, not appropriate for sparse due to m.end(...)
                size_t idx = it - m.begin(el);
                if (m(el, idx) != 1)
                    continue;
                if (!nodes[idx])
                    neighboursList[depth].push_back(idx); // write node into output
                children.push_back(idx);
                nodes[idx] = true;
                std::cout << "added node: " << idx << "\n"; // output for DEBUG, TODO remove
                std::cout << "parent: " << el << "\n";
                std::cout << "depth: " << depth << "\n";
                std::cout << "value: " << m(el, idx) << "\n";
                std::cout << "last index in row: " << m.end(el) - m.begin(el) << "\n";

                //*/
            }
        }

        depth++;
        parents.swap(children);
        children = {};
    }

    return neighboursList;
}


template <typename WeightType, bool isDense, bool isSymmetric>
template <typename T, bool denseFlag>
typename std::enable_if_t<std::is_same<T, bool>::value && denseFlag, std::vector<std::vector<size_t>>>
Graph<WeightType, isDense, isSymmetric>::getNeighbours(const size_t index, const size_t maxDeep)
{
        // similar to sparse specialization except the way m elements are accessed

        //std::cout << "isDense & default value type specialization called\n";

        std::vector<std::vector<size_t>> neighboursList(maxDeep+1);

        if (index >= m.columns())
            return neighboursList;

        std::vector<size_t> parents;
        std::vector<size_t> children;
        std::vector<bool> nodes = std::vector<bool>(m.columns(), false);

        parents.push_back(index);
        neighboursList[0].push_back(index);
        nodes[index] = true;

        size_t depth = 1;

        while (depth<=maxDeep)
        {
            typename MatrixType::Iterator it;

            for (auto el : parents)
            {
                for (it = m.begin(el); it!=m.end(el); it++)
                {
                    size_t idx = it - m.begin(el);
                    if (m(el, idx) != 1)
                        continue;
                    if (!nodes[idx])
                        neighboursList[depth].push_back(idx); // write node into output
                    children.push_back(idx);
                    nodes[idx] = true;
                    //std::cout << "added node: " << idx << "\n"; // output for DEBUG, TODO remove
                    //std::cout << "depth: " << depth << "\n";
                }
            }

            depth++;
            parents.swap(children);
            children = {};
        }

        return neighboursList;
}



template <typename WeightType, bool isDense, bool isSymmetric>
typename Graph<WeightType, isDense, isSymmetric>::MatrixType
Graph<WeightType, isDense, isSymmetric>::get_matrix()
{
    return m;
}

// end of base class implementation






// Grid4_blaze

Grid4::Grid4(size_t nodesNumber) : Graph<>(nodesNumber)
{
    int s = sqrt(nodesNumber);
    if ((s * s) != nodesNumber) {
        valid = false;
    } else {
        construct(s, s);
    }
}


Grid4::Grid4(size_t width, size_t height) : Graph<>(width * height)
{
    construct(width, height);
}

void Grid4::construct(size_t width, size_t height)
{
    unsigned long n_nodes = width*height;
    m.resize(n_nodes, n_nodes);

    std::vector<std::pair<size_t, size_t>> edgesPairs;

    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {

            int ii0 = 0, ii1 = 0;
            int jj0 = 0, jj1 = 0;

            if (i > 0) {
                ii0 = -1;
            }
            if (j > 0) {
                jj0 = -1;
            }

            if (i < height - 1) {
                ii1 = 1;
            }
            if (j < width - 1) {
                jj1 = 1;
            }

            for (int ii = ii0; ii <= ii1; ++ii) {
                for (int jj = jj0; jj <= jj1; ++jj) {
                    if ((ii == 0) or (jj == 0)) {
                        // edges[i * width + j].push_back((i + ii) * width + (j + jj));
                        //m.insert(i * width + j, (i + ii) * width + (j + jj), 1);
                        edgesPairs.emplace_back(i * width + j, (i + ii) * width + (j + jj));
                    }
                }
            }
        }
    }

    buildEdges(edgesPairs);

    valid = true;
}






// Grig6_blaze

Grid6::Grid6(size_t nodesNumber) : Graph<>(nodesNumber)
{
    int s = sqrt(nodesNumber);
    if ((s * s) != nodesNumber) {
        valid = false;
    } else {
        construct(s, s);
    }
}

Grid6::Grid6(size_t width, size_t height) : Graph<>(width * height)
{
    construct(width, height);
}

void Grid6::construct(size_t width, size_t height)
{
    unsigned long n_nodes = width*height;
    m.resize(n_nodes, n_nodes);

    std::vector<std::pair<size_t, size_t>> edgesPairs;

    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {

            bool odd = true;
            if (i % 2 == 0) {
                odd = false;
            }

            bool up = true;
            if (i == 0) {
                up = false;
            }
            bool down = true;
            if (i == height - 1) {
                down = false;
            }
            bool left = true;
            if (j == 0) {
                left = false;
            }

            bool right = true;
            if (j == width - 1) {
                right = false;
            }

            if (up) {
                //edges[i * width + j].push_back((i - 1) * width + (j + 0));
                //m.insert(i * width + j, (i - 1) * width + (j + 0), 1);
                edgesPairs.emplace_back(i * width + j, (i - 1) * width + (j + 0));
            }
            if (down) {
                //edges[i * width + j].push_back((i + 1) * width + (j + 0));
                //m.insert(i * width + j, (i + 1) * width + (j + 0), 1);
                edgesPairs.emplace_back(i * width + j, (i + 1) * width + (j + 0));
            }
            if (left) {
                //edges[i * width + j].push_back((i + 0) * width + (j - 1));
                //m.insert(i * width + j, (i + 0) * width + (j - 1), 1);
                edgesPairs.emplace_back(i * width + j, (i + 0) * width + (j - 1));
            }
            if (right) {
                //edges[i * width + j].push_back((i + 0) * width + (j + 1));
                //m.insert(i * width + j, (i + 0) * width + (j + 1), 1);
                edgesPairs.emplace_back(i * width + j, (i + 0) * width + (j + 1));
            }

            if (!odd and left) {
                if (up) {
                    //edges[i * width + j].push_back((i - 1) * width + (j - 1));
                    //m.insert(i * width + j, (i - 1) * width + (j - 1), 1);
                    edgesPairs.emplace_back(i * width + j, (i - 1) * width + (j - 1));
                }
                if (down) {
                    //edges[i * width + j].push_back((i + 1) * width + (j - 1));
                    //m.insert(i * width + j, (i + 1) * width + (j - 1), 1);
                    edgesPairs.emplace_back(i * width + j, (i + 1) * width + (j - 1));
                }
            }

            if (odd and right) {
                if (up) {
                    //edges[i * width + j].push_back((i - 1) * width + (j + 1));
                    //m.insert(i * width + j, (i - 1) * width + (j + 1), 1);
                    edgesPairs.emplace_back(i * width + j, (i - 1) * width + (j + 1));
                }
                if (down)
                    //edges[i * width + j].push_back((i + 1) * width + (j + 1));
                    //m.insert(i * width + j, (i + 1) * width + (j + 1), 1);
                    edgesPairs.emplace_back(i * width + j, (i + 1) * width + (j + 1));
            }
        }
    }

    valid = true;

    buildEdges(edgesPairs);
}





// Grid8_blaze

Grid8::Grid8(size_t nodesNumber) : Graph<>(nodesNumber)
{
    int s = sqrt(nodesNumber);
    if ((s * s) != nodesNumber) {
        valid = false;
    } else {
        construct(s, s);
    }
}

Grid8::Grid8(size_t width, size_t height) : Graph<>(width * height)
{
    construct(width, height);
}

void Grid8::construct(size_t width, size_t height)
{
    unsigned long n_nodes = width*height;
    m.resize(n_nodes, n_nodes);

    std::vector<std::pair<size_t, size_t>> edgesPairs;

    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {

            int ii0 = 0, ii1 = 0;
            int jj0 = 0, jj1 = 0;

            if (i > 0) {
                ii0 = -1;
            }
            if (j > 0) {
                jj0 = -1;
            }

            if (i < height - 1) {
                ii1 = 1;
            }
            if (j < width - 1) {
                jj1 = 1;
            }

            for (int ii = ii0; ii <= ii1; ++ii) {
                for (int jj = jj0; jj <= jj1; ++jj) {
                    if ((ii != 0) or (jj != 0)) {
//                        edges[i * width + j].push_back((i + ii) * width + (j + jj));
//                        m.insert(i * width + j, (i + ii) * width + (j + jj), 1); // TODO replace with polymorphyc metod suitable also for SymmetricMatrix!!
                        edgesPairs.emplace_back(i * width + j, (i + ii) * width + (j + jj));
                    }
                }
            }
        }
    }

    buildEdges(edgesPairs);

    valid = true;
}





// Paley_blaze

Paley::Paley(size_t nodesNumber) : Graph<>(nodesNumber)
{
    if (nodesNumber % 4 != 1) {
        return;
    }

    std::vector<std::pair<size_t, size_t>> edgesPairs;

    std::vector<size_t> squareList;

    size_t l = (nodesNumber - 1) / 2;
    squareList.reserve(l);

    for (auto i = 0; i < l; ++i) {
        squareList.push_back(i * i % nodesNumber);
    }

    for (auto i = 0; i < nodesNumber; ++i) {
        for (auto j: squareList) {
            edgesPairs.emplace_back(i, (i + j) % nodesNumber);
        }
    }

    buildEdges(edgesPairs);

    valid = true;
}






// LPS_blaze

LPS::LPS(size_t nodesNumber) : Graph<>(nodesNumber)
{
        if (!millerRabin(nodesNumber)) {
                return;
        }

        std::vector<std::pair<size_t, size_t>> edgesPairs;

        for (size_t i = 0; i < nodesNumber; ++i) {
                if (i == 0) {
                        edgesPairs.emplace_back(0, nodesNumber - 1);
                        edgesPairs.emplace_back(0, 1);
                } else {
                        edgesPairs.emplace_back(i, i - 1);
                        edgesPairs.emplace_back(i, (i + 1) % nodesNumber);
                        edgesPairs.emplace_back(i, modularPow(i, nodesNumber - 2, nodesNumber));
                }
        }

        buildEdges(edgesPairs);

        valid = true;
}

bool LPS::millerRabin(const size_t nodesNumber)
{
        srand(time(NULL));

        auto d = nodesNumber - 1;
        auto s = 0;

        while (d % 2 == 0) {
                d >>= 1;
                s += 1;
        }

        for (int repeat = 0; repeat < 20; ++repeat) {
                size_t a = 0;
                while (a == 0) {
                        a = rand() % nodesNumber;
                }

                if (!miller_rabin_pass(a, s, d, nodesNumber)) {
                        return false;
                }
        }
        return true;
}

bool LPS::miller_rabin_pass(const size_t a, const size_t s,
                                                        const size_t d, const size_t nodesNumber)
{
        auto p = size_t(pow(a, d)) % nodesNumber;
        if (p == 1) {
                return true;
        }

        for (auto i = 0; i < s - 1; ++i) {
                if (p == nodesNumber - 1) {
                        return true;
                }
                p = (p * p) % nodesNumber;
        }

        return p == nodesNumber - 1;
}






// Margulis_blaze

Margulis::Margulis(size_t nodesNumber) : Graph<>(nodesNumber) {
    int s = sqrt(nodesNumber);
    if ((s * s) != nodesNumber) {
        valid = false;
    } else {
        std::vector<std::pair<size_t, size_t>> edgesPairs;

        for (size_t i = 0; i < s; ++i) {
            for (size_t j = 0; j < s; ++j) {
                edgesPairs.emplace_back(i * s + j, ((i + 2 * j) % s) * s + j);
                edgesPairs.emplace_back(i * s + j, i * s + (2 * i + j) % s);
                edgesPairs.emplace_back(i * s + j, i * s + (2 * i + j + 1) % s);
                edgesPairs.emplace_back(i * s + j, ((i + 2 * j + 1) % s) * s + j);
            }
        }

        buildEdges(edgesPairs);

        valid = true;
    }
}






}
}
