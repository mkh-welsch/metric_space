#include <iostream>
#include "details/graph.hpp"



int main()
{



    size_t w = 3;
    size_t h = 3;
    size_t node = 1;
    size_t max_depth = 4;

    //*
    std::cout << "Testing Graph_blaze\n";


//    auto g = metric::graph::Grid4(h, w);
//    auto g = metric::graph::Grid6(h, w);
    auto g = metric::graph::Grid8(h, w);
//    auto g = metric::graph::Margulis(h*w);
//    auto g = metric::graph::Paley(h*w);
//    auto g = metric::graph::LPS(h*w); // TODO FIX: no nodes in graph
    std::cout << g.get_matrix() << "\n";

    std::cout << "new:\n";

    std::vector<std::vector<size_t>> neighbors = g.getNeighboursOld(node, max_depth);
    for (size_t i=0; i<neighbors.size(); i++)
        for (size_t j=0; j<neighbors[i].size(); j++)
            std::cout << i << " | " << neighbors[i][j] << "\n";

    std::cout << "new modified:\n";

    std::vector<std::vector<size_t>> neighborsNew = g.getNeighbours(node, max_depth);
    for (size_t i=0; i<neighborsNew.size(); i++)
        for (size_t j=0; j<neighborsNew[i].size(); j++)
            std::cout << i << " | " << neighborsNew[i][j] << "\n";


//    auto g_old = metric::graph::Grid4_old(h, w);
//    auto g_old = metric::graph::Grid6_old(h, w);
    auto g_old = metric::graph::Grid8_old(h, w);
//    auto g_old = metric::graph::Margulis_old(h*w);
//    auto g_old = metric::graph::Paley_old(h*w);
//    auto g_old = metric::graph::LPS_old(h*w); // TODO FIX: no nodes in graph

    std::cout << "old:" << "\n";

    std::vector<std::vector<size_t>> neighbors_old = g_old.getNeighbours(node, max_depth);
    for (size_t i=0; i<neighbors_old.size(); i++)
        for (size_t j=0; j<neighbors_old[i].size(); j++)
            std::cout << i << " | " << neighbors_old[i][j] << "\n";

    //*/

    // testing template parameters

    std::vector<std::pair<size_t, size_t>> edges;
    edges.emplace_back(0, 1);
    edges.emplace_back(1, 2);
    edges.emplace_back(2, 0);

    auto g_custom = metric::graph::Graph<bool, true, false>(edges); // edge value type, isDense, isSymmetric
//    g_custom.buildEdges(edges);

    std::cout << "\ncustom graph:\n";
    std::cout << g_custom.get_matrix() << "\n";

    std::vector<std::vector<size_t>> neighborsNewCustom = g_custom.getNeighbours(node, max_depth);
    for (size_t i=0; i<neighborsNewCustom.size(); i++)
        for (size_t j=0; j<neighborsNewCustom[i].size(); j++)
            std::cout << i << " | " << neighborsNewCustom[i][j] << "\n";

    return 0;

}
