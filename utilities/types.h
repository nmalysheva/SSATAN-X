//
// Created by Malysheva, Nadezhda on 26.01.21.
//

#ifndef ALGO_TYPES_H
#define ALGO_TYPES_H

#include <vector>
#include <lemon/list_graph.h>
#include "contact_network/Specie.h"

struct specieState
{
    int id;
    Specie sp;
    std::vector<int> contacts;
};
using NetworkStorage = std::vector<std::pair<double, std::vector<specieState>>>;
using Edge = lemon::ListGraph::Edge;
using Node = lemon::ListGraph::Node;

#endif //ALGO_TYPES_H
