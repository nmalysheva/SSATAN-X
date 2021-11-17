/* Created by Malysheva, Nadezhda on 2019-07-28.
 *
 *
 * Gillespie algorithm (SSA)
 */


#ifndef ALGO_SSA_H
#define ALGO_SSA_H


#include <random>
#include <vector>

#include "contact_network/ContactNetwork.h"
#include "utilities/types.h"

class SSA
{
public:
    SSA();
    void execute(double tStart, double tEnd, ContactNetwork &contNetwork,
                 NetworkStorage &nwStorage);
    ~SSA() = default;

private:

    void   executeReaction(ContactNetwork & contNetwork, const std::string &reactId, double rStart,
                           double rBound, double time,
                            std::vector<std::pair<double, lemon::ListGraph::Edge>> &propDel,
                            std::vector<std::pair<double, lemon::ListGraph::Edge>> &propAdd,
                            std::vector<std::pair<double, lemon::ListGraph::Edge>> &propTransmit,
                            std::vector<std::pair<double, Node>> &propDiagnos,
                            std::vector<std::pair<double, Node>> &propDeath,
                           NetworkStorage &nwStorage);

private:
    static std::mt19937_64 generator;
};


#endif //ALGO_SSA_H
