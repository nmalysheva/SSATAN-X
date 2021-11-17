//
// Created by Malysheva, Nadezhda on 2019-07-30.
//

#ifndef ALGO_NSA_H
#define ALGO_NSA_H

#include <random>
#include "contact_network/ContactNetwork.h"


class SSATANX
{
public:
    SSATANX();
    void execute(double tStart, double tEnd, ContactNetwork &contNetwork,
                 NetworkStorage &nwStorage, size_t &nRejections, size_t &nAcceptance, size_t &nThin);
    ~SSATANX() {};

private:

    double  getPropUpperLimit (double lookAheadTime, ContactNetwork & contNetwork,
                               double diagnosisUpperLimit, double deathUpperLimit) const;

    double  getPropUpperLimit_naive (ContactNetwork & contNetwork,
                               double diagnosisUpperLimit, double deathUpperLimit) const;

    void executeReaction(ContactNetwork & contNetwork, const std::string &reactId, double rStart,
                         double rBound, double time,
                         std::vector<std::pair<double, Edge>> &propTransmit,
                         std::vector<std::pair<double, Node>> &propDiagnos,
                         std::vector<std::pair<double, Node>> &propDeath);

private:

    static std::mt19937_64 generator;
};


#endif //ALGO_NSA_H
