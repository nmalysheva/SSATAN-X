//
// Created by Malysheva, Nadezhda on 2019-07-28.
//

#include <string>
#include <numeric>
#include <unistd.h>
#include <chrono>
#include "SSA.h"
#include "utilities/Utility.h"

std::mt19937_64 SSA::generator{std::random_device{}()};

SSA::SSA()
{
    generator.seed(::time(nullptr) * getpid()); //to change the seed for every run
}

void SSA::execute(double tStart, double tEnd, ContactNetwork &contNetwork,
                  NetworkStorage &nwStorage)
{
    std::vector<std::pair<double, Edge>> propDel;
    std::vector<std::pair<double, Edge>> propAdd;
    std::vector<std::pair<double, Edge>> propTransmit;
    std::vector<std::pair<double, Node>> propDiagnos;
    std::vector<std::pair<double, Node>> propDeath;

    std::unordered_map<std::string, double >propensities {
            {"edge_del", 0},
            {"edge_add", 0},
            {"transmission", 0},
            {"diagnosis", 0},
            {"death", 0},
            {"birth", 0},
    };


    double time = tStart;

    nwStorage.emplace_back(time, contNetwork.getNetworkState());

    while (time < tEnd)
    {
        propDel = contNetwork.getEdgeDeletionRateSum();
        propAdd = contNetwork.getEdgeAdditionRateSum();
        propTransmit = contNetwork.getTransmissionRateSum();
        propDiagnos = contNetwork.getDiagnosisRateSum();
        propDeath = contNetwork.getDeathRateSum();

        propensities.at("edge_del") = propDel.back().first;
        propensities.at("edge_add") = propAdd.back().first;

        propensities.at("transmission") = propTransmit.back().first;
        propensities.at("diagnosis") = propDiagnos.back().first;
        propensities.at("death") = propDeath.back().first;

        propensities.at("birth") = contNetwork.getBirthRateSum();


        double propensitieSum = std::accumulate(propensities.begin(), propensities.end(), 0.0, [](double value, const std::unordered_map<std::string, double >::value_type &b){return value + b.second;});

        if (propensitieSum == 0)
        {
            time = tEnd;
            nwStorage.emplace_back(time, contNetwork.getNetworkState());
            break;
        }
        double r = sampleRandUni(generator);
        double proposedTime = 1 / propensitieSum * std::log(1/r);
        if (time + proposedTime > tEnd )
        {
            time = tEnd;
            nwStorage.emplace_back(time, contNetwork.getNetworkState());
            break;
        }
        else
        {
            time += proposedTime;
            r = sampleRandUni(generator);
            double pSum = 0;
            for (auto const &it: propensities)
            {

                if (pSum + it.second >= propensitieSum * r)
                {
                    executeReaction(contNetwork, it.first, pSum, propensitieSum * r, time,
                                    propDel, propAdd, propTransmit, propDiagnos, propDeath, nwStorage);
                    break;
                }
                pSum += it.second;
            }
        }
    }
}

void SSA::executeReaction(ContactNetwork & contNetwork, const std::string &reactId, double rStart,
                          double rBound, double time,
                          std::vector<std::pair<double, Edge>> &propDel,
                          std::vector<std::pair<double, Edge>> &propAdd,
                          std::vector<std::pair<double, Edge>> &propTransmit,
                          std::vector<std::pair<double, Node>> &propDiagnos,
                          std::vector<std::pair<double, Node>> &propDeath,
                          NetworkStorage &nwStorage)
{
    if (reactId == "edge_del")
    {
        auto edgeIterator = std::lower_bound(propDel.begin(), propDel.end(), rBound - rStart, lambdaLess);
        contNetwork.removeEdge(edgeIterator->second);
    }

    else if (reactId == "edge_add")
    {
        auto edgeIterator = std::lower_bound(propAdd.begin(), propAdd.end(), rBound - rStart, lambdaLess);
        contNetwork.addEdge(edgeIterator->second);
    }
    else if (reactId == "transmission")
    {
        auto edgeIterator = std::lower_bound(propTransmit.begin(), propTransmit.end(), rBound - rStart, lambdaLess);
        contNetwork.executeTransmission(edgeIterator->second, time);

        nwStorage.emplace_back(time, contNetwork.getNetworkState());
    }

    else if (reactId == "diagnosis")
    {
        auto nodeIterator = std::lower_bound(propDiagnos.begin(), propDiagnos.end(), rBound - rStart, lambdaLess);
        contNetwork.executeDiagnosis(nodeIterator->second, time);
        nwStorage.emplace_back(time, contNetwork.getNetworkState());
    }

    else if (reactId == "death")
    {
        auto nodeIterator = std::lower_bound(propDeath.begin(), propDeath.end(), rBound - rStart, lambdaLess);
        contNetwork.executeDeath(nodeIterator->second);
        nwStorage.emplace_back(time, contNetwork.getNetworkState());
    }

    else if (reactId == "birth")
    {
        //contNetwork.executeBirth(rStart, rBound);
    }
}