//
// Created by Malysheva, Nadezhda on 2019-07-30.
//

#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <unistd.h>
#include "SSATANX.h"
#include "utilities/Utility.h"
#include "algorithms/AndersonTauLeap.h"

std::mt19937_64 SSATANX::generator{std::random_device{}()};

SSATANX::SSATANX()
{
    generator.seed(::time(nullptr) * getpid()); //to change the seed for every run
}

void SSATANX::execute(double tStart, double tEnd, ContactNetwork &contNetwork, NetworkStorage &nwStorage,
                   size_t &nRejections, size_t &nAcceptance, size_t &nThin)
{
    double time = tStart;

    nwStorage.emplace_back(time, contNetwork.getNetworkState());

    double lookAheadTime  =  0; //init look-ahead time
    double propUpperLimit = -1; //init upper limit for propensitie sum
    double networkLastUpdate = tStart;

    double proposedTime = -1;
    std::vector<std::pair<double, Edge>> propTransmit = contNetwork.getTransmissionRateSum();
    std::vector<std::pair<double, Node>> propDeath = contNetwork.getDeathRateSum();
    std::vector<std::pair<double, Node>> propDiagnos = contNetwork.getDiagnosisRateSum();

    std::unordered_map<std::string, double >propensities{
            {"transmission", propTransmit.back().first},
            {"diagnosis",propDiagnos.back().first},
            {"death", propDeath.back().first},
            {"birth", contNetwork.getBirthRateSum()},
    };

    while (time < tEnd)
    {
        //choose look-ahead time
        lookAheadTime = tEnd - time;
        propUpperLimit = getPropUpperLimit(lookAheadTime, contNetwork,
                                           propDiagnos.back().first,
                                           propDeath.back().first);
        if (propUpperLimit == 0)
        {
            time = tEnd;
            nwStorage.emplace_back(time, contNetwork.getNetworkState());
            break;
        }
        else
        {
            double r = sampleRandUni(generator);

            proposedTime = 1 / propUpperLimit * std::log(1/r);
            if (proposedTime > lookAheadTime)
            {
                nRejections ++;
                time += lookAheadTime;

                nwStorage.emplace_back(time, contNetwork.getNetworkState());
            }
            else
            {
                time += proposedTime;
                Anderson::AndersonTauLeap(networkLastUpdate, time, contNetwork, generator);
                networkLastUpdate = time;
                propTransmit = contNetwork.getTransmissionRateSum();
                propensities.at("transmission") = propTransmit.back().first;
                propDiagnos = contNetwork.getDiagnosisRateSum();
                propensities.at("diagnosis") = propDiagnos.back().first;
                propDeath = contNetwork.getDeathRateSum();
                propensities.at("death") = propDeath.back().first;
                propensities.at("birth") = contNetwork.getBirthRateSum();

                double propensitieSum = std::accumulate(propensities.begin(), propensities.end(), 0.0, [] (double value, const std::map<std::string, double>::value_type& p)
                { return value + p.second; });

                r = sampleRandUni(generator);

                double searchBound = propUpperLimit * r;

                if (propensitieSum >= searchBound)
                {
                    nAcceptance ++;
                    double pSum = 0;

                    for (auto &it: propensities)
                    {

                        if (pSum + it.second >= searchBound)
                        {
                            if (it.first == "diagnosis")
                            {
                                networkLastUpdate = time;
                            }
                            executeReaction(contNetwork, it.first, pSum, searchBound, time,
                                    propTransmit,propDiagnos,propDeath);

                            nwStorage.emplace_back(time, contNetwork.getNetworkState());

                            propTransmit = contNetwork.getTransmissionRateSum();
                            propensities.at("transmission") = propTransmit.back().first;

                            propDeath = contNetwork.getDeathRateSum();
                            propensities.at("death") = propDeath.back().first;

                            propDiagnos = contNetwork.getDiagnosisRateSum();
                            propensities.at("diagnosis") = propDiagnos.back().first;

                            propensities.at("birth") = contNetwork.getBirthRateSum();

                            break;
                        }
                        pSum += it.second;

                    }
                }
                else
                {
                   nThin++;
                }

            }
        }
    }
}

double  SSATANX::getPropUpperLimit_naive (ContactNetwork & contNetwork,
                                 double diagnosisUpperLimit, double deathUpperLimit) const
{
    return (contNetwork.countByState(Specie::I) * contNetwork.getTransmissionRateLimit() +
            contNetwork.countByState(Specie::D) * contNetwork.getTransmissionRateLimit() * 0.5) *
           contNetwork.countByState(Specie::S) + diagnosisUpperLimit +
            deathUpperLimit;
}

double  SSATANX::getPropUpperLimit (double lookAheadTime, ContactNetwork & contNetwork, double diagnosisUpperLimit, double deathUpperLimit) const
{
    //get estimation of the Max.contacts based on rates
    size_t numberOfInfected = contNetwork.countByState(Specie::I);
    size_t numberOfDiagnosed = contNetwork.countByState(Specie::D);
    size_t numberOfSusceptible = contNetwork.countByState(Specie::S);

    double maxContInfected = contNetwork.getMaxContactsLimitByState(Specie::I, lookAheadTime);
    maxContInfected  = std::min(maxContInfected, static_cast<double>(numberOfSusceptible * numberOfInfected));

    double maxContDiagnosed = contNetwork.getMaxContactsLimitByState(Specie::D, lookAheadTime);
    maxContDiagnosed  = std::min(maxContDiagnosed, static_cast<double>(numberOfSusceptible * numberOfDiagnosed));

    double maxContSusceptible = contNetwork.getMaxContactsLimitByState(Specie::S, lookAheadTime);

    maxContSusceptible = std::min(maxContSusceptible, static_cast<double>(numberOfSusceptible * (numberOfInfected + numberOfDiagnosed)));

    double limit1 = 0;

    if (maxContInfected + maxContDiagnosed <= maxContSusceptible)
    {
        limit1 = maxContInfected  * contNetwork.getTransmissionRateLimit() + maxContDiagnosed * contNetwork.getTransmissionRateLimit() * 0.5;
    }
    else if (maxContSusceptible <= maxContInfected)
    {
        limit1 = maxContSusceptible * contNetwork.getTransmissionRateLimit();
    }
    else
    {
        limit1 = maxContInfected * contNetwork.getTransmissionRateLimit() + (maxContSusceptible - maxContInfected) * contNetwork.getTransmissionRateLimit() * 0.5;
    }

    //get estimation of the Max. possible contacts based on number of species.

    size_t maxContEsteem = (numberOfInfected + numberOfDiagnosed) * numberOfSusceptible;

    double limit2 = numberOfInfected * numberOfSusceptible * contNetwork.getTransmissionRateLimit() +
             numberOfDiagnosed * numberOfSusceptible * contNetwork.getTransmissionRateLimit() * 0.5;

    double result = diagnosisUpperLimit + deathUpperLimit;
    if (std::min(maxContInfected + maxContDiagnosed, maxContSusceptible) < maxContEsteem)
    {
        result += limit1;
    }
    else
    {
        result += limit2;
    }

    return result;

}


void SSATANX::executeReaction(ContactNetwork & contNetwork, const std::string &reactId, double rStart,
                          double rBound, double time,
                          std::vector<std::pair<double, Edge>> &propTransmit,
                          std::vector<std::pair<double, Node>> &propDiagnos,
                          std::vector<std::pair<double, Node>> &propDeath)
{

    if (reactId == "transmission")
    {
        auto edgeIterator = std::lower_bound(propTransmit.begin(), propTransmit.end(), rBound - rStart, lambdaLess);
        contNetwork.executeTransmission(edgeIterator->second, time);
    }
    else if (reactId == "diagnosis")
    {
        auto nodeIterator = std::lower_bound(propDiagnos.begin(), propDiagnos.end(), rBound - rStart, lambdaLess);
        contNetwork.executeDiagnosis(nodeIterator->second, time);
    }

    else if (reactId  == "death" )
    {
        auto nodeIterator = std::lower_bound(propDeath.begin(), propDeath.end(), rBound - rStart, lambdaLess);
        contNetwork.executeDeath(nodeIterator->second);
    }
}


