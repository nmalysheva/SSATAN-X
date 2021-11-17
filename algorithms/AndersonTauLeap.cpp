//
// Created by Malysheva, Nadezhda on 21.07.20.
//
#include <numeric>
#include <chrono>
#include "AndersonTauLeap.h"
#include "utilities/Utility.h"
#include "utilities/types.h"


void Anderson::AndersonTauLeap(double &tLastNetworkUpdate, double tEnd, ContactNetwork & contNetwork,
                     std::mt19937_64 &generator)
{
    size_t M = 2; //number of reactions, size of propensity vector

    std::vector<double> T(M, 0);
    std::vector<size_t> C(M, 0);

    std::vector<std::vector<std::pair<double, size_t>>> S;
    S.reserve(1e6 + 1);
    S.push_back({{0.0, 0}});
    S.push_back({{0.0, 0}});

    std::vector<double> propensities(M, 0);
    std::vector<std::pair<double, Edge>> propDel = contNetwork.getEdgeDeletionRateSum();
    std::vector<std::pair<double, Edge>> propAdd = contNetwork.getEdgeAdditionRateSum();

    propensities.at(0) = propDel.back().first;
    propensities.at(1) = propAdd.back().first;

    std::vector<size_t> row(M, 0);

    double t = tLastNetworkUpdate;

    std::vector<size_t> X = {propDel.size() - 1, propAdd.size() - 1};
    double tau = getTau(propensities, X);
    while (t < tEnd)
    {
        if (propensities.at(0) + propensities.at(1) == 0)
        {
            t = tEnd;
            break;
        }
        tau = std::min(tau, tEnd - t);

        if (tau < 10.0 / (propensities.at(0) + propensities.at(1)))
        {
            executeSSA(100, M, tEnd, contNetwork, t, tLastNetworkUpdate, generator, T, C, S);
            propDel = contNetwork.getEdgeDeletionRateSum();
            propAdd = contNetwork.getEdgeAdditionRateSum();

            propensities.at(0) = propDel.back().first;
            propensities.at(1) = propAdd.back().first;
            X = {propDel.size() - 1, propAdd.size() - 1};
            tau = getTau(propensities, X);
        }

        else
        {
            std::vector<size_t> change = getChange(M, S, T, C, propensities, row, tau, generator);

            bool pass = change.at(0) <= std::max(epsilon * (propDel.size() - 1), 1.0) &&
                        change.at(1) <= std::max(epsilon * (propAdd.size() - 1), 1.0);

            if (pass)
            {
                if (t + tau > tEnd)
                {
                    t = tEnd;
                    break;
                }
                acceptLeap(t, tLastNetworkUpdate, tau, M, contNetwork,
                        S, T, C, row, propensities, change,
                        propAdd,propDel,generator);
            }
            else
            {
                rejectLeap(M,  tau, S, T, C, propensities, row, change);
            }
        }
    }
}

void Anderson::executeSSA(size_t n, size_t M, double tEnd, ContactNetwork & contNetwork, double &t,
                double &tLastNetworkUpdate, std::mt19937_64 &generator,
                std::vector<double> &T, std::vector<size_t> &C,
                std::vector<std::vector<std::pair<double, size_t>>> &S)

{
    std::vector<double> propensities(M, 0);
    std::vector<std::pair<double, Edge>> propDel;
    std::vector<std::pair<double, Edge>> propAdd;

    for (size_t ind = 0; ind < n; ind++)
    {
        propDel = contNetwork.getEdgeDeletionRateSum();
        propAdd = contNetwork.getEdgeAdditionRateSum();

        propensities.at(0) = propDel.back().first;
        propensities.at(1) = propAdd.back().first;

        double propensitiesSum = propensities.at(0) + propensities.at(1);

        if (propensitiesSum == 0)
        {
            t = tEnd;
            tLastNetworkUpdate = tEnd; //used to update netw.Upd.Time
            break;
        }

        double r = sampleRandUni(generator);

        double proposedTime = std::log(1 / r) * 1 / propensitiesSum ;

        if (t + proposedTime > tEnd)
        {
            t = tEnd;
            tLastNetworkUpdate = t;
            break;
        }

        t += proposedTime;
        tLastNetworkUpdate = t; //used to update netw. Upd.Time

        double rbound =  propensitiesSum * sampleRandUni(generator);
        //deletion
        if (propensities.at(0) >= rbound)
        {
            auto edgeIterator = std::lower_bound(propDel.begin(), propDel.end(), rbound, lambdaLess);

            if (edgeIterator == propDel.end() || edgeIterator == propDel.begin())
            {
                std::string msg = "ERROR: INVALID SSA del!";
                throw std::domain_error(msg);
            }
            contNetwork.removeEdge(edgeIterator->second);
            C.at(0)++;
        }
        else
        {
            auto edgeIterator = std::lower_bound(propAdd.begin(), propAdd.end(), rbound - propensities.at(0), lambdaLess);

            if (edgeIterator == propAdd.end() || edgeIterator == propAdd.begin())
            {
                std::string msg = "ERROR: INVALID SSA add!";
                throw std::domain_error(msg);
            }
            contNetwork.addEdge(edgeIterator->second);
            C.at(1)++;
        }


        for (size_t i = 0; i < M; i ++)
        {
            T.at(i) = T.at(i) + propensities.at(i) * proposedTime;

            auto new_end = std::remove_if(S.at(i).begin(), S.at(i).end(),
                                          [x = std::as_const(T.at(i))](const std::pair<double, size_t> & s)
                                          { return s.first <= x; });
            S.at(i).erase(new_end, S.at(i).end());
            S.at(i).emplace(S.at(i).begin(), std::make_pair(T.at(i), C.at(i)));
        }
    }
}
double Anderson::getTau(const std::vector<double> &props, const std::vector<size_t> &X)
{
    double gi = 1.0;

    double mu_edges = - props.at(0) + props.at(1);
    double mu_compl_edges = props.at(0) - props.at(1);

    double sigma_square = props.at(0) + props.at(1);

    double max_edges = std::max(epsilon * X.at(0) / gi, 1.0);
    double max_compl_edges = std::max(epsilon * X.at(1) / gi, 1.0);

    double tau = std::min({max_edges /abs(mu_edges),
                    max_edges * max_edges / sigma_square,
                    max_compl_edges /abs(mu_compl_edges),
                    max_compl_edges * max_compl_edges / sigma_square});
    return tau;
}


void Anderson::updateNetwork(std::vector<size_t> k, std::mt19937_64 &generator,
                   std::vector<std::pair<double, Edge>> &propAdd,
                   std::vector<std::pair<double, Edge>> &propDel,
                   ContactNetwork & contNetwork)
{
    std::vector<size_t> order;
    for (size_t ind = 0; ind < k.size(); ind++)
    {
        order.insert(order.end(), k.at(ind), ind);
    }
    std::shuffle(order.begin(), order.end(), generator);

    for (auto i : order)
    {
        if (i == 0)
        {
            double rbound = propDel.back().first * sampleRandUni(generator);
            auto edgeIterator = std::lower_bound(propDel.begin(), propDel.end(), rbound, lambdaLess);

            if (edgeIterator == propDel.end() || edgeIterator == propDel.begin())
            {
                std::string msg = "ERROR: INVALID update del!";
                throw std::domain_error(msg);
            }

            std::pair<int, int> b = contNetwork.removeEdge(edgeIterator->second);
            propDel.erase(edgeIterator);

            Edge e = contNetwork.getComplementEdge(b.first, b.second);
            propAdd.emplace_back(propAdd.back().first + contNetwork.getEdgeAdditionRate(e), e);

        }
        else if (i == 1)
        {
            double rbound = propAdd.back().first * sampleRandUni(generator);
            auto edgeIterator = std::lower_bound(propAdd.begin(), propAdd.end(), rbound, lambdaLess);

            if (edgeIterator == propAdd.end() || edgeIterator == propAdd.begin())
            {
                std::string msg = "ERROR: INVALID update add 2!";
                throw std::domain_error(msg);
            }

            std::pair<int, int> b = contNetwork.addEdge(edgeIterator->second);
            propAdd.erase(edgeIterator);

            Edge e = contNetwork.getEdge(b.first, b.second);
            propDel.emplace_back(propDel.back().first + contNetwork.getEdgeDeletionRate(e), e);
        }
    }
}

double Anderson::updateTau(double tau, size_t numOfEdgesExist, size_t numOfEdgesComplement, const std::vector<size_t>& NN)
{
    bool pass2 = NN.at(0) <= std::max(0.75 * epsilon * numOfEdgesExist, 1.0) &&
                 NN.at(1) <= std::max(0.75 * epsilon * numOfEdgesComplement, 1.0);
    double result = tau;
    if (pass2)
    {
        if (tau <= 1)
        {
            result = std::pow(result, q1);
        }
        else
        {
            result = std::pow(result, q2);
        }
    }
    else
    {
        result = result * p1;
    }
    return result;
}

std::vector<size_t> Anderson::getChange(size_t M, const std::vector<std::vector<std::pair<double, size_t>>> &S,
                         const std::vector<double> &T, const std::vector<size_t> &C,
                         const std::vector<double> &propensities, std::vector<size_t> &row,
                         double tau, std::mt19937_64 &generator)
{
    std::vector<size_t> change(M, 0);

    for (size_t i = 0; i < M; i ++)
    {
        std::vector<std::pair<double, size_t>> Sk = S.at(i);
        size_t B = Sk.size() - 1;
        if (propensities.at(i) * tau + T.at(i) >= Sk.at(B).first)
        {
            std::poisson_distribution<size_t> poiss(propensities.at(i) * tau + T.at(i) - Sk.at(B).first);
            change.at(i) = poiss(generator) + Sk.at(B).second - C.at(i);
            row.at(i) = B;
        }
        else
        {
            auto iterator = std::upper_bound(Sk.begin(), Sk.end(), propensities.at(i) * tau + T.at(i), [](double value, const std::pair<double, size_t> &a) { return a.first > value; });

            if (iterator == Sk.end() || iterator == Sk.begin())
            {
                std::string msg = "ERROR: Invalid";
                throw std::domain_error(msg);
            }

            size_t index = std::distance(Sk.begin(), iterator);
            double r = (T.at(i) + propensities.at(i) * tau - Sk.at(index - 1).first) / (Sk.at(index).first - Sk.at(index - 1).first);

            std::binomial_distribution<size_t> binom(Sk.at(index).second - Sk.at(index - 1).second, r);

            change.at(i) = binom(generator) + Sk.at(index - 1).second - C.at(i);
            row.at(i) = index - 1;
        }
    }

    return change;
}

void Anderson::acceptLeap(double &t, double &tLastNetworkUpdate, double &tau, size_t M,
                ContactNetwork & contNetwork,
                std::vector<std::vector<std::pair<double, size_t>>> &S,
                std::vector<double> &T, std::vector<size_t> &C,
                std::vector<size_t> &row, std::vector<double> &propensities,
                const std::vector<size_t> &change,
                std::vector<std::pair<double, Edge>> &propAdd,
                std::vector<std::pair<double, Edge>> &propDel,
                std::mt19937_64 &generator)
{
    for (size_t i = 0; i < M; i ++)
    {
        T.at(i) += propensities.at(i) * tau;
        C.at(i) += change.at(i);

        S.at(i).erase(S.at(i).begin(),  S.at(i).begin() + row.at(i) + 1);
        S.at(i).emplace(S.at(i).begin(), std::make_pair(T.at(i), C.at(i)));

    }

    t+= tau;
    tau = updateTau(tau, propDel.size() - 1, propAdd.size() - 1, change);


    updateNetwork(change, generator, propAdd, propDel,contNetwork);
    tLastNetworkUpdate = t;
    propDel = contNetwork.getEdgeDeletionRateSum();
    propAdd = contNetwork.getEdgeAdditionRateSum();

    propensities.at(0) = propDel.back().first;
    propensities.at(1) = propAdd.back().first;
}

void Anderson::rejectLeap(size_t M,  double &tau, std::vector<std::vector<std::pair<double, size_t>>> &S,
                const std::vector<double> &T, const std::vector<size_t> &C,
                const std::vector<double> &propensities, std::vector<size_t> &row,
                const std::vector<size_t> &change)
{
    for (size_t i = 0; i < M; i ++)
    {
        std::pair<double, size_t> toInsert(propensities.at(i) * tau + T.at(i), C.at(i) + change.at(i));
        if (row.at(i) == S.at(i).size() - 1)
        {
            S.at(i).emplace_back(toInsert);
        }
        else
        {
            S.at(i).emplace(S.at(i).begin() + row.at(i) + 1, toInsert);
        }
    }
    tau *=  p;
}