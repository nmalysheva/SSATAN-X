/*
 * Created by Malysheva, Nadezhda on 21.07.20.
*/

#ifndef ALGO_ANDERSONTAULEAP_H
#define ALGO_ANDERSONTAULEAP_H


#include "contact_network/ContactNetwork.h"
#include "utilities/types.h"

class Anderson{
public:
static void AndersonTauLeap(double &tLastNetworkUpdate, double tEnd, ContactNetwork & contNetwork, std::mt19937_64 &generator);

private:

Anderson(){};

static void updateNetwork(std::vector<size_t> k, std::mt19937_64 &generator,
                   std::vector<std::pair<double, Edge>> &propAdd,
                   std::vector<std::pair<double, Edge>> &propDel,
                   ContactNetwork & contNetwork);


static void executeSSA(size_t n, size_t M, double tEnd, ContactNetwork & contNetwork, double &t,
                double &tLastNetworkUpdate, std::mt19937_64 &generator,
                std::vector<double> &T, std::vector<size_t> &C,
                std::vector<std::vector<std::pair<double, size_t>>> &S);

static double getTau(const std::vector<double> &props, const std::vector<size_t> &X);
static double updateTau(double tau, size_t numOfEdgesExist, size_t numOfEdgesComplement, const std::vector<size_t>& NN);

static std::vector<size_t> getChange(size_t M, const std::vector<std::vector<std::pair<double, size_t>>> &S,
                                            const std::vector<double> &T, const std::vector<size_t> &C,
                                            const std::vector<double> &propensities, std::vector<size_t> &row,
                                            double tau, std::mt19937_64 &generator);

static void acceptLeap(double &t, double &tLastNetworkUpdate, double &tau, size_t M,
                    ContactNetwork & contNetwork,
                    std::vector<std::vector<std::pair<double, size_t>>> &S,
                    std::vector<double> &T, std::vector<size_t> &C,
                    std::vector<size_t> &row, std::vector<double> &propensities,
                    const std::vector<size_t> &change,
                    std::vector<std::pair<double, Edge>> &propAdd,
                    std::vector<std::pair<double, Edge>> &propDel,
                    std::mt19937_64 &generator);

static void rejectLeap(size_t M,  double &tau, std::vector<std::vector<std::pair<double, size_t>>> &S,
                    const std::vector<double> &T, const std::vector<size_t> &C,
                    const std::vector<double> &propensities, std::vector<size_t> &row,
                    const std::vector<size_t> &change);

private:

    static constexpr double p = 0.75;
    static constexpr double p1 = 0.9;
    static constexpr double q1 = 0.98;
    static constexpr double q2 = 1.02;

    static constexpr double epsilon = 0.03;

};
#endif //ALGO_ANDERSONTAULEAP_H