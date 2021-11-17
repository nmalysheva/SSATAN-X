//
// Created by Malysheva, Nadezhda on 10.07.20.
//

#ifndef ALGO_UTILITY_H
#define ALGO_UTILITY_H

#include <random>

auto lambdaLess = []<typename T>(const std::pair<double, T> &a,  double value) { return a.first < value; };

double sampleRandUni(std::mt19937_64 &generator);
#endif //ALGO_UTILITY_H
