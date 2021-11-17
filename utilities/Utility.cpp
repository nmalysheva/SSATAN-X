//
// Created by Malysheva, Nadezhda on 10.07.20.
//
#include "Utility.h"

double sampleRandUni(std::mt19937_64 &generator)
{
    std::uniform_real_distribution<> randuni;
    double r = randuni(generator);
    while (r == 0)
    {
        r = randuni(generator);
    }
    return r;
}

