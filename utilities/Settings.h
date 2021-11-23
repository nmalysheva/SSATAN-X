//
// Created by Malysheva, Nadezhda on 14.10.21.
//

#ifndef ALGO_SETTINGS_H
#define ALGO_SETTINGS_H

#include <vector>
#include <map>
#include "nlohmann/json.h"

struct rateParameters
{
    std::string distribution;
    std::vector<double> parameters;
    //double parameters;

    void from_json(const nlohmann::json& j, rateParameters& p)
    {
        j.at("distribution").get_to(p.distribution);
        j.at("parameters").get_to(p.parameters);
    }
};

struct distParameters
{
    double a;
    double b;
};

struct SpecieSettings{
    size_t amount;
    double deathRate;
};

struct Transition {
    std::string fromState;
    std::string toState;
    rateParameters rate;

};

struct Interaction {
    std::vector<std::string> fromStates;
    std::vector<std::string> toStates;
    rateParameters rate;

};

class Settings {
public:
    double getSimulationTime() const;
    size_t getNumberOfEdges() const;
    std::unordered_map<std::string, SpecieSettings> getStateSettings() const;
    double getDiagnosisRate() const;
    double getTransmissionRate() const;
    //double getBirthRate() const;
    uint getSeed() const;

    distParameters getLooseConactRateParameters()const;
    distParameters getNewConactRateParameters()const;

    void parseSettings(const std::string & configFileName);

private:
    double simulationTime;
    size_t numOfEdges;
    std::unordered_map<std::string, SpecieSettings> stateSettings;
    double diagnosisRate;
    double transmissionRate;
    //double birthRate;
    uint seed;

    distParameters looseContactParameters;
    distParameters newContactParameters;

    rateParameters parseDistribution(const nlohmann::json& distInfo);
};

#endif //ALGO_SETTINGS_H
