//
// Created by Malysheva, Nadezhda on 14.10.21.
//
#include <sstream>
#include <fstream>
#include <string>
#include "nlohmann/json.h"
#include "Settings.h"


uint Settings::getSeed() const
{
    return seed;
}

distParameters Settings::getLooseConactRateParameters()const
{
    return looseContactParameters;
}

distParameters Settings::getNewConactRateParameters()const
{
    return newContactParameters;
}

double Settings::getSimulationTime() const
{
    return simulationTime;
}

size_t Settings::getNumberOfEdges() const
{
    return numOfEdges;
}

/*double Settings::getBirthRate() const
{
    return birthRate;
}*/

double Settings::getDiagnosisRate() const
{
    return diagnosisRate;
}

double Settings::getTransmissionRate() const
{
    return transmissionRate;
}

std::unordered_map<std::string, SpecieSettings> Settings::getStateSettings() const
{
    return stateSettings;
}


void Settings::parseSettings(const std::string & configFileName)
{
    std::ifstream confFile;
    confFile.open(configFileName);
    std::stringstream jsonStr;
    jsonStr << confFile.rdbuf();

    auto jsonObj = nlohmann::json::parse(jsonStr.str());
    auto species = jsonObj.at("species");
    for (auto specie: species)
    {
        SpecieSettings st;
        st.amount = specie.at("amount").get<size_t>();

        st.deathRate = specie.at("death_rate").get<double>();
        stateSettings.emplace(specie.at("state"), st);
    }

    seed = jsonObj.at("seed").get<uint>();
    //birthRate = jsonObj.at("birth_rate").get<double>();
    numOfEdges = jsonObj.at("initial_edges").get<size_t>();
    simulationTime = jsonObj.at("simulation_time").get<double>();

    looseContactParameters.a = jsonObj.at("loose_contact_rate")[0].get<double>();
    looseContactParameters.b = jsonObj.at("loose_contact_rate")[1].get<double>();

    newContactParameters.a = jsonObj.at("new_contact_rate")[0].get<double>();
    newContactParameters.b = jsonObj.at("new_contact_rate")[1].get<double>();

    transmissionRate = jsonObj.at("transmission_rate").get<double>();
    diagnosisRate = jsonObj.at("diagnosis_rate").get<double>();
}

rateParameters Settings::parseDistribution(const nlohmann::json& distInfo)
{
    rateParameters rp;
    if (! distInfo.is_number() && distInfo.is_string())
    {
        std::string dist = distInfo.get<std::string>();
        if (dist.find("U", 0) == 0)
        {
            rp.distribution = "uniform";
        }

        size_t first = dist.find("(");
        size_t last = dist.find(")");
        std::string strNew = dist.substr (first + 1,last-first - 1);

        std::string aa = strNew.substr(0, strNew.find(","));
        double aaa = stod(strNew.substr(0, strNew.find(",") ));
        rp.parameters.push_back(aaa);
        double bbb = stod(strNew.substr(strNew.find(",") + 1, strNew.length() -  strNew.find(",") - 1));
        rp.parameters.push_back(bbb);

    }
    else if (distInfo.is_number())
    {
        rp.distribution = "none";
        double aa = distInfo.get<double>();
        rp.parameters.push_back(aa);

    }
    else
    {
        std::string msg = "Invalid parameters. Provide constant or distribution";
        throw std::domain_error(msg);
    }

    return rp;
}

