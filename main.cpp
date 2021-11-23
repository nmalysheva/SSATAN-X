#include <iostream>
#include <chrono>
#include <fstream>
#include <string>

#include "contact_network/ContactNetwork.h"
#include "algorithms/SSA.h"
#include "algorithms/SSATANX.h"
#include "utilities/types.h"
#include "utilities/Settings.h"

void saveInitialStates(nlohmann::ordered_json &output, const ContactNetwork &contNetwork, const Settings& settings)
{
    output["initial_states"][Specie::S]["S"] = contNetwork.countByState(Specie::S);
    output["initial_states"][Specie::I]["I"] = contNetwork.countByState(Specie::I);
    output["initial_states"][Specie::D]["D"] = contNetwork.countByState(Specie::D);
    output["start_edges"] = contNetwork.countEdges();

    output["rate_of_make_a_new_contact"] = {settings.getNewConactRateParameters().a, settings.getNewConactRateParameters().b};
    output["rate_of_loose_a_contact"] = {settings.getLooseConactRateParameters().a, settings.getLooseConactRateParameters().b};
    output["birth_rate"] = settings.getBirthRate();
    output["diagnosis_rate"] = settings.getDiagnosisRate();
    output["transmission_rate"] = settings.getTransmissionRate();
}

void saveOutput(nlohmann::ordered_json &output, const ContactNetwork &contNetwork, const NetworkStorage &nwStorage)
{
    output["final_states"][Specie::S]["S"] = contNetwork.countByState(Specie::S);
    output["final_states"][Specie::I]["I"] = contNetwork.countByState(Specie::I);
    output["final_states"][Specie::D]["D"] = contNetwork.countByState(Specie::D);


    size_t i = 0;
    for (const auto &item : nwStorage)
    {
        size_t j = 0;
        output["networkStates"][i]["time"] = item.first;
        for (const auto &spcs : item.second)
        {
            output["networkStates"][i]["nw_states"][j]["id"] = spcs.id;
            output["networkStates"][i]["nw_states"][j]["state"] = spcs.sp.getState();
            output["networkStates"][i]["nw_states"][j]["rate_of_make_a_new_contact"] = spcs.sp.getNewContactRate();
            output["networkStates"][i]["nw_states"][j]["rate_of_loose_a_contact"] = spcs.sp.getLooseContactRate();
            output["networkStates"][i]["nw_states"][j]["death_rate"] = spcs.sp.getDeathRate();
            output["networkStates"][i]["nw_states"][j]["diagnosis_rate"] = spcs.sp.getDiagnosisRate();
            output["networkStates"][i]["nw_states"][j]["neighbors"] = spcs.contacts;
            j ++;
        }
        i++;

    }
}

void executeSSA(const Settings& settings)
{
    ContactNetwork contNetwork(settings);

    nlohmann::ordered_json output;
    saveInitialStates(output, contNetwork, settings);

    NetworkStorage nwStorage;
    nwStorage.reserve(1e6 + 1);

    auto start_time = std::chrono::high_resolution_clock::now();
    auto const filename = std::chrono::system_clock::now().time_since_epoch().count();
    SSA().execute(0, settings.getSimulationTime(), contNetwork, nwStorage);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;

    output["duration_in_milliseconds"] = std::chrono::duration <double, std::milli> (time).count();
    saveOutput(output, contNetwork, nwStorage);

    std::string fileName = "SSA_" + std::to_string(filename) + ".txt";

    std::ofstream newFile;
    newFile.open(fileName);
    newFile <<  output << std::endl;
    newFile.close();

}

void executeSSATANX(const Settings& settings)
{
    ContactNetwork contNetwork(settings);

    nlohmann::ordered_json output;
    saveInitialStates(output, contNetwork, settings);

    NetworkStorage nwStorage;
    nwStorage.reserve(1e6 + 1);

    size_t nRejections = 0;
    size_t nAcceptance = 0;
    size_t nThin = 0;


    auto start_time = std::chrono::high_resolution_clock::now();
    auto const filename = std::chrono::system_clock::now().time_since_epoch().count();
    SSATANX().execute(0, settings.getSimulationTime(), contNetwork, nwStorage,nRejections, nAcceptance, nThin);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;

    std::string fileName = "NSA_" + std::to_string(filename) + ".txt";

    output["duration_in_milliseconds"] = std::chrono::duration <double, std::milli> (time).count();
    output["accepted"] = nAcceptance;
    output["rejected"] = nRejections;
    output["thined"] = nThin;

    saveOutput(output, contNetwork, nwStorage);

    std::ofstream newFile;
    newFile.open(fileName);
    newFile <<  output << std::endl;
    newFile.close();

}

void viralDynamics(int argc, char* argv[])
{
    std::string mode = std::string(argv[2]);
    //size_t simulationNumber = std::stoi(argv[3]);
    std::string fileName = std::string(argv[1]);
    Settings settings;
    settings.parseSettings(fileName);
    if (mode=="-SSA")
    {
        executeSSA(settings);
    }
    else if (mode=="-SSX")
    {
        executeSSATANX(settings);
    }
    else
    {
        std::string msg = "Invalid algorthm specified";
        throw std::domain_error(msg);
    }
}

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        std::string msg = "Invalid parameters";
        throw std::domain_error(msg);
    }
    else
    {
        viralDynamics(argc, argv);
    }
    return 0;
}
