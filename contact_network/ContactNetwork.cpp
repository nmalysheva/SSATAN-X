//
// Created by Malysheva, Nadezhda on 2019-07-28.
//

#include <random>
#include <unistd.h>
#include <algorithm>
#include <stdexcept>
#include <lemon/full_graph.h>
#include <lemon/adaptors.h>

#include "ContactNetwork.h"

void ContactNetwork::init(const Settings&settings)
{
    std::unordered_map<std::string, SpecieSettings> statesSettings = settings.getStateSettings();

    size_t  nPopulation = statesSettings.at("S").amount +
            statesSettings.at("D").amount + statesSettings.at("I").amount;

    // init network graph (full graph) and two filters - existing edges and edges to add
    lemon::FullGraph fullG(nPopulation);
    lemon::GraphCopy<lemon::FullGraph, lemon::ListGraph> cg(fullG, graph);
    lemon::FullGraph::EdgeMap<bool> filterFalseEdges(fullG, false);
    lemon::FullGraph::EdgeMap<bool> filterTrueEdges(fullG, true);
    cg.edgeMap(filterFalseEdges, filterExistingEdges);
    cg.edgeMap(filterTrueEdges, filterComplementEdges);
    cg.run();

    std::mt19937_64 generator;
    if (settings.getSeed() == 0)
    {
        generator = std::mt19937_64(::time(nullptr) * getpid());
    }
    else
    {
        generator = std::mt19937_64(settings.getSeed());
    }

    transmissionRate = settings.getTransmissionRate();
    birthRate = settings.getBirthRate();
    diagnosisRate = settings.getDiagnosisRate();
    deathRate.push_back(statesSettings.at("S").deathRate);
    deathRate.push_back(statesSettings.at("I").deathRate);
    deathRate.push_back(statesSettings.at("D").deathRate);

    /* rates of assemble and disassemble edges are sampled from distributions
     *
     */
    std::uniform_real_distribution<double>::param_type looseRateParameters(settings.getLooseConactRateParameters().a, settings.getLooseConactRateParameters().b);
    looseContactDistribution.param(looseRateParameters);

    std::uniform_real_distribution<double>::param_type newRateParameters(settings.getNewConactRateParameters().a, settings.getNewConactRateParameters().b);
    createContactDistribution.param(newRateParameters);

    size_t maxContacts = nPopulation - 1; //max. num of contacts for each node

    lemon::ListGraph::NodeIt nIt(graph);
    for (size_t i = 0; i <  statesSettings.at("S").amount; i ++)
    {
        double looseContRate = looseContactDistribution(generator);
        double newContRate = createContactDistribution(generator);

        Specie::State st = Specie::S;
        Specie sp = Specie(maxContacts, 0, deathRate.at(st), newContRate, looseContRate, st, 0);
        population[nIt] = sp;
        ++nIt;
    }

    for (size_t i = 0; i <  statesSettings.at("I").amount; i ++)
    {
        double looseContRate = looseContactDistribution(generator);
        double newContRate = createContactDistribution(generator);

        Specie::State st = Specie::I;
        Specie sp = Specie(maxContacts, 0, deathRate.at(st), newContRate, looseContRate, st, diagnosisRate);
        population[nIt] = sp;
        ++nIt;
    }

    for (size_t i = 0; i <  statesSettings.at("D").amount; ++i)
    {
        double looseContRate = looseContactDistribution(generator);
        double newContRate = createContactDistribution(generator);

        Specie::State st = Specie::D;

        Specie sp = Specie(maxContacts, 0, deathRate.at(st), newContRate * 0.3, looseContRate, st, 0);
        population[nIt] = sp;
        ++nIt;
    }


    std::vector<int> v(lemon::countEdges(complementEdges));
    std::iota (std::begin(v), std::end(v), 0);
    std::shuffle(v.begin(), v.end(), generator);

    for (size_t i = 0; i < settings.getNumberOfEdges(); i ++)
    {
        int edgeId = v.at(i);
        Edge cEdge = complementEdges.edgeFromId(edgeId);
        if (getEdgeAdditionRate(cEdge) > 0)
        {
            addEdge(cEdge);
        }
    }
}

size_t  ContactNetwork::countByState(Specie::State st) const
{
    size_t result = 0;

    for (lemon::ListGraph::NodeIt nIt(graph); nIt != lemon::INVALID; ++nIt)
    {
        if (population[nIt].getState() == st)
        {
            result++;
        }
    }
    return result;
}


std::vector<std::pair<double, Edge>> ContactNetwork::getTransmissionRateSum() const
{
    std::vector<std::pair<double, Edge>> propCumSum;
    propCumSum.reserve(1e+6);

    //element <0, INVALID>
    std::pair<double, Edge> invalidElem {0, Edge(lemon::INVALID)};
    propCumSum.push_back(invalidElem);

    for (lemon::FilterEdges<lemon::ListGraph>::EdgeIt eIt(existingEdges); eIt != lemon::INVALID; ++eIt)
    {
        double rate = transmissionRates[eIt];
        if (rate > 0)
        {
            propCumSum.emplace_back(propCumSum.back().first + rate, eIt);
        }
    }
    propCumSum.shrink_to_fit();
    return propCumSum;

}

std::vector<std::pair<double, Node>> ContactNetwork:: getDiagnosisRateSum()const
{
    std::vector<std::pair<double, Node>> propCumSum;
    propCumSum.reserve(1e+6);

    //element <0, INVALID>
    std::pair<double, Node> invalidElem {0, Node (lemon::INVALID)};
    propCumSum.push_back(invalidElem);

    for (lemon::FilterEdges<lemon::ListGraph>::NodeIt nIt(existingEdges); nIt != lemon::INVALID; ++nIt)
    {
        double rate = population[nIt].getDiagnosisRate();
        if (rate > 0)
        {
            propCumSum.emplace_back(propCumSum.back().first + rate, nIt);
        }
    }
    propCumSum.shrink_to_fit();
    return propCumSum;
}

std::vector<std::pair<double, Edge>> ContactNetwork::getEdgeDeletionRateSum()const
{
    std::vector<std::pair<double, Edge>> propCumSum;
    propCumSum.reserve(1e+6);

    //element <0, INVALID>
    std::pair<double, Edge> invalidElem {0, Edge (lemon::INVALID)};
    propCumSum.push_back(invalidElem);

    for (lemon::FilterEdges<lemon::ListGraph>::EdgeIt eIt(existingEdges); eIt != lemon::INVALID; ++eIt)
    {
        double rate = getEdgeDeletionRate(eIt);
        if (rate > 0)
        {
            propCumSum.emplace_back(propCumSum.back().first + rate, eIt);
        }
    }
    propCumSum.shrink_to_fit();
    return propCumSum;

}


std::vector<std::pair<double, Edge>> ContactNetwork::getEdgeAdditionRateSum()const
{
    std::vector<std::pair<double, Edge>> propCumSum;
    //int popSize = size();
    //propCumSum.reserve(popSize * (popSize - 1) / 2);
    propCumSum.reserve(1e+6);

    //element <0, INVALID>
    std::pair<double, Edge> invalidElem {0, Edge (lemon::INVALID)};
    propCumSum.push_back(invalidElem);

    for (lemon::FilterEdges<lemon::ListGraph>::EdgeIt eIt(complementEdges); eIt != lemon::INVALID; ++eIt)
    {
        double rate = getEdgeAdditionRate(eIt);
        if (rate > 0)
        {
            propCumSum.emplace_back(propCumSum.back().first + rate, eIt);
        }
    }
    propCumSum.shrink_to_fit();
    return propCumSum;
}

std::vector<std::pair<double, Node>> ContactNetwork::getDeathRateSum()const
{
    std::vector<std::pair<double, Node>> propCumSum;
    //propCumSum.reserve(size());
    propCumSum.reserve(1e+6);

    //element <0, INVALID>
    std::pair<double, Node> invalidElem {0, Node(lemon::INVALID)};
    propCumSum.push_back(invalidElem);

    for(lemon::FilterEdges<lemon::ListGraph>::NodeIt nIt(existingEdges); nIt!=lemon::INVALID; ++nIt)
    {
        double rate = population[nIt].getDeathRate();

        propCumSum.emplace_back(propCumSum.back().first + rate, nIt);
    }

    propCumSum.shrink_to_fit();
    return propCumSum;

}

double ContactNetwork::getTransmissionRateLimit() const
{
    return transmissionRate;
}


double ContactNetwork::getBirthRateSum()const
{
    return birthRate;
}

size_t ContactNetwork::size() const
{
    return lemon::countNodes(existingEdges);
}

std::pair<int, int> ContactNetwork::addEdge(Edge &complementEdge)
{

    //nodes of the given edge in a complement graph
    Node nodeU = graph.u(complementEdge);
    Node nodeV = graph.v(complementEdge);
    filterExistingEdges[complementEdge] = true;
    filterComplementEdges[complementEdge] = false;

    std::pair<int, int> result = std::make_pair(graph.id(nodeU), graph.id(nodeV));

    //for these nodes increase number of contacts
    population[nodeU].incNumberOfContacts();
    population[nodeV].incNumberOfContacts();

    //calculate transmission rate
    double trRate = 0;
    if ((population[nodeU].getState()  == Specie::S  && population[nodeV].getState() == Specie::I) ||
        (population[nodeU].getState()  == Specie::I  && population[nodeV].getState() == Specie::S))
    {
        trRate = transmissionRate;
    }

    if ((population[nodeU].getState()  == Specie::S  && population[nodeV].getState() == Specie::D) ||
        (population[nodeU].getState()  == Specie::D  && population[nodeV].getState() == Specie::S))
    {
        //for edges Diagnosed "D"- Susceptible "S" transmission rate is 50% lower than for edges "I" - "S"
        trRate = transmissionRate * 0.5;
    }

    transmissionRates[complementEdge] = trRate;
    return result;
}

std::pair<int, int> ContactNetwork::removeEdge(Edge &edge)
{
    if (!graph.valid(edge))
    {
        std::string msg = "ERROR: INVALID EDGE TO DEL!";
        throw std::domain_error(msg);
    }
    //nodes of the given edge in a graph
    Node nodeU = graph.u(edge);
    Node nodeV = graph.v(edge);

    std::pair<int, int> result = std::make_pair(graph.id(nodeU), graph.id(nodeV));

    transmissionRates[edge] = 0;

    filterExistingEdges[edge] = false;
    filterComplementEdges[edge] = true;

    // after removing an edge decrease num. of actual contacts for each of the incident nodes
    population[nodeU].decNumberOfContacts();
    population[nodeV].decNumberOfContacts();

    return result;

}

void ContactNetwork::removeNode(Node &node)
{
    //delete all edges to given node in network
    lemon::ListGraph::IncEdgeIt ieIt(graph, node);

    while (ieIt != lemon::INVALID)
    {
        lemon::ListGraph::IncEdgeIt tmpIt = ieIt;
        Node oppositeNode = graph.oppositeNode(node, ieIt);

        population[oppositeNode].decNumberOfContacts();
        ++ieIt;
        graph.erase(tmpIt);

    }

    graph.erase(node);
}

void ContactNetwork::executeTransmission(Edge & edge, double time)
{
    Node nodeU = graph.u(edge);
    Node nodeV = graph.v(edge);

    Node infectedNode;

    if (population[nodeU].getState()  == Specie::S)
    {
        population[nodeU].changeState(Specie::I, time);
        population[nodeU].setDeathRate(deathRate.at(Specie::I));
        infectedNode = nodeU;
    }

    else if (population[nodeV].getState()  == Specie::S)
    {
        population[nodeV].changeState(Specie::I, time);
        population[nodeV].setDeathRate(deathRate.at(Specie::I));
        infectedNode = nodeV;
    }

    population[infectedNode].setDiagnosisRate(diagnosisRate);
    for(lemon::FilterEdges<lemon::ListGraph>::IncEdgeIt ieIt(existingEdges, infectedNode); ieIt!=lemon::INVALID; ++ieIt)
    {
        transmissionRates[ieIt] = 0;
        Node neighbourNode = existingEdges.oppositeNode(infectedNode, ieIt);
        if(population[neighbourNode].getState() == Specie::S)
        {
            transmissionRates[ieIt] = transmissionRate;//transmitDistribution(generator);

        }
    }

}

void ContactNetwork::executeDiagnosis(Node & node, double time)
{
    population[node].changeState(Specie::D, time);

    population[node].setDiagnosisRate(0); //diagnosed can't be diagnosed anymore

    //adaptivity: as soon as diagnosed, cut all contacts and reduce
    //new contact rate to 30%
    lemon::FilterEdges<lemon::ListGraph>::IncEdgeIt ieIt(existingEdges, node);
    while (ieIt != lemon::INVALID)
    {
        lemon::FilterEdges<lemon::ListGraph>::IncEdgeIt tmpIt(ieIt);
        ++ieIt;
        removeEdge(tmpIt);
    }

    population[node].setDeathRate(deathRate.at(Specie::D));
    population[node].setNewContactRate(
            population[node].getNewContactRate() * 0.3);
}

void ContactNetwork::executeDeath(Node & node)
{
    removeNode(node);

    // after removing node from the population decrease max. number of contacts for each specie.
    for (lemon::ListGraph::NodeIt nIt(graph); nIt != lemon::INVALID; ++nIt)
    {
        population[nIt].setMaxNumberOfContacts(population[nIt].getMaxNumberOfContacts() - 1);
    }

}

void ContactNetwork::executeBirth(double rStart, double rBound)
{
    double result = rStart;
    result += birthRate;
    if (result >= rBound)
    {
    }
}

double  ContactNetwork::getEdgeAdditionRate(const Edge &complementEdge) const
{
    Node nodeU = graph.u(complementEdge);
    Node nodeV = graph.v(complementEdge);

    double sourceRate = population[nodeU].getNewContactRate();
    double targetRate = population[nodeV].getNewContactRate();

    // rate of adding an edge is a multiplication of rates of both incident nodes
    double result = sourceRate * targetRate;
    return result;
}

double  ContactNetwork::getEdgeDeletionRate(const Edge &networkEdge) const
{
    Node nodeU = graph.u(networkEdge);
    Node nodeV = graph.v(networkEdge);

    double sourceRate = population[nodeU].getLooseContactRate();
    double targetRate = population[nodeV].getLooseContactRate();

    // rate of deleting an edge is a multiplication of rates of both incident nodes
    double result = sourceRate * targetRate;
    return result;

}

size_t  ContactNetwork::countEdges() const
{
    return lemon::countEdges(existingEdges);}


std::vector<specieState> ContactNetwork::getNetworkState() const
{
    std::vector<specieState> result;
    //result.reserve(this->size()); //reserving space for vector.
    result.reserve(1e+6);
    for(lemon::FilterEdges<lemon::ListGraph>::NodeIt nIt(existingEdges); nIt!=lemon::INVALID; ++nIt)
    {
        std::vector<int> neighbors;
        neighbors.reserve(1e+6);

        for(lemon::FilterEdges<lemon::ListGraph>::IncEdgeIt e(existingEdges, nIt); e!=lemon::INVALID; ++e)
        {
            Node neighbor = existingEdges.oppositeNode(nIt, e);
            neighbors.push_back(existingEdges.id(neighbor));

        }
        neighbors.shrink_to_fit();

        specieState spState;
        spState.id = existingEdges.id(nIt);
        spState.sp = population[nIt];
        spState.contacts = neighbors;
        result.push_back(spState);
    }
    return result;
}

double   ContactNetwork::getMaxContactsLimitByState(Specie::State state, double t) const
{
    double result = 0;
    double numConMax = static_cast<double> (size() - 1);
    for (lemon::FilterEdges<lemon::ListGraph>::NodeIt nIt(existingEdges); nIt != lemon::INVALID; ++nIt)
    {
        if (population[nIt].getState() == state)
        {
            double meanTheta = getMeanEdgeDeletionRate(nIt);
            double meanLambda = getMeanEdgeAdditionRate(nIt);
            double numConStart = population[nIt].getNumberOfContacts();
            double numConEnd = numberOfContactEstimation(meanLambda, meanTheta, t, numConMax,numConStart);
            double maxCont = std::max(numConStart, numConEnd);

            double numConExtrema = getExtremaPoint(meanLambda, meanTheta, t, numConMax,numConStart);
            maxCont = std::max(numConExtrema, maxCont);

            maxCont = std::min(maxCont, numConMax);

            result+= maxCont;
        }
    }
    return result;
}



Edge ContactNetwork::getComplementEdge(int a, int b)
{
    Node nodeU = complementEdges.nodeFromId(a);
    Node nodeV = complementEdges.nodeFromId(b);

    Edge e = lemon::findEdge(complementEdges, nodeU, nodeV);
    return e;
}
Edge ContactNetwork::getEdge(int a, int b)
{

    Node nodeU = existingEdges.nodeFromId(a);
    Node nodeV = existingEdges.nodeFromId(b);
    Edge e = lemon::findEdge(existingEdges, nodeU, nodeV);
    return e;

}

double ContactNetwork::getMeanEdgeAdditionRate (const Node &complementNode) const
{
    double meanLambda = 0;
    size_t counter = 0;
    for(lemon::FilterEdges<lemon::ListGraph>::IncEdgeIt cieIt(complementEdges, complementNode); cieIt!=lemon::INVALID; ++cieIt)
    //for(lemon::ListGraph::IncEdgeIt cieIt(graph, complementNode); cieIt!=lemon::INVALID; ++cieIt)
    {
        meanLambda += getEdgeAdditionRate(cieIt);
        counter++;
    }

    if (counter > 0)
    {
        meanLambda = meanLambda / counter;
    }

    return meanLambda;
}

double ContactNetwork::getMeanEdgeDeletionRate (const Node &networkNode) const
{
    double meanTheta = 0;

    size_t counter = 0;
    for(lemon::FilterEdges<lemon::ListGraph>::IncEdgeIt ieIt(existingEdges, networkNode); ieIt!=lemon::INVALID; ++ieIt)
    //for(lemon::ListGraph::IncEdgeIt ieIt(graph, networkNode); ieIt!=lemon::INVALID; ++ieIt)
    {
        meanTheta += getEdgeDeletionRate(ieIt);
        counter++;
    }

    if (counter > 0)
    {
        meanTheta = meanTheta / counter;
    }

    return meanTheta;
}

double ContactNetwork::numberOfContactEstimation(double meanEdgeAdditionRate, double meanEdgeDeletionRate, double t, double Cmax, double C0) const
{
    double result = -1.0; //TO DO error indication
    double a = meanEdgeAdditionRate * Cmax;
    double b = meanEdgeAdditionRate + meanEdgeDeletionRate;

    if (t > 0 && b > 0 && a + C0 > 0) {
        double contactsExpectation = a / b - (a / b - C0) * exp(-b * t);
        double contactsVariance = contactsExpectation - C0 * exp(-2 * b * t);
        double numberOfContactEstimation = contactsExpectation + 3 * sqrt(contactsVariance);
        //double numberOfContactEstimation = contactsExpectation + 2 * sqrt(contactsVariance);

        result = numberOfContactEstimation;
    }
    if (result < 0)
    {
        std::string msg = "Something went wrong with expectation";
        std::cout << "meanEdgeAdditionRate=" << meanEdgeAdditionRate << "; meanEdgeDeletionRate=" << meanEdgeDeletionRate << std::endl;
        throw std::domain_error(msg);
    }
    return result;
}

double ContactNetwork::getExtremaPoint(double meanEdgeAdditionRate, double meanEdgeDeletionRate, double t,
                                       double Cmax, double C0) const
{
    double result = -1;

    double a = meanEdgeAdditionRate * Cmax;
    double b = meanEdgeAdditionRate + meanEdgeDeletionRate;

    double rootVal = sqrt(9 * pow(b, 2) * C0 + pow(a - b * C0, 2));
    //double rootVal = sqrt(4 * pow(b, 2) * C0 + pow(a - b * C0, 2));
    double underLog = (b*C0 - a) * rootVal + abs(pow(a, 2)- pow(b*C0, 2));
    underLog = underLog /(2 *b * C0 *rootVal );

    double extrema = 0;
    if (underLog > 0 && a>=0 && b>0 && C0 > a/b)
    {
        extrema = -log(underLog) / b;
    }

    if (extrema > 0 && extrema <=t)
    {
        double numConExtrema = numberOfContactEstimation(meanEdgeAdditionRate, meanEdgeDeletionRate, extrema, Cmax,C0);
        result = std::max(result,numConExtrema);

    }
    return  result;
}