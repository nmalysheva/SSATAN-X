/**
 * Created by Malysheva, Nadezhda on 2019-07-28.
 *
 * Class ContactNetwork describes a time evolving network.
 * Nodes - Species are connected to each other with edges. Graph structure change through time by
 * assembling /disassembling edges.
 * Graph is simple, i.e only one edg between each pair of nodes. Graph is undirected.
 * infection is not possible (for example, S-S, I-I edges) are undirected and ones with transmission rate > 0
 * for instance edge between Infected (I) and Susceprible (S) is directed (as transmission only possible one direction)
 * Uses LEMON library for graph representation.
 * Graph is represented by two graphs - actual network & complement network. Complement NW is used for
 * more convenient and direct addition of the edges.
*/

#ifndef ALGO_CONTACTNETWORK_H
#define ALGO_CONTACTNETWORK_H

#include <lemon/list_graph.h>
#include <lemon/maps.h>
#include <random>
#include "Specie.h"
#include "utilities/types.h"
#include "utilities/Settings.h"
#include <lemon/adaptors.h>



class ContactNetwork {

public:

    ContactNetwork(const Settings& settings) :transmissionRates(graph),
                                   population(graph),
                                   filterExistingEdges(graph, false),
                                   filterComplementEdges(graph, true),
                                   existingEdges(graph, filterExistingEdges),
                                   complementEdges(graph, filterComplementEdges)


                                   {
                                       init(settings);
                                   };


    size_t  size() const; //@return amount of nodes
    size_t  countByState(Specie::State st) const;  //@return amount of I/S/R etc. species in network
    size_t  countEdges() const;//@return amount of edges

/* calculates cumulative sum of rates of particular reactions.
 * Used in SSA & SSATANX to find reaction being executed
 * @return a vector of pairs <cumulative sum, instance> for instances of question.
 * first element of the vector is always pair <0, INVALID> for convenience
 */
    std::vector<std::pair<double, Edge>> getTransmissionRateSum() const; //transmission
    std::vector<std::pair<double, Edge>> getEdgeDeletionRateSum()const;
    std::vector<std::pair<double, Edge>> getEdgeAdditionRateSum()const;
    std::vector<std::pair<double, Node>> getDeathRateSum()const;
    std::vector<std::pair<double, Node>> getDiagnosisRateSum()const;

    double  getBirthRateSum()const;

/*
 * @return highest possible transmission rate in network. Since now rate is static,
 * just return transmission rate given by initiating Contact Network.
 */
    double  getTransmissionRateLimit() const;

/*
 * @return  max. number of contacts susceptible nodes can have during time t.
 */

    double  getMaxContactsLimitByState(Specie::State state, double t) const;


 /*
 * Adding edge to the network. input - reference to the edge from complement network
 * @return pair of ids of nodes that was connected by given edge
 */
    std::pair<int, int> addEdge(Edge & complementEdge);

    /*
    * Removing edge from the network. input - reference to the edge from actual network
    * @return pair of ids of nodes that was disconnected
    */

    std::pair<int, int> removeEdge(Edge & edge);

    /*
    * Removing node from the network. input - reference to the node from actual network
    */
    void removeNode(Node & node);


    /*
    * Executing particular reactions.
    */
    void executeTransmission(Edge & edge, double time);
    void executeDiagnosis(Node & node, double time);
    void executeDeath(Node & node);
    void executeBirth(double rStart, double rBound);

    /*
     * Gets degree distribution of the network.
     * @return vector of the degrees of each node in the network
    */
    std::vector<specieState> getNetworkState() const;

    double  getEdgeAdditionRate(const Edge &complementEdge) const;
    double  getEdgeDeletionRate(const Edge &networkEdge) const;

    Edge getComplementEdge(int a, int b); //@return complement edge by given nodes ids
    Edge getEdge(int a, int b);//@return edge of actual network by given nodes ids


private:

    void initRates(double transmRate, double diagnRate, double dRate, double bRate);

    double getMeanEdgeAdditionRate (const Node &complementNode) const;
    double getMeanEdgeDeletionRate (const Node &networkNode) const;

    double numberOfContactEstimation(double meanEdgeAdditionRate, double meanEdgeDeletionRate, double t,
                                     double Cmax, double C0) const;

    /**
 * calculates analytical expectation  for number of contacts specie can have after time t,
 * given current specie parameters - rates of establishing and loosing contacts & number of current contacts.
 * Considering specie has:
 * a - rate of establish a contact
 * b - rate of loose a contact
 * Cmax - maximum number of contacts
 * c - current number of contacts.
 *
 * then analytical ODE for expectation  of contacts of separated specie w/o external factors:
 * @see PAPER
 * @return variance
 */
    /**
    * calculates analytical expectation  for number of contacts specie can have after time t,
    * given current specie parameters - rates of establishing and loosing contacts & number of current contacts.
    * Considering specie has:
    * a - rate of establish a contact
    * b - rate of loose a contact
    * Cmax - maximum number of contacts
    * c - current number of contacts.
    *
    * then analytical ODE for expectation  of contacts of separated specie w/o external factors:
    * @see PAPER
    * @return expectation
    *
    */
    /**
 * calculates analytical solution for max. number of contacts specie can have after time t,
 * given current specie parameters - rates of establishing and loosing contacts & number of current contacts.
 * Considering specie has:
 * a - rate of establish a contact
 * b - rate of loose a contact
 * Cmax - maximum number of contacts
 * c - current number of contacts.
 *
 * then analytical ODE for number of contacts of separated specie w/o external factors:
 * @see PAPER
 *
 */
    double getExtremaPoint(double meanEdgeAdditionRate, double meanEdgeDeletionRate, double t,
                           double Cmax, double C0) const;


    void init(const Settings&settings);


    lemon::ListGraph graph;

    lemon::ListGraph::EdgeMap<double> transmissionRates;

    std::uniform_real_distribution<double> looseContactDistribution;
    std::uniform_real_distribution<double> createContactDistribution;

    double transmissionRate;
    double diagnosisRate;
    double birthRate;
    std::vector<double> deathRate;

    lemon::ListGraph::NodeMap<Specie> population;
    lemon::ListGraph::EdgeMap<bool> filterExistingEdges;
    lemon::ListGraph::EdgeMap<bool> filterComplementEdges;
    lemon::FilterEdges<lemon::ListGraph> existingEdges;
    lemon::FilterEdges<lemon::ListGraph> complementEdges;

};


#endif //ALGO_CONTACTNETWORK_H
