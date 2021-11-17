/***
 * Created by Malysheva, Nadezhda on 2019-07-28.
 * Class describes a specie (aka node in network).
***/

#ifndef ALGO_SPECIE_H
#define ALGO_SPECIE_H


#include <cstddef>

class Specie {

/**
 * Epidemic state.
 * Depending on the model, can have variety of states,
 * for example in SIR model - states are "S" - susceptible,
 * "I" - infected, "D" - diagnosed.
 */
public:
    enum State {S, I, D}; //Susceptible, Infected, Recovered, Diagnosed

public:

/**
 * Default constructor
 * Create a specie with status "S" by default and all other params set to "0"
 */
    Specie();

    Specie(size_t maxNumberOfContacts,
           size_t numberOfContacts,
           double deathRate,
           double newContactRate,
           double looseContactRate,
           State st = S,
           double diagnosisRate = 0);

    [[nodiscard]] size_t  getMaxNumberOfContacts() const; //@return max. number of contacts specie "allowed to have"
    [[nodiscard]] size_t  getNumberOfContacts() const; //@return current number of contacts specie has

    [[nodiscard]] double getNewContactRate() const;   //@return rate of establishing a new contact
    [[nodiscard]] double getLooseContactRate() const; //@return rate of loosing already existing contacts

    [[nodiscard]] double  getDeathRate() const;  //@return specie's death rate

    [[nodiscard]] double getDiagnosisRate() const;

/**
 * @return time of the last state change.
 * for instance, if infection occurred and state changed from "S" to "I",
 * this function will return the time of infection.
 */
    [[nodiscard]] double  getLastStateChangeTime() const;

/**
 * @return epidemic state of the specie
 */
    [[nodiscard]]
    State  getState() const;

    //setters

    void setMaxNumberOfContacts(unsigned int maxNumOfCont);

    void incNumberOfContacts(); // increase number of current contacts to 1
    void decNumberOfContacts(); // decrease number of current contacts to 1


    void setNewContactRate  (double newContRate);
    void setLooseContactRate(double looseContRate);

    void setDiagnosisRate(double diagnRate);

    void setDeathRate(double dRate);


    //status changes
    void changeState (State st, double time);  //change status of the specie



    //comparison operator for map
    bool operator== (const Specie &sp) const;

    //destructor
    ~Specie() = default;


private:

    void setNumberOfContacts(size_t nOfCont);

private:

    size_t  maxNumberOfContacts; //max. number of contacts specie "allowed to have"
    size_t  numberOfContacts;   //current number of contacts
    double deathRate;          // death rate of the specie
    double newContactRate;  //rate of establishing a new contact, i.e. new edge
    double looseContactRate;//rate to loose a contact, i.e. an edge
    double diagnosisRate;

    double stateChangeTime;  //last time of state change (i.e. S->I, I->R, etc.)

    State state;

};


#endif //ALGO_SPECIE_H
