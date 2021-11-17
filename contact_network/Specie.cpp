//
// Created by Malysheva, Nadezhda on 2019-07-28.
//

#include "Specie.h"

Specie::Specie()
{
    maxNumberOfContacts = 0;
    numberOfContacts = 0;

    deathRate        = 0;
    newContactRate   = 0;
    looseContactRate = 0;

    state = S;
    stateChangeTime =   0;
}


Specie::Specie(/*unsigned char age,*/ size_t maxNumberOfContacts,
                                      size_t numberOfContacts, double deathRate,
               double newContactRate, double looseContactRate, State st, double diagnosisRate)
{
    //setAge(age);
    setMaxNumberOfContacts(maxNumberOfContacts);
    setNumberOfContacts(numberOfContacts);
    setDeathRate(deathRate);
    setNewContactRate(newContactRate);
    setLooseContactRate(looseContactRate);
    setDiagnosisRate(diagnosisRate);

    state = st;
    stateChangeTime = 0;

}

size_t Specie::getMaxNumberOfContacts() const
{
    return maxNumberOfContacts;
}

size_t Specie::getNumberOfContacts() const
{
    return numberOfContacts;
}

double Specie::getNewContactRate() const
{
    double result = newContactRate;

    return result;
}

double Specie::getLooseContactRate() const
{
    double result = looseContactRate;
    return result;
}

double Specie::getDeathRate() const
{
    return deathRate;
}

Specie::State Specie::getState() const
{
    return state;
}

void Specie::setMaxNumberOfContacts(unsigned int maxNumOfCont)
{
    maxNumberOfContacts = maxNumOfCont;
}

void Specie::setNumberOfContacts(size_t nOfCont)
{
    numberOfContacts = nOfCont;
}

void Specie::incNumberOfContacts()
{
    numberOfContacts++;
}
void Specie::decNumberOfContacts()
{
    numberOfContacts--;
}

void Specie::setNewContactRate  (double newContRate)
{
    newContactRate = newContRate;
}

void Specie::setLooseContactRate(double looseContRate)
{
    looseContactRate = looseContRate;
}

void Specie::setDeathRate(double dRate)
{
    deathRate = dRate;
}


void Specie::changeState (State st, double time)
{
    state = st;
    stateChangeTime = time;
}

double Specie::getLastStateChangeTime() const
{
    return stateChangeTime;
}


bool Specie::operator== (const Specie &sp) const
{
    bool result = (
                   maxNumberOfContacts == sp.getMaxNumberOfContacts() &&
                   numberOfContacts == sp.getNumberOfContacts() &&
                   deathRate == sp.getDeathRate() &&
                   newContactRate == sp.getNewContactRate() &&
                   looseContactRate == sp.getLooseContactRate() &&
                   stateChangeTime == sp.getLastStateChangeTime() &&
                   state == sp.getState());
    return result;

}

double Specie::getDiagnosisRate() const
{
    return diagnosisRate;
}

void Specie::setDiagnosisRate(double diagnRate)
{
    diagnosisRate = diagnRate;
}

