//
// Created by jan on 4/5/19.
//

#ifndef MAGNETICSPHERES_CLUSTER_H
#define MAGNETICSPHERES_CLUSTER_H


#include <vector>
#include <memory>
#include "Dipole.h"

class Cluster {
private:
    std::shared_ptr< std::vector<Dipole> > config;
    int n;
    double diam;
    double energy;
public:
    Cluster(int n);
    Dipole* GetDipole(int id);

    void Print();

    double ComputeEnergy();
    double dummy();
    
    void MetropolisStep();
};


#endif //MAGNETICSPHERES_CLUSTER_H
