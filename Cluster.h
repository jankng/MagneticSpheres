//
// Created by jan on 4/5/19.
//

#ifndef MAGNETICSPHERES_CLUSTER_H
#define MAGNETICSPHERES_CLUSTER_H


#include <vector>
#include "Dipole.h"

class Cluster {
private:
    std::vector<Dipole> config;
    double n, diam;
    double energy;
public:
    Cluster(int n);
    Dipole* GetDipole(int id);

    void Print();

    double ComputeEnergy();
    double dummy();
};


#endif //MAGNETICSPHERES_CLUSTER_H
