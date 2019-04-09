//
// Created by jkoenig on 09/04/19.
//
/*
#ifndef MAGNETICSPHERES_METROPOLIS_H
#define MAGNETICSPHERES_METROPOLIS_H

#include <memory>
#include "cluster.h"

class Metropolis {
private:
    std::shared_ptr<cluster> cluster;
    double beta, energy_current;
    int step_count, attempt_count, dipole_count;

    //attempt a minimization step with a completely random dipole
    void minimization_attempt_random();
    void minimization_attempt_random_angle();


public:
    Metropolis(); // generates random cluster with 8 spheres

    double get_energy(){return energy_current;}
    cluster get_cluster_snapshot(){return *cluster;}

    //minimize randomly with n steps
    void minimize_random(int n);

    void test();
};


#endif //MAGNETICSPHERES_METROPOLIS_H
*/