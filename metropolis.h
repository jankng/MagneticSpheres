//
// Created by jkoenig on 09/04/19.
//

#ifndef MAGNETICSPHERES_METROPOLIS_H
#define MAGNETICSPHERES_METROPOLIS_H

#include <memory>
#include "cluster.h"

#define DEFAULT_BETA_START 1

class metropolis {
private:
    std::shared_ptr<cluster> cl;
    int cluster_size;
    double beta, energy_current;
    int step_count, attempt_count;


public:
    explicit metropolis(int n); // generates random cluster with n spheres

    double get_energy(){return energy_current;}
};


#endif //MAGNETICSPHERES_METROPOLIS_H
