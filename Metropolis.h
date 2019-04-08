//
// Created by jkoenig on 08/04/19.
//

#ifndef MAGNETICSPHERES_METROPOLIS_H
#define MAGNETICSPHERES_METROPOLIS_H

#include <memory>
#include "Cluster.h"

class Metropolis {
private:
    std::shared_ptr<Cluster> cluster;

public:
    void MinimizationStep();
};


#endif //MAGNETICSPHERES_METROPOLIS_H
