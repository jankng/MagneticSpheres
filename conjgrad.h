//
// Created by jkoenig on 29/04/19.
//

#ifndef MAGNETICSPHERES_CONJGRAD_H
#define MAGNETICSPHERES_CONJGRAD_H

#include "definitions.h"
#include "cluster.h"
#include "dipole.h"

class conjgrad {
private:
    cluster* cl;
    std::vector<double> grad;

    // gives energy in direction of gradient
    double energy_in_grad_dir(double t);

public:
    explicit conjgrad(cluster* given);
    void dosomehting();
};


#endif //MAGNETICSPHERES_CONJGRAD_H
