//
// Created by jkoenig on 08/04/19.
//

#include <cmath>
#include "Metropolis.h"
#include "misc.h"
/*
void Metropolis::MinimizationStep() {
    Cluster candidate(cluster.GetN());
    for(int i = 0; i<cluster.GetN(); i++){
        std::vector<double> r = cluster.config[i].GetR();
        candidate.config[i].SetR(r[0], r[1], r[2]);
    }

    if(candidate.ComputeEnergy() < cluster.ComputeEnergy()){
        cluster = candidate;
    }
    else if (misc::SimpleRandom() <= exp(-1*temperature*(cluster.energy - candidate.energy))){
        cluster = candidate;

    }

}
 */
