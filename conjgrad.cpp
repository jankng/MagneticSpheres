//
// Created by jkoenig on 29/04/19.
//

#include "conjgrad.h"
#include <iostream>

double conjgrad::energy_in_grad_dir(double t) {
    std::vector<double> config;
    cl->config_to_vec(&config);
    cl->compute_energy_gradient(&grad);

    for(int i = 0; i<config.size(); i++){
        config[i] = config[i] + t * grad[i];
    }

    cluster cl2 = cluster(config);

    return cl2.compute_energy();
}

conjgrad::conjgrad(cluster *given) : cl(given) {
    cl->compute_energy_gradient(&grad);
}

void conjgrad::dosomehting() {
    LOG(energy_in_grad_dir(10));
}
