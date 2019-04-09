//
// Created by jkoenig on 09/04/19.
//

#include "metropolis.h"
#include "misc.h"
#include <cmath>

metropolis::metropolis(int n) {
    cluster_size = n;
    step_count = attempt_count = 0;
    beta = DEFAULT_BETA_START;

    std::shared_ptr<cluster> c = std::make_shared<cluster>(n);
    cl = c;
}
