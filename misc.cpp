//
// Created by jan on 4/5/19.
//

#include "misc.h"
#include <iostream>
#include <gsl/gsl_rng.h>

gsl_rng* misc::class_r = nullptr;

// TODO better random_simple() implementation
double misc::random_simple() {
    return gsl_rng_uniform(class_r);
}

double misc::dot_product(const std::vector<double>& u, const std::vector<double>& v) {
    if(u.size() != v.size()){
        std::cout << "Vectors are not of same dimension!" << std::endl;
        return 0;
    }

    double ret = 0;
    for(int i = 0; i<u.size(); i++){
        ret += u[i] * v[i];
    }

    return ret;
}

void misc::new_rng() {
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    class_r = gsl_rng_alloc (T);

    // seed with time
    gsl_rng_set(class_r, time(nullptr));

}

void misc::delete_rng() {
    gsl_rng_free(class_r);
}