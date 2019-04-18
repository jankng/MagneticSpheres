//
// Created by jan on 4/5/19.
//

#include "misc.h"
#include <iostream>
#include <gsl/gsl_rng.h>
#include <sstream>
#include <iomanip>

gsl_rng* misc::static_r = nullptr;
bool misc::r_is_init = false;

double misc::random_simple() {
    return gsl_rng_uniform(static_r);
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

void misc::setup_static_rng() {
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_mt19937;
    static_r = gsl_rng_alloc (T);

    // seed with time
    gsl_rng_set(static_r, time(nullptr));

    r_is_init = true;

}

void misc::delete_static_rng() {
    gsl_rng_free(static_r);
    r_is_init = false;
}

// TODO DO NOT CLONE!
gsl_rng *misc::get_static_rng() {
    return static_r;
}

bool misc::rng_is_initialized() {
    return r_is_init;
}

std::string misc::get_time() {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::stringstream buffer;
    buffer << std::put_time(&tm, "%Y-%m-%d-%H-%M-%S");

    return buffer.str();
}
