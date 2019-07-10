//
// Created by jan on 7/9/19.
//

#include "acmetropolis.h"

#include <cmath>
#include <iostream>
#include <gsl/gsl_randist.h>

#include "misc.h"
#include "gsledits.h"

acmetropolis::acmetropolis(int n):
        owns_cluster(true), r(misc::get_static_rng()), cluster_size(n),
        params{0, ITERS_FIXED_T, STEP_SIZE, 1.0, INITIAL_T, MU_T, T_MIN},
        verbose(false){
    this->cfg = new acconfig;
    cfg->angs.reserve(n/2 + n%2 - 2);
    cfg->dips.reserve(n);
    cfg->n = n;

    for(int i = 0; i<n; i++){
        cfg->dips.emplace_back(0);

        if(i < n/2 + n%2 - 1)
            cfg->angs.emplace_back(0);
    }
}


void acmetropolis::copy_func(void *source, void *dest){
    auto* s = (acconfig*) source;
    auto* d = (acconfig*) dest;

    *d = *s;
}

void * acmetropolis::copy_constructor(void *xp){
    auto* x = (acconfig*) xp;
    auto* y = new acconfig;

    *y = *x;

    return y;
}

void acmetropolis::destroy_state(void *xp){
    auto* x = (acconfig*) xp;
    delete x;
}

double acmetropolis::energy_func(void *xp){
    auto* cfg = (acconfig*) xp;
    cluster cl;
    ac_to_cluster(*cfg, &cl);
    return cl.compute_energy_for_metropolis();
}

double acmetropolis::mod1(double x){
    x = fmod(x, 0.99);
    if(x < 0)
        x = x + 1;
    return x;
}

void acmetropolis::take_step(const gsl_rng *r, void *xp, double step_size){
    // set up boilerplate
    auto* cfg = (acconfig*) xp;
    double d_size = gsl_rng_uniform(r) * step_size;

    // delta of coordinates
    int n_of_coords = cfg->angs.size() + cfg->dips.size();
    std::vector<double> delta_coords;
    delta_coords.reserve(n_of_coords);

    //new config will be stored here
    acconfig new_cfg = {cfg->n};
    new_cfg.angs.reserve(n_of_coords);
    new_cfg.dips.reserve(n_of_coords);

    // generate distance vector and norm square
    double d_squared = 0;
    for(int i = 0; i<n_of_coords; i++) {
        double rnd = gsl_ran_gaussian(r, 1.0);
        delta_coords.emplace_back(rnd);
        d_squared += pow(rnd, 2);
    }

    // normalize distance vector
    double d = sqrt(d_squared);
    for(int i = 0; i<n_of_coords; i++){
        delta_coords[i] = (delta_coords[i] * d_size) / d;
    }

    // generate step
    for(int i = 0; i<cfg->angs.size(); i++){

        // add new coords to new config
        new_cfg.angs.emplace_back(cfg->angs[i] + delta_coords[i]);
    }

    for(int i = cfg->angs.size(); i<n_of_coords; i++){

        // add new coords to new config
        new_cfg.dips.emplace_back(cfg->dips[i] + delta_coords[i]);
    }

    //for(int i = 0; i<new_cfg.dips.size() / 2; i++){
        //new_cfg.dips[new_cfg.dips.size()- 1 - i] = -1.0*new_cfg.dips[i];
    //}

    *cfg = new_cfg;

}

void acmetropolis::print_state(void* xp){
    auto* c = (cluster*) xp;
    std::cout << " x=" << "somewhere" << " ";
}

void acmetropolis::start_siman(){

    // tries, steps/temp, max step, k, temp init, temp cooldown, temp min
    // best so far: gsl_siman_params_t params = {0, 1000, 5.0, 1.0, 10, 1.01, 0.001};
    if(verbose)
        gsl_edits::gsl_siman_solve(r, cfg,
                                   acmetropolis::energy_func, acmetropolis::take_step, acmetropolis::print_state,
                                   acmetropolis::copy_func, acmetropolis::copy_constructor, acmetropolis::destroy_state,
                                   sizeof(cluster), params);
    else
        gsl_edits::gsl_siman_solve(r, cfg,
                                   acmetropolis::energy_func, acmetropolis::take_step, nullptr,
                                   acmetropolis::copy_func, acmetropolis::copy_constructor, acmetropolis::destroy_state,
                                   sizeof(cluster), params);


}

acmetropolis::~acmetropolis() {
    if(owns_cluster)
        delete this->cfg;
}

void acmetropolis::ac_to_cluster(acconfig& cfg, cluster* cl) {
    std::vector<dipole> dp_config;
    dp_config.reserve(cfg.n);

    double x = 0;
    double y = 0;
    double z = 10;

    // first half
    for(int i = 0; i<cfg.n / 2; i++){
        dp_config.emplace_back(dipole(
                x, y, z,
                cos(cfg.dips[i]), 0, sin(cfg.dips[i])
                ));

        if(i + 1 < cfg.n / 2) {
            x += cos(cfg.angs[i]);
            z += sin(cfg.angs[i]);
        }
    }

    // if length is odd
    if(cfg.n % 2 == 1){
        x += cos(cfg.angs[cfg.angs.size() - 1]);
        z += sin(cfg.angs[cfg.angs.size() - 1]);

        dp_config.emplace_back(dipole(
                x, y, z,
                cos(cfg.dips[cfg.n/2]), 0, sin(cfg.dips[cfg.n/2])
        ));

        x += cos(-1.0*cfg.angs[cfg.angs.size() - 1]);
        z += sin(-1.0*cfg.angs[cfg.angs.size() - 1]);
    } else{
        x += 1.0;
    }

    std::reverse(cfg.angs.begin(), cfg.angs.end());

    //second half
    for(int i = cfg.n/2 + cfg.n%2; i<cfg.n; i++){
        dp_config.emplace_back(dipole(
                x, y, z,
                cos(cfg.dips[i]), 0, sin(cfg.dips[i])
        ));

        if(i+1 < cfg.n) {
            x += cos(-1.0*cfg.angs[i - cfg.n/2]);
            z += sin(-1.0*cfg.angs[i - cfg.n/2]);
        }
    }

    std::reverse(cfg.angs.begin(), cfg.angs.end());


    *cl = cluster(dp_config);

}



