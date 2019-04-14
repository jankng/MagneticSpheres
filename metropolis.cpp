//
// Created by jkoenig on 09/04/19.
//
#include <cmath>
#include <iostream>

#include <gsl/gsl_randist.h>

#include "metropolis.h"
#include "misc.h"
#include "gsledits.h"

metropolis::metropolis(int n):
cl(new cluster(8)), cluster_size(n), cluster_energy(cl->compute_energy()), step_count(0), attempt_count(0),
r(misc::make_rng()){}


void metropolis::copy_func(void *source, void *dest){
    auto* s = (cluster*) source;
    auto* d = (cluster*) dest;

    *d = *s;
}

void * metropolis::copy_constructor(void *xp){
    auto* x = (cluster*) xp;
    auto* y = new cluster;

    *y = *x;
    return y;
}

void metropolis::destroy_state(void *xp){
    auto* x = (cluster*) xp;
    delete x;
}

double metropolis::energy_func(void *xp){
    auto* c = (cluster*) xp;
    return c->compute_energy();
}

double metropolis::mod1(double x){
    x = fmod(x, 0.99);
    if(x < 0)
        x = x + 1;
    return x;
}

void metropolis::take_step(const gsl_rng *r, void *xp, double step_size){
    // set up boilerplate
    auto* c = (cluster*) xp;
    double d_size = gsl_rng_uniform(r) * step_size;

    // delta of coordinates
    int n_of_coords = c->get_size() * 5;
    std::vector<double> delta_coords;
    delta_coords.reserve(n_of_coords);

    //new config will be stored here
    std::vector<dipole> new_config;
    new_config.reserve(c->get_size());

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
    int j = 0;
    for(int i = 0; i<n_of_coords; i+=5){
        // old coords and angs
        std::vector<double> xyz = c->get_dipole_by_ref(j)->get_r();
        std::vector<double> angs = c->get_dipole_by_ref(j)->get_angles();

        //get new angs
        double phi_gen = angs[0] / (2.0*M_PI) + delta_coords[i+3];
        phi_gen = mod1(phi_gen);

        double theta_gen = 0.5*(1.0-cos(angs[1])) + delta_coords[i+4];
        theta_gen = mod1(theta_gen);

        //set up new coords
        std::vector<double> new_coords = {
                xyz[0] + delta_coords[i+0],
                xyz[1] + delta_coords[i+1],
                xyz[2] + delta_coords[i+2],
                2.0*M_PI*phi_gen,
                acos(1.0-2.0*theta_gen)
        };

        //debug block
        if(new_coords[3] > 2*M_PI || new_coords[4] > M_PI)
            std::cout << "break" << std::endl;

        // add new coords to new config
        new_config.emplace_back(dipole(new_coords));

        j++;
    }

    *c = cluster(new_config);

    //std::cout << "" << std::endl;

}

void metropolis::print_state(void* xp){
    auto* c = (cluster*) xp;
    std::cout << " e=" << c->compute_energy() << " ";
}

void metropolis::doStuff(){

    // tries, steps/temp, max step, k, temp init, temp cooldown, temp min
    // best so far: gsl_siman_params_t params = {0, 1000, 5.0, 1.0, 10, 1.01, 0.001};
    gsl_siman_params_t params = {0, 1000, 1.0, 1.0, 10, 1.001, 0.001};

    gsl_edits::gsl_siman_solve(r, cl,
                               metropolis::energy_func, metropolis::take_step, nullptr, //metropolis::print_state,
                               metropolis::copy_func, metropolis::copy_constructor, metropolis::destroy_state,
                               sizeof(cluster), params);


    //cl->print();
    //std::cout << cl->compute_energy() << std::endl;

}

metropolis::~metropolis() {
    delete this->cl;
    gsl_rng_free(this->r);
}
