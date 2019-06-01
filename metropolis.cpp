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
        owns_cluster(true), r(misc::get_static_rng()), cl(new cluster(8)), cluster_size(n),
        params{0, ITERS_FIXED_T, STEP_SIZE, 1.0, INITIAL_T, MU_T, T_MIN},
        verbose(false){}

metropolis::metropolis(cluster* cluster_given):
        owns_cluster(false), r(misc::get_static_rng()), cl(cluster_given), cluster_size(cl->get_size()),
        params{0, ITERS_FIXED_T, STEP_SIZE, 1.0, INITIAL_T, MU_T, T_MIN},
        verbose(false){}

metropolis::metropolis(cluster* cluster_given, gsl_siman_params_t params_given):
        owns_cluster(false), r(misc::get_static_rng()), cl(cluster_given), cluster_size(cl->get_size()),
        params(params_given),
        verbose(false){}


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

// original version
/*
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

        // add new coords to new config
        new_config.emplace_back(dipole(new_coords));

        j++;
    }

    *c = cluster(new_config);

}
*/

// take stept with fixed first and last dipole
void metropolis::take_step(const gsl_rng *r, void *xp, double step_size){
    // set up boilerplate
    auto* c = (cluster*) xp;
    double d_size = gsl_rng_uniform(r) * step_size;

    // delta of coordinates
    int n_of_coords = c->get_size() * 5 - 10;
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
    // start change
    dipole start = *(c->get_dipole_by_ref(0));
    new_config.emplace_back(start);
    int j = 1;
    //end change
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

        // add new coords to new config
        new_config.emplace_back(dipole(new_coords));

        j++;
    }


    // start change
    dipole end = *(c->get_dipole_by_ref(j));
    new_config.emplace_back(end);
    //end change


    *c = cluster(new_config);

}

void metropolis::print_state(void* xp){
    auto* c = (cluster*) xp;
    std::cout << " x=" << "somewhere" << " ";
}

void metropolis::start_siman(){

    // tries, steps/temp, max step, k, temp init, temp cooldown, temp min
    // best so far: gsl_siman_params_t params = {0, 1000, 5.0, 1.0, 10, 1.01, 0.001};
    if(verbose)
        gsl_edits::gsl_siman_solve(r, cl,
                                   metropolis::energy_func, metropolis::take_step, metropolis::print_state,
                                   metropolis::copy_func, metropolis::copy_constructor, metropolis::destroy_state,
                                   sizeof(cluster), params);
    else
        gsl_edits::gsl_siman_solve(r, cl,
                                   metropolis::energy_func, metropolis::take_step, nullptr,
                                   metropolis::copy_func, metropolis::copy_constructor, metropolis::destroy_state,
                                   sizeof(cluster), params);


}

metropolis::~metropolis() {
    if(owns_cluster)
        delete this->cl;
}
