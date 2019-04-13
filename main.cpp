#include <iostream>
#include "dipole.h"
#include "cluster.h"
#include "metropolis.h"
#include "misc.h"
#include <gsl/gsl_siman.h>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <cstring>


double E1(void* xp){
    cluster* c = (cluster*) xp;
    return c->compute_energy();
}

double mod1(double x){
    x = fmod(x, 0.99);
    if(x < 0)
        x = x + 1;
    return x;
}

void S1(const gsl_rng* r, void* xp, double step_size){
    // set up boilerplate
    cluster* c = (cluster*) xp;
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

void P1(void* xp){
    std::cout << " x=" << "somewhere" << " ";
}

int main() {

    //std::cout << acos(100) << std::endl;

    misc::new_rng();

    /*
    std::vector<dipole> ds;
    ds.reserve(1000);
    for(int i = 0; i<1000; i++){
        ds.emplace_back(dipole(0, 0, i, 0, 0));
    }
    cluster test(ds);
    std::cout << test.compute_energy() << std::endl;
    */

    const gsl_rng_type* T;
    gsl_rng* r;

    gsl_rng_env_setup();
    T = gsl_rng_mt19937;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(nullptr));

    cluster x_initial(8);
    gsl_siman_params_t params = {0, 1000, 5.0, 1.0, 0.01, 1.01, 0.001};
    gsl_siman_solve(r, &x_initial, E1, S1, nullptr, P1, nullptr, nullptr, nullptr, sizeof(cluster), params);
    gsl_rng_free(r);

    x_initial.print();

    misc::delete_rng();
    return 0;
}