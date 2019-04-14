//
// Created by jkoenig on 09/04/19.
//

#ifndef MAGNETICSPHERES_METROPOLIS_H
#define MAGNETICSPHERES_METROPOLIS_H

#include <memory>

#include <gsl/gsl_siman.h>

#include "cluster.h"

#define DEFAULT_BETA_START 1

class metropolis {
private:
    cluster* cl;
    int cluster_size;
    double cluster_energy;
    int step_count, attempt_count;
    gsl_rng* r;

    // helper functions needed for GSL
    static void copy_func(void *source, void *dest);
    static void * copy_constructor(void *xp);
    static void destroy_state(void *xp);
    static double energy_func(void *xp);
    static double mod1(double x);
    static void take_step(const gsl_rng *r, void *xp, double step_size);
    static void print_state(void* xp);



public:
    explicit metropolis(int n); // generates random cluster with n spheres
    ~metropolis();

    double get_energy(){return cluster_energy;}

    void doStuff();
};


#endif //MAGNETICSPHERES_METROPOLIS_H
