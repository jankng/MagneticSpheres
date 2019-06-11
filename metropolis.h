//
// Created by jkoenig on 09/04/19.
//

#ifndef MAGNETICSPHERES_METROPOLIS_H
#define MAGNETICSPHERES_METROPOLIS_H

#include <memory>

#include <gsl/gsl_siman.h>

#include "cluster.h"

#define ITERS_FIXED_T 1000
#define STEP_SIZE 0.1
#define INITIAL_T 10
#define MU_T 1.001
#define T_MIN 0.0005

class metropolis {
private:
    bool owns_cluster;
    gsl_rng* r;
    cluster* cl;
    int cluster_size;
    gsl_siman_params_t params;
    bool verbose;


    // helper functions needed for GSL
    static void copy_func(void *source, void *dest);
    static void * copy_constructor(void *xp);
    static void destroy_state(void *xp);
    static double energy_func(void *xp);
    static double mod1(double x);
    static void print_state(void* xp);

    // step functions
    static void take_step(const gsl_rng *r, void *xp, double step_size);
    static void take_step_no_constraints(const gsl_rng *r, void *xp, double step_size);
    static void take_step_fixed_ends(const gsl_rng *r, void *xp, double step_size);
    static void take_step_rigid_body(const gsl_rng *r, void *xp, double step_size);



public:
    explicit metropolis(int n); // generates random cluster with n spheres
    explicit metropolis(cluster* cluster_given); // takes pregenerated cluster
    metropolis(cluster* cluster_given, gsl_siman_params_t params_given);
    ~metropolis();

    cluster* get_cluster(){return cl;}

    void start_siman();

    void enable_verbose_mode(){verbose = true;}
    void disable_verbose_mode(){verbose = false;}
};


#endif //MAGNETICSPHERES_METROPOLIS_H
