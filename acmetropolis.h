//
// Created by jan on 7/9/19.
//

#ifndef MAGNETICSPHERES_ACMETROPOLIS_H
#define MAGNETICSPHERES_ACMETROPOLIS_H


#include <memory>
#include <algorithm>

#include <gsl/gsl_siman.h>

#include "cluster.h"
#include "dipole.h"

//Best parameters for simultaneous

#define ITERS_FIXED_T 1000
#define STEP_SIZE 64
#define INITIAL_T 10
#define MU_T 1.001
#define T_MIN 0.0005


//Parameters for independent
/*
#define ITERS_FIXED_T 1000
#define STEP_SIZE 0.1
#define INITIAL_T 0.001
#define MU_T 1.001
#define T_MIN 0.0005
 */

struct acconfig{
    int n;
    std::vector<double> angs;
    std::vector<double> dips;
};


class acmetropolis {
private:
    bool owns_cluster;
    gsl_rng* r;
    acconfig* cfg;
    int cluster_size;
    gsl_siman_params_t params;
    bool verbose;


    // helper functions needed for GSL
    static void copy_func(void *source, void *dest);
    static void * copy_constructor(void *xp);
    static void destroy_state(void *xp);
    static double mod1(double x);
    static void print_state(void* xp);
    static double energy_func(void *xp);

    // step functions
    static void take_step(const gsl_rng *r, void *xp, double step_size);




public:
    explicit acmetropolis(int n); // generates random cluster with n spheres
    acmetropolis(int n, double g, double step_size, bool symmetric_dipoles, bool constraints);
    ~acmetropolis();

    acconfig* get_cluster(){return cfg;}

    void start_siman();

    void enable_verbose_mode(){verbose = true;}
    void disable_verbose_mode(){verbose = false;}

    // conversion to cluster
    static void ac_to_cluster(acconfig& cfg, cluster* cl);
};


#endif //MAGNETICSPHERES_ACMETROPOLIS_H
