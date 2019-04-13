//
// Created by jkoenig on 09/04/19.
//

#ifndef MAGNETICSPHERES_METROPOLIS_H
#define MAGNETICSPHERES_METROPOLIS_H

#include <memory>
#include "cluster.h"

#define DEFAULT_BETA_START 1

class metropolis {
private:
    std::shared_ptr<cluster> cl;
    int cluster_size;
    double beta, energy_current;
    int step_count, attempt_count;


public:
    explicit metropolis(int n); // generates random cluster with n spheres

    double get_energy(){return energy_current;}
};


#include <gsl/gsl_siman.h>

    void
    gsl_siman_solve_debug (const gsl_rng * r, void *x0_p, gsl_siman_Efunc_t Ef,
                     gsl_siman_step_t take_step,
                     gsl_siman_metric_t distance,
                     gsl_siman_print_t print_position,
                     gsl_siman_copy_t copyfunc,
                     gsl_siman_copy_construct_t copy_constructor,
                     gsl_siman_destroy_t destructor,
                     size_t element_size,
                     gsl_siman_params_t params);


#endif //MAGNETICSPHERES_METROPOLIS_H
