#include <iostream>
#include "dipole.h"
#include "cluster.h"
#include "metropolis.h"
#include "misc.h"
#include <gsl/gsl_siman.h>
#include <cmath>


double E1(void* xp){
    double x = *((double*) xp);
    return -1.0*sin(x) / x;
}

double M1(void *xp, void* yp){
    double x = *((double*) xp);
    double y = *((double*) yp);
    return fabs(x - y);
}

void S1(const gsl_rng* r, void* xp, double step_size){
    double* xptr = (double*) xp;
    double oldx = *((double*) xp);
    double newx;

    double u = misc::random_simple();
    newx = oldx + (1.0 - 2.0*u) * step_size;

    *xptr = newx;
}

void P1(void* xp){
    std::cout << " x=" << *((double*) xp) << " ";
}

int main() {

    misc::new_rng();

    const gsl_rng_type* T;
    gsl_rng* r;

    gsl_rng_env_setup();
    T = gsl_rng_mt19937;
    r = gsl_rng_alloc(T);

    double x_initial = 3;
    gsl_siman_params_t params = {0, 1000, 1.0, 1.0, 0.008, 1.003, 2.0e-6};
    gsl_siman_solve(r, &x_initial, E1, S1, nullptr, P1, nullptr, nullptr, nullptr, sizeof(double), params);
    gsl_rng_free(r);

    misc::delete_rng();
    return 0;
}