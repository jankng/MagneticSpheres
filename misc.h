//
// Created by jan on 4/5/19.
//

#ifndef MAGNETICSPHERES_MISC_H
#define MAGNETICSPHERES_MISC_H

#include <vector>
#include <string>
#include <gsl/gsl_rng.h>

class misc {
private:
    // rng variables
    static gsl_rng* static_r;
    static bool r_is_init;

public:
    // get time string
    static std::string get_time();

    // math functions
    static double dot_product(const std::vector<double>& u, const std::vector<double>& v);
    static double mod1(double d);
    static double modn(double d, double n);

    // functions for random numbers
    static void setup_static_rng();
    static void delete_static_rng();
    static gsl_rng* get_static_rng();
    static bool rng_is_initialized();

    // uniformly distributed number in [0, 1]
    static double random_simple();
};
#endif //MAGNETICSPHERES_MISC_H
