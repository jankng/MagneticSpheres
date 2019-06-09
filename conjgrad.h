//
// Created by jkoenig on 29/04/19.
//

#ifndef MAGNETICSPHERES_CONJGRAD_H
#define MAGNETICSPHERES_CONJGRAD_H

#include "definitions.h"
#include "cluster.h"
#include "dipole.h"

#include <gsl/gsl_multimin.h>

class conjgrad {
private:
    std::vector<double> config;
    std::vector<double> grad;

    void compute_gradient(int i, int mode);
    double compute_energy();
    double compute_energy_in_direction(double t, const std::vector<double>& dir);
    double minimize_in_direction(const std::vector<double>& dir);
    void go_in_direction(double t, const std::vector<double>& dir);

public:
    explicit conjgrad(std::vector<double>* conf);
    explicit conjgrad(cluster* cl);
    explicit conjgrad(int n);

    void minimize_simultaneous();

    void print_energy_in_direction(std::vector<double>* dir);

    void dosomething();
};


#endif //MAGNETICSPHERES_CONJGRAD_H
