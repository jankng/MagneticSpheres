#include <iostream>
#include <cstdlib>
#include "dipole.h"
#include "cluster.h"
#include "metropolis.h"
#include <gsl/gsl_sf_bessel.h>

int main() {
    srand(time(nullptr));

    std::vector<dipole> config = {};
    int N = 1000;
    for(int i = 0; i<N; i++){
        dipole d(0, 0, i, 0, 0);
        config.emplace_back(d);
    }

    cluster c(config);
    c.print();

    std::cout << "Energy of c: " << c.compute_energy() << std::endl;

    std::cout << gsl_sf_bessel_J0(5.0) << std::endl;
    

    return 0;
}