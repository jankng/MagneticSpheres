// C++ Standard Library
#include <iostream>
#include <cmath>
#include <cstring>
#include <thread>
#include <string>

// GSL Headers
#include <gsl/gsl_randist.h>
#include <iomanip>

// Project Headers
#include "dipole.h"
#include "cluster.h"
#include "metropolis.h"
#include "misc.h"
#include "gsledits.h"

void minimize(metropolis* m){
    m->doStuff();
}


int main() {
    misc::setup_static_rng();

/*
    std::vector<dipole> ds;
    ds.reserve(8);
    for(int i = 0; i<8; i++){
        ds.emplace_back(dipole(0, 0, i, 0, 0));
    }
    cluster test(ds);
    std::cout << test.compute_energy() << std::endl;
*/

    metropolis x(8);
    metropolis y(8);

    std::thread one(minimize, &x);
    std::thread two(minimize, &y);

    one.join();
    two.join();

    auto* a = x.get_cluster();
    auto* b = y.get_cluster();

    a->print();
    std::cout << a->compute_energy() << std::endl;
    b->print();
    std::cout << b->compute_energy() << std::endl;

    if(a->compute_energy() < b->compute_energy())
        a->write_to_file();
    else
        b->write_to_file();

    misc::delete_static_rng();
    return 0;
}