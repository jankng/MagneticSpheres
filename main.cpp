// C++ Standard Library
#include <iostream>
#include <cmath>
#include <cstring>
#include <thread>

// GSL Headers
#include <gsl/gsl_randist.h>

// Project Headers
#include "dipole.h"
#include "cluster.h"
#include "metropolis.h"
#include "misc.h"
#include "gsledits.h"


int main() {
    misc::setup_static_rng();

///*
    std::vector<dipole> ds;
    ds.reserve(8);
    for(int i = 0; i<8; i++){
        ds.emplace_back(dipole(0, 0, i, 0, 0));
    }
    cluster test(ds);
    std::cout << test.compute_energy() << std::endl;
//*/

    //std::thread one(doStuff);
    //std::thread two(doStuff);
    //std::thread three(doStuff);

    //one.join();
    //two.join();
    //three.join();

    metropolis min(8);
    min.doStuff();

    misc::delete_static_rng();
    return 0;
}