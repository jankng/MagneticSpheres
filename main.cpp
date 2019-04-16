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

void minimize(int n){
    double e_min = 0;

    // TODO CHANGE 2 TO 250
    for(int i = 0; i<250; i++){
        std::cout << "Thread " << n << ", i=" << i << std::endl;
        metropolis candidate(8);
        candidate.start_siman();

        double e_candidate = candidate.get_cluster()->compute_energy();
        if(e_candidate < e_min){
            e_min = e_candidate;
            std::string filename = "t" + std::to_string(n) + "i" + std::to_string(i) +
                    "e" + std::to_string(e_min) + ".txt";
            candidate.get_cluster()->write_to_file(filename);
        }
    }

    std::cout << "Thread " << n << " ended." << std::endl;
}


int main() {
    misc::setup_static_rng();


    std::vector<dipole> ds;
    ds.reserve(8);
    for(int i = 0; i<8; i++){
        ds.emplace_back(dipole(0, 0, i, 0, 0));
    }
    cluster cl(16);
    std::cout << cl.compute_energy() << std::endl;


/*
    std::cout << "stepsize5" << std::endl;
    std::cout << "starting threads..." << std::endl;
    std::thread one(minimize, 1);
    std::thread two(minimize, 2);
    std::thread three(minimize, 3);
    //std::thread four(minimize, 4);

    one.join();
    two.join();
    three.join();
    //four.join();

    std::cout << "threads joined." << std::endl;

    gsl_rng* test = misc::get_static_rng();
    std::cout << misc::random_simple() << std::endl;
    std::cout << gsl_rng_uniform(test) << std::endl;
    */

    metropolis test(&cl);
    test.enable_verbose_mode();
    test.start_siman();

    test.get_cluster()->print();

    misc::delete_static_rng();
    return 0;
}