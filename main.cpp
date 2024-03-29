// C++ Standard Library
#include <iostream>
#include <cmath>
#include <cstring>
#include <thread>
#include <string>
#include <iomanip>

// GSL Headers
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>

// Project Headers
#include "definitions.h"
#include "dipole.h"
#include "cluster.h"
#include "metropolis.h"
#include "acmetropolis.h"
#include "misc.h"
#include "gsledits.h"
#include "conjgrad.h"

void startMetropolis(int id){
    int n = 0;
    bool symmetric = true;
    bool constraints = true;

    //determine parameters
    gsl_siman_params_t params_given;
    switch(id) {
        case 1:
            n=8;
            break;
        case 2:
            n=9;
            break;
        case 3:
            n = 10;
            break;
        default:
            n = 5;
    }

    //start iterations
    for (int i = 0; i < 1000; i++) {
        acmetropolis met(n, 100, 1, n, symmetric, constraints);
        met.start_siman();
        //double energy = met.compute_energy();
        int G = (int) 1;
        std::string filename =
                "n" + std::to_string(n) + "g" + std::to_string(G) + "i" + std::to_string(i) + "d5t.txt";

        met.write_to_file(filename);
        std::cout << "Thread " << id << " n" << n << " i" << i << " ended." << std::endl;
    }

    std::cout << "Thread " << id << " ended." << std::endl;
}

void startMetropolisThreads(){
    std::thread t1(startMetropolis, 1);
    std::thread t2(startMetropolis, 2);
    //std::thread t3(startMetropolis, 3);

    t1.join();
    t2.join();
    //t3.join();

    std::cout << "Threads joined successfully";
}

cluster* make_perfect_chain(int n){
    std::vector<dipole> dps;
    for(int i = 0; i<n; i++){
        dipole d(i, 0, 10, 0, M_PI / 2.0);
        dps.emplace_back(d);
    }

    cluster* c = new cluster(dps);
    return c;
}

cluster* make_random_chain(int n){
    std::vector<dipole> dps;
    for(int i = 0; i<n; i++){
        double angle = misc::random_simple() * 2.0 * M_PI;

        dipole d(i, 0, 5, cos(angle), 0, sin(angle));
        dps.emplace_back(d);
    }

    cluster* c = new cluster(dps);
    return c;
}

int main() {
    misc::setup_static_rng();


    //acmetropolis(int n, double h, double g, double step_size, bool symmetric_dipoles, bool constraints);
    //acmetropolis m(8, 8, 1, 1, true, false);
    //m.enable_verbose_mode();
    //m.start_siman();
    //std::cout << m.to_string() << std::endl;


    // SET OUTPUT FILE DIRECTORY BEFORE UNCOMMENTING
    //startMetropolisThreads();

    cluster* cl = make_perfect_chain(8);
    conjgrad cg(cl);
    cg.minimize_simultaneous();

    misc::delete_static_rng();
    return 0;
}