#include <iostream>
#include "dipole.h"
#include "cluster.h"
#include "metropolis.h"
#include "misc.h"

double efunc(void* xp){
    return ((cluster*) xp)->compute_energy();
}

int main() {

    misc::new_rng();

    std::vector<dipole> config = {};
    int N = 1000;
    for(int i = 0; i<N; i++){
        dipole d(0, 0, i, 0, 0);
        config.emplace_back(d);
    }

    cluster ideal(config);
    cluster rnd(N, chain);
    cluster rndc(N);

    std::cout << "Energy of ideal chain of length " << N << ": " << ideal.compute_energy() << std::endl;
    std::cout << "Energy of random chain of length " << N << ": " << rnd.compute_energy() << std::endl;
    std::cout << "Energy of random cluster of size " << N << ": " << rndc.compute_energy() << std::endl;

    rndc.print();


    misc::delete_rng();
    return 0;
}