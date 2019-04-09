#include <iostream>
#include <cstdlib>
#include "Dipole.h"
#include "Cluster.h"

int main() {
    srand(time(NULL));

    int N = 8;
    Cluster c(N);
    for(int i = 0; i<N; i++){
        Dipole* d = c.GetDipole(i);
        d->SetR(0, 0, i);
    }

    double e = c.ComputeEnergy();
    std::cout << "Energy of the configuration in the beginning: " << e << std::endl;

    for(int i = 0; i<10000; i++){
        c.MetropolisStep();
    }

    std::cout << "Energy of the configuration after 1000 Steps: " << c.ComputeEnergy() << std::endl;

    c.Print();

    return 0;
}