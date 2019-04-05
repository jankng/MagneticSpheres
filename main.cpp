#include <iostream>
#include <cstdlib>
#include "Dipole.h"
#include "Cluster.h"

int main() {
    srand(time(NULL));

    int N = 1000;
    Cluster c(N);
    for(int i = 0; i<N; i++){
        Dipole* d = c.GetDipole(i);
        d->SetM(0, 0);
        d->SetR(0, 0, i);
    }

    double e = c.ComputeEnergy();
    std::cout << e << std::endl;

    std::cout << "Hello, World!" << std::endl;
    return 0;
}