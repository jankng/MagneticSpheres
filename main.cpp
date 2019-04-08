#include <iostream>
#include <cstdlib>
#include "Dipole.h"
#include "Cluster.h"

int main() {
    srand(time(NULL));

    int N = 8;
    std::shared_ptr<Cluster> cluster = std::make_shared<Cluster>(N);
    for(int i = 0; i<N; i++){
        Dipole* d = cluster->GetDipole(i);
        d->SetR(0, 0, i);
    }

    double e = cluster->ComputeEnergy();
    std::cout << "Energy of the configuration in the beginning: " << e << std::endl;

    return 0;
}