#include <iostream>
#include <cstdlib>
#include "dipole.h"
#include "cluster.h"
#include "Metropolis.h"

int main() {
    srand(time(NULL));

    cluster c(10, chain);
    c.print();

    return 0;
}