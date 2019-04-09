#include <iostream>
#include <cstdlib>
#include "dipole.h"
#include "cluster.h"
#include "metropolis.h"

int main() {
    srand(time(NULL));

    cluster c(10, chain);
    c.print();

    metropolis m(10);

    return 0;
}