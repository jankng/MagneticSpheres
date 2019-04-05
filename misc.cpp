//
// Created by jan on 4/5/19.
//

#include "misc.h"
#include <cstdlib>
#include <iostream>

double misc::SimpleRandom() {
    return (double) rand() / RAND_MAX;
}

double misc::DotProduct(std::vector<double> u, std::vector<double> v) {
    if(u.size() != v.size()){
        std::cout << "Vectors are not of same dimension!" << std::endl;
        return 0;
    }

    double ret = 0;
    for(int i = 0; i<u.size(); i++){
        ret += u[i] * v[i];
    }

    return ret;
}
