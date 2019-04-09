//
// Created by jan on 4/5/19.
//

#ifndef MAGNETICSPHERES_MISC_H
#define MAGNETICSPHERES_MISC_H

#include <vector>

class misc {
public:

    static double dot_product(const std::vector<double>& u, const std::vector<double>& v);

    // functions for random numbers


    // uniformly distributed number in [0, 1]
    static double random_simple();
};


#endif //MAGNETICSPHERES_MISC_H
