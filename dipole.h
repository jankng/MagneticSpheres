//
// Created by jan on 4/5/19.
//

#ifndef MAGNETICSPHERES_DIPOLE_H
#define MAGNETICSPHERES_DIPOLE_H

#include <vector>

#define MAX_RANDOM_R 1
#define DEFAULT_M 1

class dipole {
private:
    //double x, y, z, phi, theta;
    std::vector<double> r, m;

    //setters for internal use in constructor only
    void set_r(double x, double y, double z);
    void set_r_random();
    void set_m(double phi, double theta, double m_length = DEFAULT_M);
    void set_m_random();
public:
    // generates dipole with random r and m
    dipole();
    // generates dipole with random m
    dipole(double x, double y, double z);
    explicit dipole(std::vector<double> r); // TODO look up why 'explicit'

    //getters
    std::vector<double> get_r();
    std::vector<double> get_m();

    // mathematical properties
    double distance_to(dipole v);
    std::vector<double> vector_to(dipole v);

    //output
    void print();

};


#endif //MAGNETICSPHERES_DIPOLE_H
