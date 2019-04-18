//
// Created by jan on 4/5/19.
//

#ifndef MAGNETICSPHERES_DIPOLE_H
#define MAGNETICSPHERES_DIPOLE_H

#include <vector>
#include <string>

#define DIPOLE_MAX_RANDOM_R 10
#define DIPOLE_DEFAULT_M 1

class dipole {
private:
    double x, y, z, phi, theta;
    std::vector<double> r, m;

    //setters for internal use in constructor only
    void set_r(double x, double y, double z);
    void set_r_random();
    void set_m(double phi, double theta, double m_length = DIPOLE_DEFAULT_M);
    void set_m_random();

public:
    // generates dipole with random r and m
    dipole();
    // generates dipole with random m
    dipole(double x, double y, double z);

    // creates dipole with all DOFs defined
    dipole(double x, double y, double z, double phi, double theta);
    explicit dipole(const std::vector<double>& coords);

    //getters
    std::vector<double> get_r();
    std::vector<double> get_m();
    std::vector<double> get_angles();

    // mathematical properties
    double distance_to(const dipole& v);
    std::vector<double> vector_to(const dipole& v);
    dipole dipole_in_direction(const std::vector<double>& dir);
    bool is_in_bounds();

    //output
    std::string to_string(char sep = ';');
    void print();

};


#endif //MAGNETICSPHERES_DIPOLE_H
