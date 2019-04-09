//
// Created by jan on 4/5/19.
//

#include <iostream>
#include <cmath>
#include "dipole.h"
#include "misc.h"

dipole::dipole() {
    set_r_random();
    set_m_random();
}

void dipole::set_r(double x, double y, double z) {
    r = {x, y, z};
}

void dipole::set_m(double phi, double theta, double m_length) {
    if(phi < 0 || phi >= 2 * M_PI){
        std::cout << "Phi has to be in [0, 2Pi)" << std::endl;
        return;
    }

    if(theta < 0 || theta >= M_PI){
        std::cout << "Phi has to be in [0, Pi)" << std::endl;
        return;
    }


    m = {
            m_length * (sin(theta) * cos(phi)),
            m_length * (sin(theta) * sin(phi)),
            m_length * (cos(theta))
    };

}

void dipole::print() {
    std::cout << r[0] << " " << r[1] << " " << r[2] << " " << m[0] << " " << m[1] << " " << m[2] << std::endl;

}

double dipole::distance_to(const dipole& v) {
    double x = this->r[0] - v.r[0];
    double y = this->r[1] - v.r[1];
    double z = this->r[2] - v.r[2];

    return sqrt(x*x + y*y + z*z);
}

std::vector<double> dipole::get_m() {
    return this->m;

}
std::vector<double> dipole::get_r() {
    return this->r;

}

std::vector<double> dipole::vector_to(const dipole& v) {
    std::vector<double> ret;
    ret.reserve(3);
    ret.push_back(v.r[0] - this->r[0]);
    ret.push_back(v.r[1] - this->r[1]);
    ret.push_back(v.r[2] - this->r[2]);

    return ret;
}

void dipole::set_r_random() {
    double x = DIPOLE_MAX_RANDOM_R * misc::random_simple();
    double y = DIPOLE_MAX_RANDOM_R * misc::random_simple();
    double z = DIPOLE_MAX_RANDOM_R * misc::random_simple();
    set_r(x, y, z);
}

void dipole::set_m_random() {
    double phi = 2*M_PI * misc::random_simple();
    double theta = M_PI * misc::random_simple();
    set_m(phi, theta);
}

dipole::dipole(const std::vector<double> &r) {
    set_m_random();
    this->r = r;
}

dipole::dipole(double x, double y, double z){
    set_m_random();
    set_r(x, y, z);
}

dipole::dipole(double x, double y, double z, double phi, double theta) {
    set_r(x, y, z);
    set_m(phi, theta);

}

