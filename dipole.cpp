//
// Created by jan on 4/5/19.
//

#include <iostream>
#include <cmath>
#include <sstream>
#include "dipole.h"
#include "misc.h"

dipole::dipole() {
    set_r_random();
    set_m_random();
}

void dipole::set_r(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;
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

    this->phi = phi;
    this->theta = theta;

    m = {
            m_length * (sin(theta) * cos(phi)),
            m_length * (sin(theta) * sin(phi)),
            m_length * (cos(theta))
    };

}

void dipole::print() {
    std::cout << to_string(' ') << std::endl;

}

double dipole::distance_to(const dipole& v) {
    double x1 = this->r[0] - v.r[0];
    double y1 = this->r[1] - v.r[1];
    double z1 = this->r[2] - v.r[2];

    return sqrt(x1*x1 + y1*y1 + z1*z1);
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
    double x1 = (int) (DIPOLE_MAX_RANDOM_R * misc::random_simple());
    double y1 = (int) (DIPOLE_MAX_RANDOM_R * misc::random_simple());
    double z1 = (int) (DIPOLE_MAX_RANDOM_R * misc::random_simple());
    set_r(x1, y1, z1);
}

void dipole::set_m_random() {
    double phi1 = 2*M_PI * misc::random_simple();
    double theta1 = acos(1 - 2 * misc::random_simple());
    set_m(phi1, theta1);
}

dipole::dipole(const std::vector<double> &r) {
    set_r(r[0], r[1], r[2]);
    set_m(r[3], r[4]);
}

dipole::dipole(double x, double y, double z){
    set_m_random();
    set_r(x, y, z);
}

dipole::dipole(double x, double y, double z, double phi, double theta) {
    set_r(x, y, z);
    set_m(phi, theta);

}

dipole dipole::dipole_in_direction(const std::vector<double> &dir) {
    double new_x = this->x + dir[0];
    double new_y = this->y + dir[1];
    double new_z = this->z + dir[2];
    double new_phi = this->phi + dir[3];
    double new_theta = this->theta + dir[4];

    std::vector<double> coords = {new_x, new_y, new_z, new_phi, new_theta};

    return dipole(coords);
}

std::vector<double> dipole::get_angles() {
    std::vector<double> ang = {phi, theta};
    return ang;
}

bool dipole::is_in_bounds() {
    bool flag = true;


    if(r[0] < 0 || r[0] > 50)
        flag = false;
    if(r[1] < 0 || r[1] > 50)
        flag = false;
    if(r[2] < 0 || r[2] > 50)
        flag = false;


    return flag;
}

std::string dipole::to_string(char sep) {
    std::stringstream ret;
    ret << r[0] << sep << r[1] << sep << r[2] << sep
        << m[0] << sep << m[1] << sep << m[2];

    return ret.str();
}

dipole::dipole(double x, double y, double z, double m1, double m2, double m3) {
    this->r = {x, y, z};
    this->m = {m1, m2, m3};

}

