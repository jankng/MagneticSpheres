//
// Created by jan on 4/5/19.
//

#include <iostream>
#include <cmath>
#include "Dipole.h"
#include "misc.h"

Dipole::Dipole() {
    r.reserve(3);
    m.reserve(3);
    for(int i = 0; i<3; i++){
        r.push_back(0);
        m.push_back(0);
    }

    double x = 10 * misc::SimpleRandom();
    double y = 10 * misc::SimpleRandom();
    double z = 10 * misc::SimpleRandom();
    double phi = 2*M_PI * misc::SimpleRandom();
    double theta = M_PI * misc::SimpleRandom();

    this->SetR(x, y, z);
    this->SetM(phi, theta);

}

void Dipole::SetR(double x, double y, double z) {
    this->r[0] = x;
    this->r[1] = y;
    this->r[2] = z;

}

void Dipole::SetM(double phi, double theta, double m) {
    if(phi < 0 || phi >= 2 * M_PI){
        std::cout << "Phi has to be in [0, 2Pi)" << std::endl;
        return;
    }

    if(theta < 0 || theta >= M_PI){
        std::cout << "Phi has to be in [0, Pi)" << std::endl;
        return;
    }

    this->m[0] = m * (sin(theta) * cos(phi));
    this->m[1] = m * (sin(theta) * sin(phi));
    this->m[2] = m * (cos(theta));

}

void Dipole::Print() {
    std::cout << r[0] << " " << r[1] << " " << r[2] << " " << m[0] << " " << m[1] << " " << m[2] << std::endl;

}

double Dipole::DistanceTo(Dipole a) {
    double x = this->r[0] - a.r[0];
    double y = this->r[1] - a.r[1];
    double z = this->r[2] - a.r[2];

    return sqrt(x*x + y*y + z*z);
}

std::vector<double> Dipole::GetM() {
    return this->m;

}
std::vector<double> Dipole::GetR() {
    return this->r;

}

std::vector<double> Dipole::VectorTo(Dipole v) {
    std::vector<double> ret;
    ret.reserve(3);
    ret.push_back(v.r[0] - this->r[0]);
    ret.push_back(v.r[1] - this->r[1]);
    ret.push_back(v.r[2] - this->r[2]);

    return ret;
}
