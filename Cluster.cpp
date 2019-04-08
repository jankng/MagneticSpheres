//
// Created by jan on 4/5/19.
//
#include <iostream>
#include <cmath>
#include "Cluster.h"
#include "misc.h"

Cluster::Cluster(int n) {
    this->n = n;
    this->diam = 1;

    config = std::make_shared< std::vector<Dipole> >();
    config->reserve(n);
    for(int i = 0; i<n; i++){
        Dipole tmp;
        this->config->emplace_back(tmp);
    }

}

void Cluster::Print() {
    std::cout << "# Rx Ry Rz Mx My Mz" << std::endl;
    for(int i = 0; i<this->config->size(); i++){
        std::cout << i << " ";
        (*config)[i].Print();

    }

}

// TODO dipole moment has to be UNIT SIZE!!!
double Cluster::ComputeEnergy() {
    double ret = 0;

    for(int i = 0; i<this->n; i++){
        for(int j = 0; j<i; j++){
            std::vector<double> mi = (*config)[i].GetM();
            std::vector<double> mj = (*config)[i].GetM();

            std::vector<double> ri = (*config)[i].GetR();
            std::vector<double> rj = (*config)[i].GetR();

            std::vector<double> rij = (*config)[i].VectorTo((*config)[i]);
            double r = (*config)[i].DistanceTo((*config)[i]);

            ret += (misc::DotProduct(mi, mj) - 3*misc::DotProduct(mi, rij)*misc::DotProduct(mj, rij) / pow(r, 2))
                    / (pow(r / this->diam, 3) * this->n);
        }
    }

    this->energy = ret;

    return ret;
}

Dipole* Cluster::GetDipole(int id) {
    return &((*config)[id]);
}

double Cluster::dummy() {
    double ret = 0;
    for(int i = 0; i<this->n; i++){
        std::vector<double> tmp = (*config)[i].GetR();
        ret += tmp[2];
    }

    return ret;
}