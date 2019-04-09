//
// Created by jan on 4/5/19.
//
#include <iostream>
#include <cmath>
#include "cluster.h"
#include "misc.h"

cluster::cluster(int n) {
    cluster_shape = other;
    this->cluster_size = n;
    this->diameter = DEFAULT_DIAMETER;
    make_other();
}

cluster::cluster(int n, shape s){
    cluster_shape = s;
    this->cluster_size = n;
    this->diameter = DEFAULT_DIAMETER;
    // TODO implement more shapes
    if(s == chain)
        make_chain();
}

void cluster::make_chain(){
    config.reserve(cluster_size);
    config = {};
    for(int i = 0; i<cluster_size; i++){
        dipole d(0, 0, i);
        config.emplace_back(d);
    }
}

void cluster::make_other(){
    this->config.reserve(cluster_size);
    for(int i = 0; i<cluster_size; i++){
        dipole tmp;
        this->config.emplace_back(tmp);
    }
}

void cluster::print() {
    std::cout << "# Rx Ry Rz Mx My Mz" << std::endl;
    for(int i = 0; i<this->config.size(); i++){
        std::cout << i << " ";
        this->config[i].print();
    }

}

// TODO dipole moment has to be UNIT SIZE!!!
double cluster::compute_energy() {
    double ret = 0;

    for(int i = 0; i<cluster_size; i++){
        for(int j = 0; j<i; j++){
            std::vector<double> mi = this->config[i].get_m();
            std::vector<double> mj = this->config[j].get_m();

            std::vector<double> ri = this->config[i].get_r();
            std::vector<double> rj = this->config[j].get_r();

            std::vector<double> rij = this->config[i].vector_to(this->config[j]);
            double r = this->config[i].distance_to(this->config[j]);

            ret += (misc::dot_product(mi, mj) - 3* misc::dot_product(mi, rij)* misc::dot_product(mj, rij) / pow(r, 2))
                    / (pow(r / this->diameter, 3) * this->cluster_size);
        }
    }

    return ret;
}

dipole* cluster::get_dipole_by_ref(int id) {
    return &(this->config[id]);
}

cluster::cluster(std::vector<dipole> config) {
    cluster_shape = other;
    this->cluster_size = config.size();
    this->diameter = DEFAULT_DIAMETER;
    this->config = config;

}
