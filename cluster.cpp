//
// Created by jan on 4/5/19.
//
#include <iostream>
#include <cmath>
#include <limits>
#include "cluster.h"
#include "misc.h"

cluster::cluster(int n) {
    cluster_shape = other;
    this->cluster_size = n;
    this->diameter = CLUSTER_DEFAULT_DIAMETER;
    make_other();
}

cluster::cluster(int n, shape s){
    cluster_shape = s;
    this->cluster_size = n;
    this->diameter = CLUSTER_DEFAULT_DIAMETER;
    // TODO implement more shapes
    if(s == chain)
        make_chain();
    else if (s == cube)
        make_cube();
    else if (s == plane)
        make_plane();
}

void cluster::make_chain(){
    config.reserve(cluster_size);
    config = {};
    for(int i = 0; i<cluster_size; i++){
        dipole d(0, 0, i);
        config.emplace_back(d);
    }
}

// TODO prevent spheres form crashing
void cluster::make_other(){
    this->config.reserve(cluster_size);
    for(int i = 0; i<cluster_size; i++){

        // prevent crashing
        bool is_accepted = false;
        while(!is_accepted){
            bool flag = true;
            dipole d;
            for(int j = 0; j<i; j++) {
                double dist = d.distance_to(config[j]);
                if (dist < 1) {
                    flag = false;
                    break;
                }
            }

            if(flag){
                is_accepted = true;
                config.emplace_back(d);
            }
        }
    }
}

void cluster::print() {
    std::cout << "# Rx Ry Rz Mx My Mz" << std::endl;
    for(int i = 0; i<this->config.size(); i++){
        std::cout << i << " ";
        this->config[i].print();
    }

}

double cluster::compute_energy() {
    double ret = 0;

    for(int i = 0; i<cluster_size; i++){
        for(int j = 0; j<i; j++){
            std::vector<double> mi = config[i].get_m();
            std::vector<double> mj = config[j].get_m();

            std::vector<double> ri = config[i].get_r();
            std::vector<double> rj = config[j].get_r();

            std::vector<double> rij = config[i].vector_to(config[j]);
            double r = config[i].distance_to(config[j]);

            //TODO define out of bounds properly
            if(r < diameter || r > 1000) {
                //std::cout << "Spheres crashed into each other" << std::endl;
                return std::numeric_limits<double>::max();
            }

            ret += (misc::dot_product(mi, mj) - 3* misc::dot_product(mi, rij)* misc::dot_product(mj, rij) / pow(r, 2))
                    / pow(r, 3);
        }
    }

    // multiply by 1/(U_up_up * n)
    ret *= pow(diameter, 3) / (cluster_size * pow(DIPOLE_DEFAULT_M, 2));
    return ret;
}

dipole* cluster::get_dipole_by_ref(int id) {
    return &(this->config[id]);
}

cluster::cluster(const std::vector<dipole> &config) {
    cluster_shape = other;
    this->cluster_size = config.size();
    this->diameter = CLUSTER_DEFAULT_DIAMETER;
    this->config = config;

}

void cluster::make_cube() {
    config.reserve(cluster_size);
    config = {};
    for(int i = 0; i<cluster_size; i++){
        for(int j = 0; j<cluster_size; j++){
            for(int k = 0; k<cluster_size; k++){
                dipole d(i, j, k);
                config.emplace_back(d);
            }
        }
    }

    // not using pow() because cluster_size is of type int.
    cluster_size = cluster_size * cluster_size * cluster_size;

}

// TODO implement correctly as of rn make_plane == make_other
void cluster::make_plane() {
    config.reserve(cluster_size);
    config = {};
    for(int i = 0; i<cluster_size; i++){
        for(int j = 0; j<cluster_size; j++){
            dipole d;
            config.emplace_back(d);
        }
    }
}

cluster::cluster(const cluster &ori) {
    config = ori.config;
    cluster_shape = ori.cluster_shape;
    cluster_size = ori.cluster_size;
    diameter = ori.diameter;

}
