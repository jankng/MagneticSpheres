//
// Created by jan on 4/5/19.
//
#include <iostream>
#include <cmath>
#include <limits>
#include <sstream>
#include <fstream>
#include "cluster.h"
#include "misc.h"

#define RESULTS_DIR "/bigwork/jkoenig/results/"
//#define RESULTS_DIR "/home/jan/Desktop/"

#define PENALTY 3
#define GRAVITATION 0.6

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
    double penalty = PENALTY;
    double ret = 0;

    for(int i = 0; i<cluster_size; i++){
        std::vector<double> ri = config[i].get_r();
        std::vector<double> mi = config[i].get_m();

        for(int j = 0; j<i; j++){
            std::vector<double> mj = config[j].get_m();
            std::vector<double> rj = config[j].get_r();

            std::vector<double> rij = config[i].vector_to(config[j]);
            double r = config[i].distance_to(config[j]);

            if(!(config[i].is_in_bounds()) || !(config[j].is_in_bounds())) {
                //std::cout << "Spheres crashed into each other" << std::endl;
                //return std::numeric_limits<double>::max();
            }

            if(r - diameter < 0)
                ret -= penalty*(r-diameter);

            // TODO check diameter implementation
            ret += (misc::dot_product(mi, mj) - 3* misc::dot_product(mi, rij)* misc::dot_product(mj, rij) / pow(r, 2))
                    / pow(r, 3);
        }

        ret += GRAVITATION*ri[2];

        for(int i = 0; i<3; i++){
            if(ri[i] > 50)
                ret += penalty*(ri[i] - 50);
            if(ri[i] < 0)
                ret -= penalty*ri[i];
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

std::string cluster::to_string(char sep) {
    std::stringstream ret;

    for(int i = 0; i<cluster_size; i++){
        ret << config[i].to_string(sep) << "\n";
    }

    return ret.str();
}

void cluster::write_to_file(const std::string& filename) {
    std::string path = RESULTS_DIR + filename;
    std::ofstream handler(path, std::ios_base::app);
    if(handler.is_open()){
        handler << to_string();
        handler.close();
    }
}

void cluster::compute_energy_gradient(std::vector<double>* ret, int index) {
    int components = cluster_size *5;
    ret->reserve(components);

    for(int i = 0; i<cluster_size; i++){
        if(i == 0 || i == cluster_size - 1){
            for(int j = 0; j<5; j++){
                ret->emplace_back(0.0);
            }
        }
        else if(i == index || index < 0) {
            ret->emplace_back(gradient_dx(i, 0));
            ret->emplace_back(gradient_dx(i, 1));
            ret->emplace_back(gradient_dx(i, 2));
            ret->emplace_back(gradient_dphi(i));
            ret->emplace_back(gradient_dtheta(i));
            //ret->emplace_back(0.0);
            //ret->emplace_back(0.0);
        } else{
            for(int j = 0; j<5; j++){
                ret->emplace_back(0.0);
            }
        }
    }
}

double cluster::gradient_dx(int i, int x) {
    double penalty = PENALTY;
    double res = 0;
    std::vector<double> ri = config[i].get_r();
    std::vector<double> mi = config[i].get_m();
    double rix = ri[x];
    double mix = mi[x];

    for(int j = 0; j<cluster_size; j++){
        if(i == j)
            continue;

        double r = config[i].distance_to(config[j]);
        std::vector<double> rj = config[j].get_r();


        std::vector<double> rij = config[i].vector_to(config[j]);
        std::vector<double> mj = config[i].get_m();
        double mjx = mj[x];

        double nom = misc::dot_product(mi, mj)*pow(r, 2) - 3*misc::dot_product(mi, rij)* misc::dot_product(mj, rij);
        double denom = pow(r, -5);

        double d_nom = 2.0*(rj[x] - ri[x])*misc::dot_product(mi, mj) - 3.0*(mix*misc::dot_product(mj, rij) + mjx*misc::dot_product(mj, rij));
        double d_denom = -5.0*(rj[x] - ri[x])*pow(r, -7);

        res += (nom*d_denom + d_nom*denom) / (2*cluster_size*pow(diameter, 3));

        if(r < diameter)
            res += -1.0*penalty * (rj[x] - ri[x]) / r;
    }

    //res = res / (cluster_size*pow(diameter, 3));

    if(ri[x] > 50)
        res += penalty;
    if(ri[x] < 0)
        res -= penalty;

    //gravity
    if(x == 2)
        res += GRAVITATION;

    return res;
}

double cluster::gradient_dphi(int i) {
    double ret = 0;
    std::vector<double> ri = config[i].get_r();
    std::vector<double> angsi = config[i].get_angles();
    std::vector<double> dmi = {-1.0*sin(angsi[1])*sin(angsi[0]),
                               sin(angsi[1]) * cos(angsi[0]),
                               0};

    for(int j = 0; j<cluster_size; j++){
        if(i == j)
            continue;

        std::vector<double> angsj = config[j].get_angles();
        std::vector<double> dmj = {-1.0*sin(angsj[1])*sin(angsj[0]),
                                   sin(angsj[1]) * cos(angsj[0]),
                                   0};


        std::vector<double> mj = config[j].get_m();

        std::vector<double> rj = config[j].get_r();

        std::vector<double> rij = config[i].vector_to(config[j]);
        double r = config[i].distance_to(config[j]);

        ret += (misc::dot_product(dmi, mj) - 3* misc::dot_product(dmi, rij)* misc::dot_product(mj, rij) / pow(r, 2))
               / pow(r, 3);
    }

    return ret;
}

double cluster::gradient_dtheta(int i){
    double ret = 0;
    std::vector<double> ri = config[i].get_r();
    std::vector<double> angsi = config[i].get_angles();
    std::vector<double> dmi = {cos(angsi[1])*cos(angsi[0]),
                               cos(angsi[1])*sin(angsi[0]),
                               -1.0*sin(angsi[1])};

    for(int j = 0; j<cluster_size; j++){
        if(i == j)
            continue;

        std::vector<double> angsj = config[j].get_angles();
        std::vector<double> dmj = {cos(angsj[1])*cos(angsj[0]),
                                   cos(angsj[1])*sin(angsj[0]),
                                   -1.0*sin(angsj[1])};


        std::vector<double> mj = config[j].get_m();

        std::vector<double> rj = config[j].get_r();

        std::vector<double> rij = config[i].vector_to(config[j]);
        double r = config[i].distance_to(config[j]);

        ret += (misc::dot_product(dmi, mj) - 3* misc::dot_product(dmi, rij)* misc::dot_product(mj, rij) / pow(r, 2))
               / pow(r, 3);
    }

    return ret;
}

void cluster::config_to_vec(std::vector<double>* target) {
    int size = cluster_size*5;
    target->clear();
    target->reserve(size);
    for(int i = 0; i<cluster_size; i++){
        std::vector<double> r = config[i].get_r();
        std::vector<double> angs = config[i].get_angles();
        target->emplace_back(r[0]);
        target->emplace_back(r[1]);
        target->emplace_back(r[2]);
        target->emplace_back(angs[0]);
        target->emplace_back(angs[1]);
    }
}

cluster::cluster(const std::vector<double> &config) :
cluster_size(config.size() / 5), diameter(CLUSTER_DEFAULT_DIAMETER), cluster_shape(other){
    this->config.reserve(cluster_size);
    for(int i = 0; i<cluster_size; i++){
        std::vector<double> dp = {config[5*i+0],
                                  config[5*i+1],
                                  config[5*i+2],
                                  misc::modn(config[5*i+3], 2*M_PI),
                                  misc::modn(config[5*i+4], M_PI)};
        this->config.emplace_back(dipole(dp));
    }
}
