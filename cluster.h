//
// Created by jan on 4/5/19.
//

#ifndef MAGNETICSPHERES_CLUSTER_H
#define MAGNETICSPHERES_CLUSTER_H


#include <vector>
#include "dipole.h"

#define DEFAULT_DIAMETER 1

typedef enum{
    cube,
    chain,
    ring,
    plane,
    other
} shape;

class cluster {
private:
    std::vector<dipole> config;
    shape cluster_shape;
    int cluster_size;
    double diameter;

    // makers for shapes
    void make_chain();
    void make_cube();
    void make_plane();
    void make_other();
public:
    // create cluster from given config
    explicit cluster(const std::vector<dipole> &config);

    // generate random cluster of cluster_size n
    explicit cluster(int n);

    //generate certain shape with random m
    cluster(int n, shape s);

    //getters
    std::vector<dipole> get_config(){return config;}
    double get_size(){return cluster_size;}
    double get_shape(){return cluster_shape;}
    dipole* get_dipole_by_ref(int id);

    void print();

    // computes energy and sets private 'energy' variable
    double compute_energy();
};


#endif //MAGNETICSPHERES_CLUSTER_H
