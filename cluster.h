//
// Created by jan on 4/5/19.
//

#ifndef MAGNETICSPHERES_CLUSTER_H
#define MAGNETICSPHERES_CLUSTER_H


#include <vector>
#include "dipole.h"
#include "misc.h"

#define CLUSTER_DEFAULT_DIAMETER 1
#define CLUSTER_DEFAULT_CHAIN_HEIGHT 5

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

    // helper functions for gradient
    double gradient_dx(int i, int x);
    double gradient_dtheta(int i);
    double gradient_dphi(int i);

public:
    cluster(){}
    // create cluster from given config
    explicit cluster(const std::vector<dipole> &config);
    explicit cluster(const std::vector<double> &config);

    // generate random cluster of cluster_size n
    explicit cluster(int n);

    //generate certain shape with random m
    cluster(int n, shape s);

    //copy constructor
    cluster(const cluster& ori);

    //getters
    void config_to_vec(std::vector<double>* target);
    int get_size(){return cluster_size;}
    shape get_shape(){return cluster_shape;}
    dipole* get_dipole_by_ref(int id);

    //output
    void print();
    std::string to_string(char sep = ' ');
    void write_to_file(const std::string& filename = misc::get_time() + ".txt");

    // computes energy and sets private 'energy' variable
    double compute_energy();
    double compute_energy_for_metropolis();
    double compute_energy_for_gradient();
    void compute_energy_gradient(std::vector<double>* ret, int index);
    void compute_angle_gradient(std::vector<double>* ret, int index);
    void compute_coordinate_gradient(std::vector<double>* ret, int index);

    // other
    bool is_valid();
};


#endif //MAGNETICSPHERES_CLUSTER_H
