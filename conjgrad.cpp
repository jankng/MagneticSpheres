//
// Created by jkoenig on 29/04/19.
//

#include "conjgrad.h"
#include <iostream>
#include <cmath>

conjgrad::conjgrad(std::vector<double>* conf) : config(*conf) {

}

conjgrad::conjgrad(int n) {
    cluster cl(n, chain);
    cl.config_to_vec(&config);
    compute_gradient(-1);
}

void conjgrad::compute_gradient(int i) {
    cluster cl(config);
    cl.compute_energy_gradient(&grad, i);
}

double conjgrad::compute_energy() {
    cluster cl(config);
    return cl.compute_energy();
}

double conjgrad::compute_energy_in_direction(double t, const std::vector<double>& dir) {
    std::vector<double> temp;
    temp.reserve(config.size());

    for(int i = 0; i<config.size(); i++){
        temp.emplace_back(config[i] - t * dir[i]);
    }

    cluster clst(temp);
    return clst.compute_energy();
}

void conjgrad::print_energy_in_direction(std::vector<double>* dir) {
    LOG("t \t E(t)");
    for(int t = 0; t < 10; t++)
        std::cout << (double) t / 10 << "\t" << compute_energy_in_direction((double) t / 10, *dir) << std::endl;
    //double trash = minimize_in_direction(grad);
}

double conjgrad::minimize_in_direction(const std::vector<double>& dir) {
    int steps = 0;
    static double w = (3.0 - sqrt(5.0)) / 2.0;
    double x = 0;
    double a = 0;
    double b = 0.5;
    double c = 1;

    while(fmax(c - b, b - a) > 0.01 && steps < 1000){
        if(b - a >= c - b){
            x = a + (1.0 - w)*(b - a);
            if(compute_energy_in_direction(b, dir) < compute_energy_in_direction(x, dir)){
                a = x;
            } else{
                c = b;
                b = x;
            }
        } else{
            x = b + w * (c - b);
            if(compute_energy_in_direction(b, dir) < compute_energy_in_direction(x, dir)){
                c = x;
            } else{
                a = b;
                b = x;
            }
        }

        steps++;
    }

    std::cout << "steps, t, nrg\t" << steps << " " << x << " " << compute_energy_in_direction(x, grad) << std::endl;
    return x;
}

void conjgrad::go_in_direction(double t, const std::vector<double>& dir) {
    for(int i = 0; i < config.size(); i++){
        config[i] = config[i] - t * dir[i];
    }

}

void conjgrad::minimize_simultaneous() {
    int j = 1;
    std::vector<double> r  = grad;

    while(j < 1000 && misc::dot_product(grad, grad) > 0.001){
        print_energy_in_direction(&r);
        double t = minimize_in_direction(r);

        go_in_direction(t, r);
        std::vector<double> grad_old = grad;
        compute_gradient(-1);

        double gamma = misc::dot_product(grad, grad) / misc::dot_product(grad_old, grad_old);
        for(int i = 0; i<r.size(); i++){
            r[i] = grad[i] + gamma * r[i];
        }


        j++;
    }

    std::cout << "Simultaneous minimization done after ... steps: " << j << std::endl;

}

void conjgrad::dosomething() {

    print_energy_in_direction(&grad);
    double t = minimize_in_direction(grad);
    go_in_direction(t, grad);
    compute_gradient(-1);
    print_energy_in_direction(&grad);
    t = minimize_in_direction(grad);

}
