//
// Created by jkoenig on 29/04/19.
//

#include "conjgrad.h"
#include <iostream>
#include <cmath>

conjgrad::conjgrad(std::vector<double>* conf) : config(*conf) {
    compute_gradient(-1, 0);
}

conjgrad::conjgrad(int n) {
    cluster cl(n);
    cl.config_to_vec(&config);
    compute_gradient(-1, 0);
}

conjgrad::conjgrad(cluster* cl){
    cl->config_to_vec(&config);
    compute_gradient(-1, 0);
}

void conjgrad::compute_gradient(int i, int mode) {
    cluster cl(config);
    if(mode == 0)
        cl.compute_energy_gradient(&grad, i);
    else if(mode == 1)
        cl.compute_coordinate_gradient(&grad, i);
    else if(mode == 2)
        cl.compute_angle_gradient(&grad, i);
}

double conjgrad::compute_energy() {
    cluster cl(config);
    return cl.compute_energy_for_gradient();
}

double conjgrad::compute_energy_in_direction(double t, const std::vector<double>& dir) {
    std::vector<double> temp;
    temp.reserve(config.size());

    for(int i = 0; i<config.size(); i++){
        temp.emplace_back(config[i] - t * dir[i]);
    }

    cluster clst(temp);
    return clst.compute_energy_for_gradient();
}

void conjgrad::print_energy_in_direction(std::vector<double>* dir) {
    LOG("t \t E(t)");
    for(int t = -5; t < 6; t++)
        std::cout << (double) t / 100 << "\t" << compute_energy_in_direction((double) t / 100, *dir) << std::endl;
    //double trash = minimize_in_direction(grad);
}

double conjgrad::minimize_in_direction(const std::vector<double>& dir) {
    int steps = 0;
    static double w = (3.0 - sqrt(5.0)) / 2.0;
    double x = 0;
    double a = -0.1;
    double b = 0;
    double c = 0;

    double s = 0.01;
    double d = s;
    while(true) {
        if (compute_energy_in_direction(d, dir) >= compute_energy_in_direction(a, dir)) {
            c = d;
            b = c / 2;
            break;
        } else {
            d += s;
        }
    }


    while(fmax(c - b, b - a) > 0.0001 && steps < 1000){
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
    int j = 0;
    std::vector<double> r  = grad;
    double min = 9001;
    double nrg = 0;

    while(j < 1000 && misc::dot_product(r, r) > 0.1){
        //print_energy_in_direction(&r);
        double t = minimize_in_direction(r);

        go_in_direction(t, r);

        cluster cl(config);
        nrg = cl.compute_energy_for_gradient();
        std::cout << "nrg, valid? " << cl.compute_energy() << " " << cl.is_valid() << std::endl;
        if(nrg < min){
            min = nrg;
            cl.print();
        }

        std::vector<double> grad_old = grad;
        compute_gradient(-1, 0);

        double gamma = misc::dot_product(grad, grad) / misc::dot_product(grad_old, grad_old);
        for(int i = 0; i<r.size(); i++){
            r[i] = grad[i] + gamma * r[i];
        }


        j++;
    }

    std::cout << "Simultaneous minimization done after ... steps: " << j << std::endl;

}

void conjgrad::minimize_single_dipoles() {

    cluster cl(config);
    cl.print();

    double min = 200;
    for(int j = 0; j<100; j++) {
        for (int i = 1; i < 7; i++) {
            compute_gradient(i, 1);
            print_energy_in_direction(&grad);
            double t = minimize_in_direction(grad);
            go_in_direction(t, grad);


            compute_gradient(i, 2);
            //print_energy_in_direction(&grad);
            t = minimize_in_direction(grad);
            go_in_direction(t, grad);

            cl = cluster(config);
            std::cout << "energy, valid? " << cl.compute_energy() << " " << cl.is_valid() << std::endl;
            double nrg = cl.compute_energy_for_gradient();
            if(nrg < min){
                min = nrg;
                cl.print();
            }
        }
    }

    cl = cluster(config);
    cl.print();

    //print_energy_in_direction(&grad);


}
