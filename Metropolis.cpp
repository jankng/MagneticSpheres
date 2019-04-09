//
// Created by jkoenig on 09/04/19.
//
/*
#include "Metropolis.h"
#include "misc.h"
#include <cmath>

Metropolis::Metropolis() {
    this->cluster = std::make_shared<cluster>(8);
    energy_current = cluster->compute_energy();
    dipole_count = 8;
    beta = 1;

}

void Metropolis::minimization_attempt_random() {
    std::shared_ptr<cluster> candidate = std::make_shared<cluster>(this->dipole_count);

    double energy_candidate = candidate->ComputeEnergy();
    if(energy_candidate < energy_current) {
        cluster = candidate;
        energy_current = energy_candidate;
        step_count++;
    }
    else if(misc::random_simple() < exp(-1*beta * (energy_current - energy_candidate) )) {
        cluster = candidate;
        energy_current = energy_candidate;
        step_count++;
    }

    attempt_count++;

}

void Metropolis::minimize_random(int n) {
    for(int i = 0; i<n; i++){
        minimization_attempt_random();
    }

}

void Metropolis::test() {
    for(int i = 0; i<dipole_count; i++){
        dipole* d = cluster->get_dipole_by_ref(i);
        d->SetR(0, 0, i);
    }

    minimize_random(100);

}

void Metropolis::minimization_attempt_random_angle() {
    std::shared_ptr<cluster> candidate = std::make_shared<cluster>(this->dipole_count);

    double energy_candidate = candidate->ComputeEnergy();
    if(energy_candidate < energy_current) {
        cluster = candidate;
        energy_current = energy_candidate;
        step_count++;
    }
    else if(misc::random_simple() < exp(-1*beta * (energy_current - energy_candidate) )) {
        cluster = candidate;
        energy_current = energy_candidate;
        step_count++;
    }

    attempt_count++;

}

*/
