// C++ Standard Library
#include <iostream>
#include <cmath>
#include <cstring>
#include <thread>
#include <string>
#include <iomanip>

// GSL Headers
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>

// Project Headers
#include "definitions.h"
#include "dipole.h"
#include "cluster.h"
#include "metropolis.h"
#include "acmetropolis.h"
#include "misc.h"
#include "gsledits.h"
#include "conjgrad.h"

void startMetropolis(int n){
    int iterations= 250;
    double e_min = 0;

    //determine parameters
    gsl_siman_params_t params_given;
    switch(n) {
        case 1:
            params_given = {0, 1000, 0.1, 1.0, 10, 1.001, 0.0005};
            break;
        case 2:
            params_given = {0, 1000, 0.1, 1.0, 0.001, 1.001, 0.0005};
            break;
        case 3:
            params_given = {0, 1000, 0.1, 1.0, 10, 1.001, 0.00000005};
            break;
        default:
            params_given = {0, ITERS_FIXED_T, STEP_SIZE, 1.0, INITIAL_T, MU_T, T_MIN};
    }

    //start iterations
    for(int i = 0; i<iterations; i++){
        std::cout << "Thread " << n << ", i=" << i << std::endl;
        auto cl = new cluster(8, chain);
        metropolis m(cl, params_given);
        m.start_siman();

        double e_candidate = m.get_cluster()->compute_energy_for_metropolis();
        if(e_candidate < e_min){
            e_min = e_candidate;
            std::string filename = "t" + std::to_string(n) + "i" + std::to_string(i) +
                    "e" + std::to_string(e_min) + ".txt";
            m.get_cluster()->write_to_file(filename);
        }

        delete cl;
    }

    std::cout << "Thread " << n << " ended." << std::endl;
}

void startMetropolisThreads(){
    std::thread t1(startMetropolis, 1);
    std::thread t2(startMetropolis, 2);
    std::thread t3(startMetropolis, 3);

    t1.join();
    t2.join();
    t3.join();

    std::cout << "Threads joined successfully";
}

cluster* make_perfect_chain(int n){
    std::vector<dipole> dps;
    for(int i = 0; i<n; i++){
        dipole d(i, 0, 5, 0, M_PI / 2.0);
        dps.emplace_back(d);
    }

    cluster* c = new cluster(dps);
    return c;
}

cluster* make_random_chain(int n){
    std::vector<dipole> dps;
    for(int i = 0; i<n; i++){
        double angle = misc::random_simple() * 2.0 * M_PI;

        dipole d(i, 0, 5, cos(angle), 0, sin(angle));
        dps.emplace_back(d);
    }

    cluster* c = new cluster(dps);
    return c;
}


double someFunc(double x){
    return (x-5)*(x-4)*(x-3)*(x-2);
}

void shft3(double& a, double& b, double& c, const double& d){
    a = b;
    b = c;
    c = d;
}

void SWAP(double& x, double& y){
    double temp = x;
    x = y;
    y = temp;
}

double SIGN(double x, double y){
    double sgn;
    if (x == 0)
        sgn = 0;
    else
        sgn = (y < 0) ? -1.0 : 1.0;

    return abs(x) * sgn;
}

void mnbrak(double& ax, double& bx, double& cx, double& fa, double& fb, double& fc, double func(double)){
    const double GOLD=1.618034,GLIMIT=100.0,TINY=1.0e-20;
    double ulim,u,r,q,fu;

    fa=func(ax);
    fb=func(bx);
    if (fb > fa) {
        SWAP(ax,bx);
        SWAP(fb,fa);
    }
    cx=bx+GOLD*(bx-ax);
    fc=func(cx);
    while (fb > fc) {
        r=(bx-ax)*(fb-fc);
        q=(bx-cx)*(fb-fa);
        u=bx-((bx-cx)*q-(bx-ax)*r)/
             (2.0*SIGN(fmax(fabs(q-r),TINY),q-r));
        ulim=bx+GLIMIT*(cx-bx);
        if ((bx-u)*(u-cx) > 0.0) {
            fu=func(u);
            if (fu < fc) {
                ax=bx;
                bx=u;
                fa=fb;
                fb=fu;
                return;
            } else if (fu > fb) {
                cx=u;
                fc=fu;
                return;
            }
            u=cx+GOLD*(cx-bx);
            fu=func(u);
        } else if ((cx-u)*(u-ulim) > 0.0) {
            fu=func(u);
            if (fu < fc) {
                shft3(bx,cx,u,cx+GOLD*(cx-bx));
                shft3(fb,fc,fu,func(u));
            }
        } else if ((u-ulim)*(ulim-cx) >= 0.0) {
            u=ulim;
            fu=func(u);
        } else {
            u=cx+GOLD*(cx-bx);
            fu=func(u);
        }
        shft3(ax,bx,cx,u);
        shft3(fa,fb,fc,fu);
    }

}

void test(){
    double min = 100;

    for(int i = 0; i<1000000; i++){
        cluster* candidate = make_random_chain(8);
        conjgrad cand(candidate);
        cand.minimize_simultaneous();

        cluster cl(cand.getConfig());
        double nrg = cl.compute_energy_for_gradient();
        if (nrg < min){
            min = nrg;
            std::cout << "New optimal cluster. E, is valid?: " << nrg << "  " << cl.is_valid() << std::endl;
            cl.print();
        }

        delete candidate;
    }
}

int main() {
    misc::setup_static_rng();

    LOG("Hello World!");

    std::vector<dipole> conf = {
        dipole(2.0889, 2.98698, 2.91663, 0.346189, 0.716947, 0.605095),
        dipole(2.60734, 3.7125, 3.37009, 0.609099, 0.763548, 0.214461),
        dipole(3.3155, 4.41318, 3.27233, 0.807048, 0.334958, -0.486288),
        dipole(2.608, 1.89024, 1.56034, -0.943785, 0.0137649, 0.330274),
        dipole(3.47682, 2.38705, 1.55014, -0.586287, -0.786966, -0.192231),
        dipole(3.84866, 3.23522, 1.9278, -0.177555, -0.847542, -0.500148),
        dipole(1.96104, 2.29123, 2.20942, -0.215821, 0.566178, 0.795527),
        dipole(3.85499, 4.04629, 2.51362, 0.175769, -0.672323, -0.719088)
    };

    cluster cl(conf);

    cluster* fixed = make_perfect_chain(8);
    std::vector<double> grad;
    fixed->compute_energy_gradient(&grad, -1);

    for(int i = 0; i<grad.size(); i+=5){
        std::cout << grad[i] << " " << grad[i+1] << " " << grad[i+2] << " " << grad[i+3] << " " << grad[i+4] << std::endl;
    }

    auto* rnd = make_random_chain(8);
    //rnd->print();
    //conjgrad c(fixed);
    //c.print_energy_in_direction(nullptr);
    //c.minimize_simultaneous();
    //c.minimize_single_dipoles();
    //c.minimize_simultaneous();

    //conjgrad c2(chain);
    //c.minimize_simultaneous();


    acmetropolis m(8);
    m.enable_verbose_mode();
    m.start_siman();
    cluster res;
    acmetropolis::ac_to_cluster(*(m.get_cluster()), &res);
    res.print();
    std::cout << res.compute_energy() << std::endl;

    //startMetropolisThreads();

    //test();

    misc::delete_static_rng();
    return 0;
}