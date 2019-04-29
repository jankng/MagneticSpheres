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
#include "misc.h"
#include "gsledits.h"
#include "conjgrad.h"

void minimize(int n){
    double e_min = 0;

    // TODO CHANGE 2 TO 250
    for(int i = 0; i<250; i++){
        std::cout << "Thread " << n << ", i=" << i << std::endl;
        metropolis candidate(8);
        candidate.start_siman();

        double e_candidate = candidate.get_cluster()->compute_energy();
        if(e_candidate < e_min){
            e_min = e_candidate;
            std::string filename = "t" + std::to_string(n) + "i" + std::to_string(i) +
                    "e" + std::to_string(e_min) + ".txt";
            candidate.get_cluster()->write_to_file(filename);
        }
    }

    std::cout << "Thread " << n << " ended." << std::endl;
}

double mod1(double x){
    x = fmod(x, 0.99);
    if(x < 0)
        x = x + 1;
    return x;
}

double
my_f (const gsl_vector *v, void *params) {
    std::vector<dipole> config;
    config.reserve(v->size / 5);

    for(int i = 0; i<v->size; i+=5){
        double phi = 2*M_PI*mod1(gsl_vector_get(v, i+3));
        double theta = M_PI*mod1(gsl_vector_get(v, i+4));

        config.emplace_back(dipole(gsl_vector_get(v, i+0),
                                   gsl_vector_get(v, i+1),
                                   gsl_vector_get(v, i+2),
                                   phi,
                                   theta));
    }

    cluster cl(config);
    return cl.compute_energy();
}

void
my_df (const gsl_vector *v, void *params,
       gsl_vector *df)
{
    std::vector<dipole> config;
    config.reserve(v->size / 5);

    for(int i = 0; i<v->size; i+=5){
        double phi = 2*M_PI*mod1(gsl_vector_get(v, i+3));
        double theta = M_PI*mod1(gsl_vector_get(v, i+4));

        config.emplace_back(dipole(gsl_vector_get(v, i+0),
                                   gsl_vector_get(v, i+1),
                                   gsl_vector_get(v, i+2),
                                   phi,
                                   theta));
    }

    cluster cl(config);
    //TODO consider heap allocation
    std::vector<double> grad;
    cl.compute_energy_gradient(&grad);

    for(int i = 0; i<grad.size(); i++){
        gsl_vector_set(df, i, grad[i]);
    }

    //debug
    //std::cout << gsl_vector_get(df, 20) << std::endl;
}

void
my_fdf (const gsl_vector *x, void *params,
        double *f, gsl_vector *df)
{
    *f = my_f(x, params);
    my_df(x, params, df);
}

void dosomething(){
    size_t iter = 0;
    int status;

    int spheres = 5;

    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;

    gsl_vector *x;
    gsl_multimin_function_fdf my_func;

    int trash = 0;

    my_func.n = 5*spheres;
    my_func.f = my_f;
    my_func.df = my_df;
    my_func.fdf = my_fdf;
    my_func.params = &trash;

    cluster cl(spheres, chain);
    cl.print();

    x = gsl_vector_alloc (5*spheres);
    for(int i = 0; i<spheres; i++){
        gsl_vector_set(x, 5*i+0, cl.get_dipole_by_ref(i)->get_r().at(0));
        gsl_vector_set(x, 5*i+1, cl.get_dipole_by_ref(i)->get_r().at(1));
        gsl_vector_set(x, 5*i+2, cl.get_dipole_by_ref(i)->get_r().at(2));


        gsl_vector_set(x, 5*i+3, cl.get_dipole_by_ref(i)->get_angles().at(0) / (2*M_PI));
        gsl_vector_set(x, 5*i+4, cl.get_dipole_by_ref(i)->get_angles().at(1) / M_PI);
    }

    T = gsl_multimin_fdfminimizer_conjugate_fr;
    s = gsl_multimin_fdfminimizer_alloc (T, 5*spheres);

    gsl_multimin_fdfminimizer_set (s, &my_func, x, 1, 1e-9);

    do
    {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate (s);

        if (status){
            std::cout << "error code" << status << std::endl;
            break;
        }

        status = gsl_multimin_test_gradient (s->gradient, 1e-3);

        if (status == GSL_SUCCESS)
            printf ("Minimum found at:\n");

        printf ("%5d %.5f %.5f %10.5f\n", (int) iter,
                gsl_vector_get (s->x, 0),
                gsl_vector_get (s->x, 1),
                s->f);

    }
    while (status == GSL_CONTINUE && iter < 100);


    for(int i = 0; i<s->x->size; i++)
        std::cout << gsl_vector_get(s->x, i) << std::endl;

    gsl_multimin_fdfminimizer_free (s);

    gsl_vector_free (x);
}

int main() {
    misc::setup_static_rng();


/*
    std::cout << "stepsize5" < std::endl;
    std::cout << "starting threads..." << std::endl;
    std::thread one(minimize, 1);
    std::thread two(minimize, 2);
    std::thread three(minimize, 3);
    //std::thread four(minimize, 4);

    one.join();
    two.join();
    three.join();
    //four.join();

    std::cout << "threads joined." << std::endl;

    gsl_rng* test = misc::get_static_rng();
    std::cout << misc::random_simple() << std::endl;
    std::cout << gsl_rng_uniform(test) << std::endl;
    */

    int N = 100;
    std::vector<dipole> conf;
    for(int i = 0; i<N; i++){
        dipole dp(0, 0, i, 0, 0);
        conf.emplace_back(dp);
    }


    std::vector<double> buffer;
    cluster cl(conf);
    LOG(cl.compute_energy());

    cl.compute_energy_gradient(&buffer);
    conjgrad cg(&cl);
    cg.dosomehting();

    misc::delete_static_rng();
    return 0;
}