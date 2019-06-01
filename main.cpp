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

/* Paraboloid centered on (p[0],p[1]), with
   scale factors (p[2],p[3]) and minimum p[4] */

double
my_f (const gsl_vector *v, void *params)
{
    std::vector<double> coords;
    coords.reserve(v->size);
    for(int i = 0; i<v->size; i++){
        coords.emplace_back(gsl_vector_get(v, i));
    }

    cluster cl(coords);
    return cl.compute_energy();
}

/* The gradient of f, df = (df/dx, df/dy). */
void
my_df (const gsl_vector *v, void *params,
       gsl_vector *df)
{
    int *p = (int *)params;

    std::vector<double> coords;
    coords.reserve(v->size);
    for(int i = 0; i<v->size; i++){
        coords.emplace_back(gsl_vector_get(v, i));
    }

    cluster cl(coords);
    std::vector<double> grad;
    cl.compute_energy_gradient(&grad, p[0]);

    for(int i = 0; i<df->size; i++){
        gsl_vector_set(df, i, grad[i]);
    }
}

/* Compute both f and df together. */
void
my_fdf (const gsl_vector *x, void *params,
        double *f, gsl_vector *df)
{
    *f = my_f(x, params);
    my_df(x, params, df);
}

void conjmin(cluster* init){
    int N = 8;
    size_t iter = 0;
    int status;

    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;

    /* Position of the minimum (1,2), scale factors
       10,20, height 30. */
    int par[1] = { 2 };

    gsl_vector *x;
    gsl_multimin_function_fdf my_func;

    my_func.n = N*5;
    my_func.f = my_f;
    my_func.df = my_df;
    my_func.fdf = my_fdf;
    my_func.params = par;

    /* Starting point, x = (5,7) */
    x = gsl_vector_alloc (N*5);
    std::vector<double> conf;
    init->config_to_vec(&conf);
    for(int i = 0; i<5*N; i++){
        gsl_vector_set(x, i, conf[i]);
    }

    T = gsl_multimin_fdfminimizer_conjugate_fr;
    s = gsl_multimin_fdfminimizer_alloc (T, 5*N);

    gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-4);

    do
    {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate (s);

        if (status) {
            std::cout << "status: " << status << std::endl;
            break;
        }

        status = gsl_multimin_test_gradient (s->gradient, 1e-3);

        if (status == GSL_SUCCESS)
            printf ("Minimum found at:\n");

        printf ("%5d %.5f %.5f %10.5f\n", iter,
                gsl_vector_get (s->x, 0),
                gsl_vector_get (s->x, 1),
                s->f);

    }
    while (status == GSL_CONTINUE && iter < 100);

    gsl_multimin_fdfminimizer_free (s);
    gsl_vector_free (x);
}

cluster* make_perfect_chain(int n){
    std::vector<dipole> dps;
    for(int i = 0; i<n; i++){
        dipole d(0, 0, i, 0, 0);
        dps.emplace_back(d);
    }

    cluster* c = new cluster(dps);
    return c;
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

    //conjgrad c(8);
    //c.print_energy_in_direction(nullptr);
    //c.minimize_simultaneous();
    //c.dosomething();

    cluster* perfect_chain = new cluster(8);
    perfect_chain->print();

    metropolis m(perfect_chain);
    m.enable_verbose_mode();
    m.start_siman();

    perfect_chain->print();


    delete perfect_chain;

    misc::delete_static_rng();
    return 0;
}