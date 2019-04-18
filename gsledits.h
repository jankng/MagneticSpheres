//
// Created by jan on 4/14/19.
//

#ifndef MAGNETICSPHERES_GSLEDITS_H
#define MAGNETICSPHERES_GSLEDITS_H

#include <gsl/gsl_machine.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_siman.h>

namespace gsl_edits{
    void
    gsl_siman_solve (const gsl_rng * r, void *x0_p, gsl_siman_Efunc_t Ef,
                     gsl_siman_step_t take_step,
                     gsl_siman_print_t print_position,
                     gsl_siman_copy_t copyfunc,
                     gsl_siman_copy_construct_t copy_constructor,
                     gsl_siman_destroy_t destructor,
                     size_t element_size,
                     gsl_siman_params_t params);
}

#endif //MAGNETICSPHERES_GSLEDITS_H
