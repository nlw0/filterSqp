//
// Created by nlw on 25/09/15.
//

#ifndef FILTERSQP_LM_OPTIM_H
#define FILTERSQP_LM_OPTIM_H

#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix_double.h>

double *levenberg_marquadt(void (*function)(double *, double *, double *, double *), int dim, double xx0[],
                           double step_size=0.001, double step_limit=1e-10, int max_iterations=10000);

gsl_vector *rec_levenberg_marquadt(void (*function)(double *, double *, double *, double *),
                                   gsl_vector *xx,
                                   gsl_vector *grad,
                                   gsl_matrix *hess,
                                   gsl_matrix *Hu,
                                   gsl_matrix *Hv,
                                   gsl_vector *Hs,
                                   gsl_vector *Hwrk,
                                   gsl_vector *xxStep,
                                   int iteration,
                                   double step_size, double step_limit, int max_iterations);

#endif //FILTERSQP_LM_OPTIM_H
