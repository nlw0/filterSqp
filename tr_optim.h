//
// Created by nlw on 25/09/15.
//

#ifndef FILTERSQP_LM_OPTIM_H
#define FILTERSQP_LM_OPTIM_H

#include <iostream>
#include <iomanip>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix_double.h>
#include "find_root.h"

double *trust_region_optimization(void (*function)(double *, double *, double *, double *), int dim, double *xx0,
                                  double step_size = 0.001, double step_limit = 1e-10, int max_iterations = 10000);

gsl_vector *rec_trust_region_optimization(void (*function)(double *, double *, double *, double *),
                                          gsl_vector *xx,
                                          gsl_vector *yy,
                                          gsl_vector *grad,
                                          gsl_matrix *hess,
                                          gsl_vector *xx_hess,
                                          gsl_vector *yy_hess,
                                          gsl_vector *grad_hess,
                                          gsl_matrix *hess_hess,
                                          int iteration,
                                          double step_size, double step_limit, int max_iterations);

void trust_region_step(gsl_vector *grad, gsl_matrix *hess, gsl_vector *xxStep, double rho);

ostream & operator<<(ostream& os, const gsl_vector *vec);

#endif //FILTERSQP_LM_OPTIM_H
