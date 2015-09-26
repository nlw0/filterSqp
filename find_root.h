//
// Created by nlw on 26/09/15.
//

#ifndef FILTERSQP_FIND_ROOT_H
#define FILTERSQP_FIND_ROOT_H

#include <iostream>
#include <iomanip>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector_double.h>
#include <math.h>

#include <experimental/optional>

using namespace std;
using namespace experimental;

optional<double> find_root(int dim, double *alpha, double *lam, double rho, double rho_tol=1e-5);

optional<double> rec_find_root(int dim, double *alpha, double *lam, double rho, double rho_tol,
                               double nu_min, double nu_a, double nu_b, double f_a, double f_b, int ki = 0);

double calculate_length(int dim, const double *alpha, const double *lam, double nu);

#endif //FILTERSQP_FIND_ROOT_H
