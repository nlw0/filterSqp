//
// Created by nlw on 30/09/15.
//

#ifndef FILTERSQP_SNAKEA1_H
#define FILTERSQP_SNAKEA1_H

#include "fadbad/badiff.h"
#include "fadbad/fadiff.h"
using namespace fadbad;

B<F<double>> snakeA_function(int points, B<F<double>> *x);

void target_snakeA_hess(int points, double *x_val, double *y_val, double *gradient, double *hessian);

#endif //FILTERSQP_SNAKEA1_H
