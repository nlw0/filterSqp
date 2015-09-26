//
// Created by nlw on 26/09/15.
//

#ifndef FILTERSQP_HIMMELBLAU_H
#define FILTERSQP_HIMMELBLAU_H

#include "fadbad/badiff.h"
#include "fadbad/fadiff.h"
using namespace fadbad;

B<F<double>> sq(const B<F<double>> &x);

B<F<double>> himmelblau_function(B<F<double>> x[2]);

void target_himmelblau_hess(double *x_val, double *y_val, double *gradient, double *hessian);

#endif //FILTERSQP_HIMMELBLAU_H
