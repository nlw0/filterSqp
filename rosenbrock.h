#ifndef FILTERSQP_ROSENBROCK_H
#define FILTERSQP_ROSENBROCK_H

#include "rosenbrock.h"
#include "fadbad/badiff.h"
#include "fadbad/fadiff.h"
using namespace fadbad;

F<double> sq(const F<double> &x);
F<double> rosenbrock_function(F<double> x[2]);
void target_rosenbrock_grad(double *x_val, double *y_val, double *dy_dx);
void target_rosenbrock_hess(double *x_val, double *y_val, double *gradient, double *hessian);

#endif //FILTERSQP_ROSENBROCK_H
