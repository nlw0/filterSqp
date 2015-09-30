#ifndef FILTERSQP_PAPERSHEET_H
#define FILTERSQP_PAPERSHEET_H

#include "fadbad/badiff.h"
#include "fadbad/fadiff.h"
using namespace fadbad;

B<F<double>> sq(const B<F<double>> &x);
B<F<double>> papersheet_function(int lines, int columns, B<F<double>> x[]);

void target_papersheet_hess(int lines, int columns, double *x_val, double *y_val, double *gradient, double *hessian);

#endif //FILTERSQP_PAPERSHEET_H
