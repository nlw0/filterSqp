#ifndef FILTERSQP_ROSENBROCK_H
#define FILTERSQP_ROSENBROCK_H

#include "rosenbrock.h"
#include "fadbad/badiff.h"
#include "fadbad/fadiff.h"
using namespace fadbad;

F<double> sq(const F<double> &x);
F<double> rosenbrock_function(F<double> x[2]);
void target_rosenbrock(double x_val[2], double y_val[1], double dy_dx[2]);

#endif //FILTERSQP_ROSENBROCK_H
