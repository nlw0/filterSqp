#ifndef FILTERSQP_ROSENBROCK_H
#define FILTERSQP_ROSENBROCK_H

adept::adouble sq(adept::adouble x);
adept::adouble rosenbrock_function(adept::adouble x[2]);
void target_rosenbrock(double x_val[2], double y_val[1], double dy_dx[2]);

#endif //FILTERSQP_ROSENBROCK_H
