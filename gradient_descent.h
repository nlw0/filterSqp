//
// Created by nlw on 24/09/15.
//

#ifndef FILTERSQP_GRADIENT_DESCENT_H
#define FILTERSQP_GRADIENT_DESCENT_H

double * gradient_descent(void (*function) (double*, double*, double*), int dim, double xx0[], double step_size);


gsl_vector *rec_gradient_descent(void (*function)(double *, double *, double *),
                                 gsl_vector *xx,
                                 gsl_vector *grad,
                                 int iteration = 0,
                                 double step_size = 0.001, double steplimit = 1e-10, int max_iterations = 10000);


#endif //FILTERSQP_GRADIENT_DESCENT_H
