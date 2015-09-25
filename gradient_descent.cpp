#include <iostream>
#include <iomanip>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector_double.h>

using namespace std;

#include "gradient_descent.h"

double *gradient_descent(void (*function)(double *, double *, double *), int dim, double xx0[],
                         double step_size, double steplimit, int max_iterations) {
    cout << setprecision(15);

    gsl_vector *gradV = gsl_vector_alloc((const size_t) dim);
    gsl_vector *xx = gsl_vector_alloc((const size_t) dim);
    gsl_vector_view xx0V = gsl_vector_view_array(xx0, (size_t) dim);

    gsl_vector_memcpy(xx, &xx0V.vector);

    gsl_vector *sol = rec_gradient_descent(function, xx, gradV, 0, step_size, steplimit, max_iterations);

    gsl_vector_free(gradV);
    //gsl_vector_free(xx);

    return sol->data;
}

gsl_vector *rec_gradient_descent(void (*function)(double *, double *, double *),
                                 gsl_vector *xx,
                                 gsl_vector *grad,
                                 int iteration,
                                 double step_size, double steplimit, int max_iterations) {

    double yy;
    function(xx->data, &yy, grad->data);

    cout << iteration << "\t" << xx->data[0] << "\t" << xx->data[1] << "\t" << yy
    << "\t" << grad->data[0] << "\t" << grad->data[1] << endl;

    gsl_vector_scale(grad, -step_size);
    gsl_vector_add(xx, grad);

    double cc = max(gsl_vector_max(grad), -gsl_vector_min(grad));

    bool stop_criterion = iteration >= max_iterations || cc < steplimit;

    if (stop_criterion)
        return xx;
    else
        return rec_gradient_descent(function, xx, grad, iteration + 1, step_size, steplimit, max_iterations);
}
