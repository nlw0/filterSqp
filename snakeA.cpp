#include "papersheet.h"
#include <iostream>

using namespace std;

B<F<double>> sq(const B<F<double>> &x) { return x * x; }

B<F<double>> snakeA_function(int points, B<F<double>> *x) {
    int dims = 2;
    int n_vars = points * dims;

    int i, j, pt;

    B<F<double>> sphere_error = 0.0;
    for (j = 0; j < points; j++) {
        pt = dims * j;
        sphere_error += sqrt(sq(sqrt(sq(x[pt]) + sq(x[pt + 1])) - 1.0));
    }

    B<F<double>> fixation_error = 0.0;
    i = points / 2;
    pt = 2 * i;
    fixation_error += sq(x[pt]);

    B<F<double>> scale_error = 0.0;
    double seglen = 0.2;

    for (i = 0; i < points - 1; i++) {
        int pta = dims * i;
        int ptb = dims * (i + 1);
        scale_error += sq(sqrt(sq(x[pta] - x[ptb]) + sq(x[pta + 1] - x[ptb + 1])) - seglen);
    }

    pt = points / 2;
    B<F<double>> fix_error = sq(x[dims * pt]);

    // return 0.1 * sphere_error + 10.0 * scale_error + fix_error;
    return 10.0 * sphere_error + 0.1 * scale_error + fix_error;

}


void target_snakeA_hess(int points, double *x_val, double *y_val, double *gradient, double *hessian) {
    int dims = 2;
    int n_vars = 2 * points;

    B<F<double>> x[n_vars];

    unsigned int j;
    unsigned int i;
    for (i = 0; i < n_vars; i++) {
        x[i] = x_val[i];
        x[i].x().diff(i, (unsigned int) n_vars);
    }

    B<F<double>> y = snakeA_function(points, x);

    y.diff(0, 1);
    *y_val = y.x().x();

    for (i = 0; i < n_vars; i++)
        gradient[i] = x[i].d(0).x();

    for (i = 0; i < n_vars; i++)
        for (j = 0; j < n_vars; j++)
            hessian[i * n_vars + j] = x[i].d(0).d(j);

}
