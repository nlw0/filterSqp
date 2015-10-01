#include "papersheet.h"
#include <iostream>

using namespace std;

B<F<double>> sq(const B<F<double>> &x) { return x * x; }

B<F<double>> bend(const B<F<double>> *a, const B<F<double>> *b, const B<F<double>> *c) {
    B<F<double>> va[2];
    B<F<double>> vb[2];

    va[0] = b[0] - a[0];
    va[1] = b[1] - a[1];
    vb[0] = c[0] - b[0];
    vb[1] = c[1] - b[1];

    B<F<double>> na2 = sq(va[0]) + sq(va[1]);
    B<F<double>> nb2 = sq(vb[0]) + sq(vb[1]);

    return sq(va[0] * vb[1] - va[1] * vb[0]) / (na2 * nb2);
}

B<F<double>> snakeA_function(int points, B<F<double>> *x) {
    int dims = 2;
    int n_vars = points * dims;

    int i, j, pt;

    B<F<double>> sphere_error = 0.0;
    for (j = 0; j < points; j++) {
        pt = dims * j;
        sphere_error += sqrt(sq(sqrt(sq(x[pt]) + sq(x[pt + 1])) - 1.0));
    }


    B<F<double>> elastic_error = 0.0;
    for (i = 1; i < points - 1; i++) {
        int ptA = dims * (i - 1);
        int ptB = dims * i;
        int ptC = dims * (i + 1);
        elastic_error += bend(x + ptA, x + ptB, x + ptC);
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

    return 1.0 * sphere_error + 0.00 * elastic_error + 0.1 * scale_error;

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
