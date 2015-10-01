#include "papersheet.h"

B<F<double>> sq(const B<F<double>> &x) { return x * x; }

B<F<double>> coplanar(const B<F<double>> *a, const B<F<double>> *b, const B<F<double>> *c) {
    B<F<double>> x[3];
    x[0] = (a[1] * c[2] - a[2] * c[1]);
    x[1] = (a[2] * c[0] - a[0] * c[2]);
    x[2] = (a[0] * c[1] - a[1] * c[0]);
    B<F<double>> n2x = sq(x[0]) + sq(x[1]) + sq(x[2]);
    B<F<double>> n2b = sq(b[0]) + sq(b[1]) + sq(b[2]);

    return sq(b[0] * x[0] + b[1] * x[1] + b[2] * x[2]) / (n2x * n2b);
}


B<F<double>> papersheet_function(int lines, int columns, B<F<double>> *x) {
    int dims = 3;
    int n_vars = lines * columns * dims;

    int i, j, pt;

    B<F<double>> sphere_error = 0.0;
    for (i = 0; i < lines; i++)
        for (j = 0; j < columns; j++) {
            pt = 3 * (i * columns + j);
            double scale = 1.0;
            if (i == lines / 2 || j == columns / 2) scale = 100.0;
            sphere_error += scale * sq(sqrt(sq(x[pt]) + sq(x[pt + 1]) + sq(x[pt + 2])) - 1.0);
        }

    B<F<double>> angular_error = 0.0;
    for (i = 0; i < lines; i++)
        for (j = 1; j < columns - 1; j++) {
            pt = 3 * (i * columns + j);
            angular_error += coplanar(x + pt - 3, x + pt, x + pt + 3);
        }
    for (i = 1; i < lines - 1; i++)
        for (j = 0; j < columns; j++) {
            pt = 3 * (i * columns + j);
            angular_error += coplanar(x + pt - 3 * columns, x + pt, x + pt + 3 * columns);
        }

    B<F<double>> fixation_error = 0.0;
    i = lines / 2;
    j = columns / 2;
    pt = 3 * (i * columns + j);
    fixation_error += sq(x[pt]);
    fixation_error += sq(x[pt + 1]);
    fixation_error += sq(x[pt + 2] + 1);

    i = lines / 2;
    for (j = 0; j < columns; j++) {
        int pt = 3 * (i * columns + j);
        fixation_error += sq(x[pt + 1]);
    }
    j = columns / 2;
    for (i = 0; i < lines; i++) {
        int pt = 3 * (i * columns + j);
        fixation_error += sq(x[pt]);
    }

    B<F<double>> scale_error = 0.0;
    B<F<double>> lenA = 0.0;
    B<F<double>> lenB = 0.0;
    double seglen = 0.3;
    i = lines / 2;
    for (j = 0; j < columns - 1; j++) {
        int pta = 3 * (i * columns + j);
        int ptb = 3 * (i * columns + (j + 1));
        lenA += sq(sqrt(sq(x[pta] - x[ptb]) + sq(x[pta + 1] - x[ptb + 1]) + sq(x[pta + 2] - x[ptb + 2])) - seglen);
    }

    j = columns / 2;
    for (i = 0; i < lines - 1; i++) {
        int pta = 3 * (i * columns + j);
        int ptb = 3 * ((i + 1) * columns + j);
        lenB += sq(sqrt(sq(x[pta] - x[ptb]) + sq(x[pta + 1] - x[ptb + 1]) + sq(x[pta + 2] - x[ptb + 2])) - seglen);
    }

    scale_error += lenA + lenB;

    return sphere_error + 10.0 * scale_error + 1.0 * angular_error + 1000.0 * fixation_error;

}

void target_papersheet_hess(int lines, int columns, double *x_val, double *y_val, double *gradient, double *hessian) {
    int dimensions = 3;
    int n_vars = lines * columns * dimensions;

    B<F<double>> x[n_vars];

    unsigned int j;
    unsigned int i;
    for (i = 0; i < n_vars; i++) {
        x[i] = x_val[i];
        x[i].x().diff(i, (unsigned int) n_vars);
    }

    B<F<double>> y = papersheet_function(lines, columns, x);

    y.diff(0, 1);
    *y_val = y.x().x();

    for (i = 0; i < n_vars; i++)
        gradient[i] = x[i].d(0).x();

    for (i = 0; i < n_vars; i++)
        for (j = 0; j < lines * columns; j++)
            hessian[-i*j] = x[i].d(0).d(j); //WTF

}
