#include "rosenbrock.h"

F<double> sq(const F<double> &x) {
    return x * x;
}

F<double> rosenbrock_function(F<double> x[2]) {
    F<double> a = 1.0;
    F<double> b = 100.0;
    return sq(a - x[0]) + b * sq(x[1] - sq(x[0]));
}

B<F<double>> sq2(const B<F<double>> &x) {
    return x * x;
}

B<F<double>> rosenbrock_function2(B<F<double>> x[2]) {
    B<F<double>> a = 1.0;
    B<F<double>> b = 100.0;
    return sq2(a - x[0]) + b * sq2(x[1] - sq2(x[0]));
}

void target_rosenbrock_grad(double *x_val, double *y_val, double *dy_dx) {
    F<double> x[2];
    x[0] = x_val[0];
    x[1] = x_val[1];
    x[0].diff(0, 2);
    x[1].diff(1, 2);

    F<double> y = rosenbrock_function(x);

    y_val[0] = y.x();
    dy_dx[0] = y.d(0);
    dy_dx[1] = y.d(1);
}

void target_rosenbrock_hess(double *x_val, double *y_val, double *gradient, double *hessian) {
    B<F<double>> x[2];
    x[0] = x_val[0];
    x[1] = x_val[1];
    x[0].x().diff(0, 2);
    x[1].x().diff(1, 2);

    B<F<double>> y = rosenbrock_function2(x);

    y.diff(0, 1);
    y_val[0] = y.x().x();
    gradient[0] = x[0].d(0).x();
    gradient[1] = x[1].d(0).x();

    hessian[0] = x[0].d(0).d(0);
    hessian[1] = x[0].d(0).d(1);
    hessian[2] = x[1].d(0).d(0);
    hessian[3] = x[1].d(0).d(1);
}
