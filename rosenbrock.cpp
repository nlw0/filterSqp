#include "rosenbrock.h"

F<double> sq(const F<double> &x) {
    return x * x;
}

F<double> rosenbrock_function(F<double> x[2]) {
    F<double> a = 1.0;
    F<double> b = 100.0;
    return sq(a - x[0]) + b * sq(x[1] - sq(x[0]));
}

void target_rosenbrock(double x_val[2], double y_val[1], double dy_dx[2]) {
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
