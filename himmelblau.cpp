#include "himmelblau.h"

B<F<double>> sq(const B<F<double>> &x) { return x * x; }

B<F<double>> himmelblau_function(B<F<double>> x[2]) {
    return sq(sq(x[0]) + x[1] - 11) + sq(x[0] + sq(x[1]) - 7);
}

void target_himmelblau_hess(double *x_val, double *y_val, double *gradient, double *hessian) {
    B<F<double>> x[2];
    x[0] = x_val[0];
    x[1] = x_val[1];
    x[0].x().diff(0, 2);
    x[1].x().diff(1, 2);

    B<F<double>> y = himmelblau_function(x);

    y.diff(0, 1);
    y_val[0] = y.x().x();
    gradient[0] = x[0].d(0).x();
    gradient[1] = x[1].d(0).x();

    hessian[0] = x[0].d(0).d(0);
    hessian[1] = x[0].d(0).d(1);
    hessian[2] = x[1].d(0).d(0);
    hessian[3] = x[1].d(0).d(1);
}
