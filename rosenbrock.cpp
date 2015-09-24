#include "adept.h"
#include "rosenbrock.h"

adept::adouble sq(adept::adouble x) { return x * x; }

adept::adouble rosenbrock_function(adept::adouble x[2]) {
    adept::adouble a = 1.0;
    adept::adouble b = 100.0;
    return sq(a - x[0]) + b * sq(x[1] - sq(x[0]));
}

void target_rosenbrock(double x_val[2], double y_val[1], double dy_dx[2]) {
    adept::Stack stack;

    adept::adouble x[2];
    set_values(x, 2, x_val);

    stack.new_recording(); // Start recording

    adept::adouble y = rosenbrock_function(x);

    y.set_gradient(1.0);
    stack.compute_adjoint();
    dy_dx[0] = x[0].get_gradient(); // Store the first gradient
    dy_dx[1] = x[1].get_gradient(); // Store the second gradient
    *y_val = y.value();
}
