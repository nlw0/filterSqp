#include <iostream>
#include <fstream>
#include "adept.h"

using namespace std;
using namespace adept;


adouble sq(adouble x) { return x * x; }


adouble rosenbrock_function(adouble x[2]) {
    adouble a = 1.0;
    adouble b = 100.0;
    return sq(a - x[0]) + b * sq(x[1] - sq(x[0]));
}

void target_rosenbrock(double x_val[2], double y_val[1], double dy_dx[2]) {
    adept::Stack stack;

    adouble x[2];
    adept::set_values(x, 2, x_val);

    stack.new_recording(); // Start recording

    adouble y = rosenbrock_function(x);

    y.set_gradient(1.0);
    stack.compute_adjoint();
    dy_dx[0] = x[0].get_gradient(); // Store the first gradient
    dy_dx[1] = x[1].get_gradient(); // Store the second gradient
    *y_val = y.value();
}

int main(int argc, const char *argv[]) {

    cout << "Hello, World!" << endl;

    double xx[2] = {0, 0};
    fstream fs;
    fs.open(argv[1], fstream::in);

    double yy[1];
    double grad[2];

    int n = 0;
    int i;

    fs >> n;
    for (i = 0; i < n; i++) {
        fs >> xx[0] >> xx[1];

        target_rosenbrock(xx, yy, grad);

        cout << xx[0] << " " << xx[1] << endl;
        cout << yy[0] << endl;
        cout << grad[0] << " " << grad[1] << endl;
        cout << endl;
    }
    return 0;
}
