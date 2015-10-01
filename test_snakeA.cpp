#include <fstream>
#include <iostream>
#include <iomanip>
#include "snakeA.h"

#include "tr_optim.h"

using namespace std;


int points = 11;
int dims = 2;

void my_target_snakeA_hess(double *x_val, double *y_val, double *gradient, double *hessian) {
    return target_snakeA_hess(points, x_val, y_val, gradient, hessian);
}


int main(int argc, const char *argv[]) {

    int n_vars = points * dims;
    double xx[n_vars];

    double seglen = 0.2;

    int i, j;
    for (i = 0; i < points; i++) {
        int pt = 2 * i;
        xx[pt] = (0.1+seglen) * (i - points / 2);
        xx[pt + 1] = 0.5;
    }

    for (i = 0; i < points; i++) {
        int pt = 2 * i;
        cout << "x:" << xx[pt + 0] << "\t" << xx[pt + 1] << endl;
    }

    double yy[1];
    double gg[dims * points];
    double hh[dims * points * dims * points];

    my_target_snakeA_hess(xx, yy, gg, hh);

    double *res = trust_region_optimization(my_target_snakeA_hess, n_vars, xx, 0.1, 1e-4, 50);

    for (i = 0; i < points; i++) {
        int pt = 2 * i;
        cout << res[pt + 0] << "\t" << res[pt + 1] << endl;
    }

    fstream fs;
    fs.open(argv[1], std::fstream::out);
    for (i = 0; i < points; i++) {
        int pt = 2 * i;
        fs << res[pt + 0] << "\t" << res[pt + 1] << endl;
    }

    return 0;
}
