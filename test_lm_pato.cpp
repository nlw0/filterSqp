#include <fstream>

#include "tr_optim.h"

using namespace std;


void saddle_function(double *xx, double *yy, double *grad, double *hess) {
    // *yy = xx[0] + 0.1*xx[1];
    *yy = xx[0] * xx[0] - xx[1] * xx[1];

    grad[0] = xx[0];
    grad[1] = -xx[1];
    hess[0] = 1.0;
    hess[1] = 0.0;
    hess[2] = 0.0;
    hess[3] = -1.0;
}

void bowl_function(double *xx, double *yy, double *grad, double *hess) {
    *yy = 2 * xx[0] * xx[0] * xx[1] * xx[1];

    grad[0] = 2 * xx[0];
    grad[1] = xx[1];
    hess[0] = 2.0;
    hess[1] = 0.0;
    hess[2] = 0.0;
    hess[3] = 1.0;
}

void biggie_function(double *xx, double *yy, double *grad, double *hess) {
    *yy = 2 * xx[0] * xx[0] * xx[1] * xx[1];

    grad[0] = -10;
    grad[1] = 6;
    grad[2] = 10;
    grad[3] = 9;

    hess[0] = -10;
    hess[1] = 0;
    hess[2] = 2;
    hess[3] = -2;
    hess[4] = 0;
    hess[5] = -3;
    hess[6] = -4;
    hess[7] = 2;
    hess[8] = 2;
    hess[9] = -4;
    hess[10] = 7;
    hess[11] = -4;
    hess[12] = -2;
    hess[13] = 2;
    hess[14] = -4;
    hess[15] = 1;
}


int main(int argc, const char *argv[]) {

    double xx[4] = {2, 0.0, 0, 0};

    double rho = 0.0;
    double rho_incr = 0.05;
    for (rho = rho_incr; rho < 2.1; rho += rho_incr)
        trust_region_optimization(saddle_function, 2, xx, rho, 1e-5, 1);
//        trust_region_optimization(biggie_function, 4, xx, rho, 1e-5, 1);
    return 0;
}
