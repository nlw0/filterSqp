#include <fstream>

#include "lm_optim.h"

using namespace std;


void plane_function(double *xx, double *yy, double *grad, double *hess) {
    // *yy = xx[0] + 0.1*xx[1];

    *yy = xx[0] * xx[1];

    grad[0] = xx[1];
    grad[1] = xx[0];
    hess[0] = 0.0;
    hess[1] = 1.0;
    hess[2] = 1.0;
    hess[3] = 0.0;
}

int main(int argc, const char *argv[]) {

    double xx[2] = {0, 0};

    levenberg_marquadt(plane_function, 2, xx, 0.1, 1e-5, 10);

    return 0;
}
