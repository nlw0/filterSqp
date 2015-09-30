#include <fstream>
#include <iostream>
#include <iomanip>
#include "papersheet.h"

#include "tr_optim.h"

using namespace std;


int points = 11;
int dims = 2;

void my_target_papersheet_hess(double *x_val, double *y_val, double *gradient, double *hessian) {
    return target_papersheet_hess(lines, columns, x_val, y_val, gradient, hessian);
}


int main(int argc, const char *argv[]) {

    int n_vars = points * dims;
    double xx[n_vars];

    double seglen = 0.3;

    int i, j;
    for (i = 0; i < lines; i++)
        for (j = 0; j < columns; j++) {
            int pt = 3 * (i * columns + j);
            xx[pt + 0] = seglen * (j - columns / 2);
            xx[pt + 1] = seglen * (lines / 2 - i);
            xx[pt + 2] = -1.0;
        }

    for (i = 0; i < lines; i++)
        for (j = 0; j < columns; j++) {
            int pt = 3 * (i * columns + j);
            cout << xx[pt + 0] << "\t" << xx[pt + 1] << "\t" << xx[pt + 2] << endl;
        }

    // Read initial position.
    //    fstream fs;
    //    fs.open(argv[1], std::fstream::in);
    //    fs >> xx[0] >> xx[1];

    // trust_region_optimization(target_rosenbrock_hess, 2, xx, 0.1);

    double yy[1];
    double gg[lines * columns * dims];
    double hh[lines * columns * dims * lines * columns * dims];

    my_target_papersheet_hess(xx, yy, gg, hh);

    double *res = trust_region_optimization(my_target_papersheet_hess, n_vars, xx, 0.1, 1e-4, 500);

    for (i = 0; i < lines; i++)
        for (j = 0; j < columns; j++) {
            int pt = 3 * (i * columns + j);
            cout << res[pt + 0] << "\t" << res[pt + 1] << "\t" << res[pt + 2] << endl;
        }

    fstream fs;
    fs.open(argv[1], std::fstream::out);
    for (i = 0; i < lines; i++)
        for (j = 0; j < columns; j++) {
            int pt = 3 * (i * columns + j);
            fs << res[pt + 0] << "\t" << res[pt + 1] << "\t" << res[pt + 2] << endl;
        }


    return 0;
}
