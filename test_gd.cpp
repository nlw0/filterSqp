#include <fstream>
#include <iostream>
#include <iomanip>
#include "adept.h"
#include "rosenbrock.h"
#include <gsl/gsl_linalg.h>

#include "gradient_descent.h"

using namespace std;
using namespace adept;


int main(int argc, const char *argv[]) {

    double xx[2] = {0, 0};
    fstream fs;
    fs.open(argv[1], std::fstream::in);

    // Read initial position.
    fs >> xx[0] >> xx[1];

    double *dd = gradient_descent(target_rosenbrock, 2, xx, 0.001);

    // cout << "sol:" << dd[0] << " " << dd[1];

    return 0;
}
