#include <fstream>
#include <iostream>
#include <iomanip>
#include "rosenbrock.h"

#include "lm_optim.h"

using namespace std;


int main(int argc, const char *argv[]) {

    double xx[2] = {0, 0};
    fstream fs;
    fs.open(argv[1], std::fstream::in);

    // Read initial position.
    fs >> xx[0] >> xx[1];

    double *dd = levenberg_marquadt(target_rosenbrock_hess, 2, xx, 0.1);

    // Output solution
    cout << "sol:" << dd[0] << " " << dd[1];

    return 0;
}
