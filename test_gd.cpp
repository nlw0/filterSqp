#include <fstream>
#include <iostream>
#include <iomanip>
#include "rosenbrock.h"

#include "gradient_descent.h"

using namespace std;


int main(int argc, const char *argv[]) {

    double xx[2] = {0, 0};
    fstream fs;
    fs.open(argv[1], std::fstream::in);

    // Read initial position.
    fs >> xx[0] >> xx[1];

    double *dd = gradient_descent(target_rosenbrock, 2, xx);

    // Output solution
    cout << "sol:" << dd[0] << " " << dd[1];

    return 0;
}
