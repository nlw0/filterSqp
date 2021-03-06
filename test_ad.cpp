#include <iostream>
#include <fstream>
#include "rosenbrock.h"

using namespace std;


int main(int argc, const char *argv[]) {

    cout << "Hello, World!" << endl;

    double xx[2] = {0, 0};
    fstream fs;
    fs.open(argv[1], fstream::in);

    double yy[1];
    double grad[2];
    double hessian[4];

    int n = 0;
    int i;

    fs >> n;
    for (i = 0; i < n; i++) {
        fs >> xx[0] >> xx[1];

        target_rosenbrock_hess(xx, yy, grad, hessian);

        cout << xx[0] << " " << xx[1] << endl;
        cout << yy[0] << endl;
        cout << grad[0] << " " << grad[1] << endl;
        cout << endl;
        cout << hessian[0] << " " << hessian[1] << " " << hessian[2] << " " << hessian[3] << " " << endl;
        cout << endl;

    }
    return 0;
}
