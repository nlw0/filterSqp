#include <iostream>
#include <fstream>
#include "adept.h"
#include "rosenbrock.h"

using namespace std;
using namespace adept;


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
