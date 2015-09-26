#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include "find_root.h"


using namespace std;


int main(int argc, const char *argv[]) {

    vector<double> xx(10);
    vector<double> yy(10);
    fstream fs;
    fs.open(argv[1], std::fstream::in);

    // Read initial position.
    double rho;
    int n, i;
    fs >> rho >> n;
    for (i = 0; i < n; i++) fs >> xx[i];
    for (i = 0; i < n; i++) fs >> yy[i];

    experimental::optional<double> dd = find_root(n, xx.data(), yy.data(), rho, 1e-10);

    // Output solution
    cout << "sol:" << *dd;

    return 0;
}
