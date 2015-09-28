//
// Created by nlw on 26/09/15.
//

#include "find_root.h"

optional<double> find_root(int dim, double *alpha, double *lam, double rho, double rho_tol) {

    double ln = lam[0];
    int i;
    for (i = 1; i < dim; i++) if (lam[i] < ln) ln = lam[i];

    // The lower bound for nu, and the initialization of the two points from the secant method.
    double nu_min, nu_a, nu_b;
    nu_min = (ln < 0) ? -ln : 0.0;
    nu_a = (ln > 0) ? 0.0 : (ln < 0) ? nu_min : 1.0;
    nu_b = (ln > 0) ? ln : (ln < 0) ? nu_min * 2.0 : 2.0;

    if (ln >= 0.0 && calculate_length(dim, alpha, lam, nu_min) < rho)
        return optional<double>(nu_min); // Perform Newton step.
    else {
        double f_a = calculate_length(dim, alpha, lam, nu_min + nu_a) - rho;
        double f_b = calculate_length(dim, alpha, lam, nu_min + nu_b) - rho;
        return rec_find_root(dim, alpha, lam, rho, rho_tol, nu_min, nu_a, nu_b, f_a, f_b);
    }
}

optional<double> rec_find_root(int dim, double *alpha, double *lam, double rho, double rho_tol,
                               double nu_min, double nu_a, double nu_b, double f_a, double f_b, int ki) {

    int MAX_ITERATIONS = 200;

    double nu_new, f_new;

//    cout << " xx " << ki << "\t" << rho << "\t" << nu_min << "\t" << nu_a << "\t" << nu_b << "\t" << f_a << "\t" << f_b << endl;

    if (ki >= (MAX_ITERATIONS - 1) || isinf(f_a) || isinf(f_b)) return nullopt;
    else if (f_a >= 0 && f_b <= 0) {

        double nu_secant = nu_b - f_b * (nu_b - nu_a) / (f_b - f_a);
        double nu_bisection = (nu_a + nu_b) / 2.0;

        nu_new = nu_secant < (nu_a * 0.8 + nu_b * 0.2) || (nu_a * 0.2 + nu_b * 0.8) < nu_secant ?
                 nu_bisection : nu_secant;

        f_new = calculate_length(dim, alpha, lam, nu_min + nu_new) - rho;

        return fabs(f_new) < rho_tol ? optional<double>(nu_min + nu_new) :
               f_new > 0 ? rec_find_root(dim, alpha, lam, rho, rho_tol, nu_min, nu_new, nu_b, f_new, f_b, ki + 1)
                         : rec_find_root(dim, alpha, lam, rho, rho_tol, nu_min, nu_a, nu_new, f_a, f_new, ki + 1);

    } else if (f_a > 0 && f_b > 0 || f_a < 0 && f_b < 0) {
        // If root is outside the interval, look to the adjoining interval.
        bool move_right = fabs(f_a) > fabs(f_b);

        nu_new = move_right ? nu_b * 2.0 : nu_a / 2.0;

        f_new = calculate_length(dim, alpha, lam, nu_min + nu_new) - rho;

        return move_right ?
               rec_find_root(dim, alpha, lam, rho, rho_tol, nu_min, nu_b, nu_new, f_b, f_new, ki + 1) :
               rec_find_root(dim, alpha, lam, rho, rho_tol, nu_min, nu_new, nu_a, f_new, f_a, ki + 1);

    } else return nullopt; // Function non-decreasing.

}

// Calculates a step size for different nu values.
double calculate_length(int dim, const double *alpha, const double *lam, double nu) {
    double output = 0.0;

    double y = 0.0;
    int i;
    for (i = 0; i < dim; ++i) {
        y = alpha[i] / (lam[i] + nu);
        output += y * y;
    }
    return sqrt(output);
}
