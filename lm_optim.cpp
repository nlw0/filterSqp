#include <iostream>
#include <iomanip>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector_double.h>
#include <math.h>

#include "lm_optim.h"
#include "find_root.h"

using namespace std;

double *levenberg_marquadt(void (*function)(double *, double *, double *, double *), int dim, double xx0[],
                           double step_size, double step_limit, int max_iterations) {
    cout << setprecision(5);

    gsl_vector *xx = gsl_vector_alloc((const size_t) dim);
    gsl_vector *grad = gsl_vector_alloc((const size_t) dim);
    gsl_matrix *hess = gsl_matrix_alloc((const size_t) dim, (const size_t) dim);

    gsl_matrix *hu = gsl_matrix_alloc((const size_t) dim, (const size_t) dim);
    gsl_matrix *hv = gsl_matrix_alloc((const size_t) dim, (const size_t) dim);
    gsl_vector *hs = gsl_vector_alloc((const size_t) dim);
    gsl_vector *hwrk = gsl_vector_alloc((const size_t) dim);

    gsl_vector *xxNew = gsl_vector_alloc((const size_t) dim);

    gsl_vector_view xx0V = gsl_vector_view_array(xx0, (size_t) dim);
    gsl_vector_memcpy(xx, &xx0V.vector);

    rec_levenberg_marquadt(function, xx, grad, hess, hu, hv, hs, hwrk, xxNew, 0, step_size, step_limit, max_iterations);

    gsl_vector_free(grad);
    gsl_vector_free(hs);
    gsl_vector_free(hwrk);
    gsl_vector_free(xxNew);
    gsl_matrix_free(hess);
    gsl_matrix_free(hv);

    return xx->data;
}


void rec_findstep(gsl_vector *grad,
                  gsl_matrix *hess,

                  gsl_matrix *Hu,
                  gsl_matrix *Hv,
                  gsl_vector *Hs,
                  gsl_vector *Hwrk,

                  gsl_vector *xxStep,

                  int iteration = 0,
                  double step_size = 0.1, double step_limit = 0.01, int max_iterations = 10000,
                  double lambda = 0.0,
                  double delta = 0.0) {

    do {
        gsl_matrix_memcpy(Hu, hess);
        gsl_vector_view hess_diag = gsl_matrix_diagonal(Hu);
        gsl_vector_add_constant(&hess_diag.vector, lambda);
        gsl_linalg_SV_decomp(Hu, Hv, Hs, Hwrk);
        gsl_linalg_SV_solve(Hu, Hv, Hs, grad, xxStep);
        lambda += lambda * 2 + 1.0;
        delta = gsl_blas_dnrm2(xxStep);
    } while (delta > step_size);
}

gsl_vector *rec_levenberg_marquadt(void (*function)(double *, double *, double *, double *),
                                   gsl_vector *xx,
                                   gsl_vector *grad,
                                   gsl_matrix *hess,

                                   gsl_matrix *Hu,
                                   gsl_matrix *Hv,
                                   gsl_vector *Hs,
                                   gsl_vector *Hwrk,

                                   gsl_vector *xxStep,

                                   int iteration,
                                   double step_size, double step_limit, int max_iterations) {

    double yy;
    function(xx->data, &yy, grad->data, hess->data);

    double lambda = 0.0;
    double delta = 0.0;

    rec_findstep(grad, hess, Hu, Hv, Hs, Hwrk, xxStep);
    // cout << "**" << step_size << " - " << delta << " - " << lambda << endl;

//    cout << iteration << "\t" << xx->data[0] << "\t" << xx->data[1] << "\t" << yy
//    << "\t" << grad->data[0] << "\t" << grad->data[1] << endl;
    cout << iteration << "\t" << xx->data[0] << "\t" << xx->data[1] << "\t" << yy
    << "\t" << xxStep->data[0] << " " << xxStep->data[1] << endl;

    // gsl_vector_scale(grad, -step_size);
    // gsl_vector_add(xx, grad);

    gsl_vector_scale(xxStep, -1.0);
    gsl_vector_add(xx, xxStep);


    double cc = max(gsl_vector_max(grad), -gsl_vector_min(grad));

    bool stop_criterion = iteration >= max_iterations || cc < step_limit;

    if (stop_criterion)
        return xx;
    else
        return rec_levenberg_marquadt(function, xx, grad, hess, Hu, Hv, Hs, Hwrk, xxStep, iteration + 1, step_size,
                                      step_limit, max_iterations);
}
