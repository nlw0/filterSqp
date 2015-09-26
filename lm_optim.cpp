#include <gsl/gsl_vector_double.h>
#include "lm_optim.h"

using namespace std;

double *levenberg_marquadt(void (*function)(double *, double *, double *, double *), int dim, double xx0[],
                           double step_size, double step_limit, int max_iterations) {
    cout << setprecision(15);

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

    rec_levenberg_marquadt(function, xx, grad, hess, xxNew, 0, step_size, step_limit, max_iterations);

    gsl_vector_free(grad);
    gsl_vector_free(hs);
    gsl_vector_free(hwrk);
    gsl_vector_free(xxNew);
    gsl_matrix_free(hess);
    gsl_matrix_free(hv);

    return xx->data;
}


bool vector_too_small(gsl_vector *vec, double tol = 1e-10) {
    double vmin, vmax;
    gsl_vector_minmax(vec, &vmin, &vmax);
    return vmin > -tol && vmax < tol;
};

void find_step(gsl_vector *grad,
               gsl_matrix *hess,
               gsl_vector *xxStep,
               double rho) {

    const size_t dim = hess->size1;

    gsl_vector *eval = gsl_vector_alloc(dim);
    gsl_matrix *evec = gsl_matrix_alloc(dim, dim);
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(dim);

    gsl_eigen_symmv(hess, eval, evec, w);

    if (vector_too_small(eval)) {
        double gradient_norm = gsl_blas_dnrm2(grad);
        if (vector_too_small(grad)) {
            // Handle creep case where both gradient and Hessian are zero.
            gsl_vector_memcpy(xxStep, &gsl_matrix_row(evec, 0).vector);
            gsl_vector_scale(xxStep, rho);
        } else {
            // Flat function with gradient, just follow the gradient.
            gsl_vector_memcpy(xxStep, grad);
            gsl_vector_scale(xxStep, rho / gradient_norm);
        }
    } else {
        gsl_vector *alpha = gsl_vector_alloc(dim);
        gsl_vector *beta = gsl_vector_alloc(dim);

        gsl_blas_dgemv(CblasTrans, 1.0, evec, grad, 0.0, alpha);

        double nu_opt = *(find_root(dim, alpha->data, eval->data, rho, 1e-5));

        bool saddle_in = false;
        double beta_norm = 0.0;

        if (vector_too_small(alpha)) {
            double min_eval = gsl_vector_min(eval);
            int i;
            for (i = 0; i < dim; i++)
                beta->data[i] = (eval->data[i] > min_eval) ? -alpha->data[i] / (eval->data[i] - min_eval) : 0;
            beta_norm = gsl_blas_dnrm2(beta);
            saddle_in = beta_norm <= rho;
        }

        if (saddle_in) {
            gsl_vector *mv = &gsl_matrix_column(evec, gsl_vector_min_index(eval)).vector;
            gsl_vector_memcpy(xxStep, mv);
            gsl_blas_dgemv(CblasNoTrans, 1.0, evec, beta, sqrt(rho * rho - beta_norm * beta_norm), xxStep);
        } else {
            // cout << "***" << endl;
            // cout << "grad: " << grad->data[0] << " " << grad->data[1] << endl;
            // cout << "hess: " << hess->data[0] << " " << hess->data[1] << " " << hess->data[2] << " " << hess->data[3] << endl;
            // cout << "eval: " << eval->data[0] << " " << eval->data[1] << endl;
            // cout << "evec: " << evec->data[0] << " " << evec->data[1] << " " << evec->data[2] << " " << evec->data[3] << endl;
            // cout << "alpha: " << alpha->data[0] << " " << alpha->data[1] << endl;
            // cout << "nu_opt: " << nu_opt << endl;

            gsl_vector_memcpy(beta, alpha);
            gsl_vector_scale(beta, 1.0);
            gsl_vector_add_constant(eval, nu_opt);
            gsl_vector_div(beta, eval);
            gsl_blas_dgemv(CblasNoTrans, 1.0, evec, beta, 0.0, xxStep);

        }

        gsl_vector_free(alpha);
        gsl_matrix_free(beta);
    }

    gsl_eigen_symmv_free(w);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
}

gsl_vector *rec_levenberg_marquadt(void (*function)(double *, double *, double *, double *),
                                   gsl_vector *xx,
                                   gsl_vector *grad,
                                   gsl_matrix *hess,

                                   gsl_vector *xxStep,

                                   int iteration,
                                   double step_size, double step_limit, int max_iterations) {

    double yy;
    function(xx->data, &yy, grad->data, hess->data);

    find_step(grad, hess, xxStep, step_size);
    // cout << "**" << step_size << " - " << delta << " - " << lambda << endl;

    //    cout << iteration << "\t" << xx->data[0] << "\t" << xx->data[1] << "\t" << yy
    //    << "\t" << grad->data[0] << "\t" << grad->data[1] << endl;
    cout << iteration << "\t" << xx->data[0] << "\t" << xx->data[1] << "\t" << yy
    << "\t" << xxStep->data[0] << " " << xxStep->data[1] << endl;

    gsl_vector_scale(xxStep, -1.0);
    gsl_vector_add(xx, xxStep);


    double cc = max(gsl_vector_max(grad), -gsl_vector_min(grad));

    bool stop_criterion = iteration >= max_iterations || cc < step_limit;

    if (stop_criterion)
        return xx;
    else
        return rec_levenberg_marquadt(function, xx, grad, hess, xxStep, iteration + 1, step_size,
                                      step_limit, max_iterations);
}
