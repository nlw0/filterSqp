#include <gsl/gsl_vector_double.h>
#include "tr_optim.h"

using namespace std;


double *trust_region_optimization(void (*function)(double *, double *, double *, double *), int dim, double *xx0,
                                  double step_size, double step_limit, int max_iterations) {
    cout << setprecision(15);

    gsl_vector *xx = gsl_vector_alloc((const size_t) dim);
    gsl_vector *yy = gsl_vector_alloc((const size_t) 1);
    gsl_vector *grad = gsl_vector_alloc((const size_t) dim);
    gsl_matrix *hess = gsl_matrix_alloc((const size_t) dim, (const size_t) dim);


    gsl_vector_view xx0V = gsl_vector_view_array(xx0, (size_t) dim);
    gsl_vector_memcpy(xx, &xx0V.vector);

    function(xx->data, yy->data, grad->data, hess->data);
    gsl_vector *result = rec_trust_region_optimization(function, xx, yy, grad, hess, 0, step_size, step_limit,
                                                       max_iterations);

    return result->data;
}

gsl_vector *rec_trust_region_optimization(
        void (*function)(double *, double *, double *, double *),
        gsl_vector *xx,
        gsl_vector *yy,
        gsl_vector *grad,
        gsl_matrix *hess,
        int iteration,
        double step_size, double step_limit, int max_iterations) {

    int i;

    gsl_vector *xxNew = gsl_vector_alloc(xx->size);

    trust_region_step(grad, hess, xxNew, step_size);
    gsl_vector_scale(xxNew, -1.0);

    double cc = max(gsl_vector_max(xxNew), -gsl_vector_min(xxNew));

    bool stop_criterion = iteration >= (max_iterations - 1) || cc < step_limit;

    gsl_vector_add(xxNew, xx);

    if (stop_criterion) {
        cout << iteration << "\t" << yy << "\t" << xx << "\t" << xxNew << endl;
        return xxNew;

    } else {
        gsl_vector *yyNew = gsl_vector_alloc(1);
        gsl_vector *gradNew = gsl_vector_alloc(xx->size);
        gsl_matrix *hessNew = gsl_matrix_alloc(xx->size, xx->size);
        function(xxNew->data, yyNew->data, gradNew->data, hessNew->data);

        // Test if predicted function value is close enough to yyNew. Reduce step size if not.
        bool approximation_was_good = false;
        if (approximation_was_good) {
            gsl_vector_free(xxNew);
            gsl_vector_free(yyNew);
            gsl_vector_free(gradNew);
            gsl_matrix_free(hessNew);
            double step_size_new = step_size / 2.0;
            return rec_trust_region_optimization(function, xx, yy, grad, hess, iteration, step_size_new,
                                                 step_limit, max_iterations);
        } else {
            cout << iteration << "\t" << yy << "\t" << xx << "\t" << xxNew << endl;
            gsl_vector_free(xx);
            gsl_vector_free(yy);
            gsl_vector_free(grad);
            gsl_matrix_free(hess);
            return rec_trust_region_optimization(function, xxNew, yyNew, gradNew, hessNew, iteration + 1, step_size,
                                                 step_limit, max_iterations);
        }

    }

}

bool vector_too_small(gsl_vector *vec, double tol = 1e-7) {
    double vmin, vmax;
    gsl_vector_minmax(vec, &vmin, &vmax);
    return vmin > -tol && vmax < tol;
}

bool hard_case(gsl_vector *vec, gsl_vector *lam, double tol = 1e-7) {
    double lambda_one = gsl_vector_min(lam);
    if (lambda_one <= 0) {
        int i;
        for (i = 0; i < vec->size; i++) if (fabs(vec->data[i]) < tol && lam->data[i] == lambda_one) return true;
    }
    return false;
}

void trust_region_step(gsl_vector *grad,
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
            gsl_vector_view firstev = gsl_matrix_row(evec, 0);
            gsl_vector_memcpy(xxStep, &(firstev.vector));
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

        bool fork_inside_region = false;
        double beta_norm = 0.0;

        if (hard_case(alpha, eval)) {
            double min_eval = gsl_vector_min(eval);
            int i;
            for (i = 0; i < dim; i++)
                beta->data[i] = (eval->data[i] > min_eval) ? alpha->data[i] / (eval->data[i] - min_eval) : 0;
            beta_norm = gsl_blas_dnrm2(beta);
            fork_inside_region = beta_norm <= rho;
        }

        if (fork_inside_region) {
            gsl_vector_view mv = gsl_matrix_column(evec, gsl_vector_min_index(eval));
            gsl_vector_memcpy(xxStep, &mv.vector);
            gsl_blas_dgemv(CblasNoTrans, 1.0, evec, beta, sqrt(rho * rho - beta_norm * beta_norm), xxStep);
        } else {
            double nu_opt = *(find_root(dim, alpha->data, eval->data, rho, 1e-5));
            gsl_vector_memcpy(beta, alpha);
            gsl_vector_scale(beta, 1.0);
            gsl_vector_add_constant(eval, nu_opt);
            gsl_vector_div(beta, eval);
            gsl_blas_dgemv(CblasNoTrans, 1.0, evec, beta, 0.0, xxStep);
        }

        gsl_vector_free(alpha);
        gsl_vector_free(beta);
    }

    gsl_eigen_symmv_free(w);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
}

ostream &operator<<(ostream &os, const gsl_vector *vec) {
    int i;
    cout << vec->data[0];
    for (i = 1; i < vec->size; i++) cout << "\t" << vec->data[i];
    return os;
}

