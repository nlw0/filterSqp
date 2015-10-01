#include <gsl/gsl_vector_double.h>
#include "tr_optim.h"

using namespace std;

double quadratic_eval(gsl_matrix *M, gsl_vector *g, gsl_vector *x) {

    int dim = x->size;
    int i, j;
    double result = 0.0;
    for (i = 0; i < dim; i++) {
        result += g->data[i] * x->data[i];
        for (j = 0; j < dim; j++)
            result += 0.5 * x->data[i] * M->data[i * dim + j] * x->data[j];
    }
    return result;
}

double *trust_region_optimization(void (*function)(double *, double *, double *, double *), int dim, double *xx0,
                                  double step_size, double step_limit, int max_iterations) {
    cout << setprecision(15);

    gsl_vector *xx = gsl_vector_alloc((const size_t) dim);
    gsl_vector *yy = gsl_vector_alloc((const size_t) 1);
    gsl_vector *grad = gsl_vector_alloc((const size_t) dim);
    gsl_matrix *hess = gsl_matrix_alloc((const size_t) dim, (const size_t) dim);
    gsl_vector *xx_new = gsl_vector_alloc((const size_t) dim);
    gsl_vector *yy_new = gsl_vector_alloc((const size_t) 1);
    gsl_vector *grad_new = gsl_vector_alloc((const size_t) dim);
    gsl_matrix *hess_new = gsl_matrix_alloc((const size_t) dim, (const size_t) dim);

    gsl_vector_view xx0V = gsl_vector_view_array(xx0, (size_t) dim);
    gsl_vector_memcpy(xx, &xx0V.vector);

    function(xx->data, yy->data, grad->data, hess->data);
    gsl_vector *result = rec_trust_region_optimization(function, xx, yy, grad, hess,
                                                       xx_new, yy_new, grad_new, hess_new,
                                                       0, step_size, step_limit,
                                                       max_iterations);

    double *result_copy = (double *) malloc(sizeof(double) * dim);
    int i;
    for (i = 0; i < dim; i++)
        result_copy[i] = result->data[i];

    gsl_vector_free(xx);
    gsl_vector_free(yy);
    gsl_vector_free(grad);
    gsl_matrix_free(hess);
    gsl_vector_free(xx_new);
    gsl_vector_free(yy_new);
    gsl_vector_free(grad_new);
    gsl_matrix_free(hess_new);

    return result_copy;
}

gsl_vector *rec_trust_region_optimization(
        void (*function)(double *, double *, double *, double *),
        gsl_vector *xx, gsl_vector *yy, gsl_vector *grad, gsl_matrix *hess,
        gsl_vector *xx_new, gsl_vector *yy_new, gsl_vector *grad_new, gsl_matrix *hess_new,
        int iteration, double step_size, double step_limit, int max_iterations) {

    trust_region_step(grad, hess, xx_new, step_size);

    double cc = max(gsl_vector_max(xx_new), -gsl_vector_min(xx_new));
    bool stop_criterion = iteration >= (max_iterations - 1) || cc < step_limit;

    if (stop_criterion) {
        gsl_vector_add(xx_new, xx);
        cout << iteration << "\t" << step_size << "\t" << yy << "\t" << xx_new << endl;
        return xx_new;

    } else {
        double pred = quadratic_eval(hess, grad, xx_new);
        double was_newton = gsl_blas_dnrm2(xx_new) < step_size;

        gsl_vector_add(xx_new, xx);
        function(xx_new->data, yy_new->data, grad_new->data, hess_new->data);

        double actual = yy_new->data[0] - yy->data[0];

        // Test if predicted function value is close enough to yyNew. Reduce step size if not.
        bool approximation_was_good = actual < 0 && pred < 0 && pred / actual < 1.2;

        if (approximation_was_good) {
            cout << iteration << "\t" << step_size << "\t" << yy << "\t" << xx << endl;

            double step_size_new = step_size * 1.5;
            return rec_trust_region_optimization(function, xx_new, yy_new, grad_new, hess_new, xx, yy, grad, hess,
                                                 iteration + 1, step_size_new, step_limit, max_iterations);
        } else {
            double step_size_new = step_size / 4.0;
            return rec_trust_region_optimization(function, xx, yy, grad, hess, xx_new, yy_new, grad_new, hess_new,
                                                 iteration, step_size_new, step_limit, max_iterations);
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

    gsl_vector_scale(xxStep, -1.0);

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

