#ifndef LM_H_
#define LM_H_

#include <float.h>
#include <gsl/gsl_blas.h>

/*
 * Ordinary Least Square: all matrices are transposed as input
 * n: data points, p: variables, k: regressions
 * X_t: p*n, Y_t: k*n, beta_t, sderr_t, t_t, pv_t: k*p
*/
#define OLS_SUCCESS 0
#define OLS_FAIL -1
#define OLS_COLINEAR -2



// workspace for OLS regression
// X: n*p, Y:n*k
typedef struct OLS_space
{
	size_t p, n, k;

	gsl_matrix *I, *T, *e_t;

} OLS_space;

OLS_space * OLS_space_alloc(const size_t p, const size_t n, const size_t k);

// resize OLS space by simply changing the dimension
void OLS_space_resize(OLS_space *r, const size_t p, const size_t n, const size_t k);

void OLS_space_free(OLS_space *r);



int OLS_matrix(
		OLS_space *space,

		const gsl_matrix *X_t,
		const gsl_matrix *Y_t,

		gsl_matrix *beta_t,
		gsl_matrix *sderr_t,
		gsl_matrix *t_t,
		gsl_matrix *pv_t,

		gsl_vector *sigma2_v);

int OLS_vector(
		OLS_space *space,

		const gsl_matrix *X_t,
		const gsl_vector *Y,

		gsl_vector *beta,
		gsl_vector *sderr,
		gsl_vector *t,
		gsl_vector *pv,
		double *sigma2);


/*
 * Frisch Waugh Lovell theorem to test R effect with background B
 * Regression: X = [B R] ~ Y
 * All matrices are transposed as input
 * n: data points, p: variables, m: X factors, k: regressions
 * B_t: p*n, R_t: m*n, Y_t: k*n, beta_t, sderr_t, t_t, pv_t, RSS_t: k*m
*/

// workspace for FWL regression
typedef struct FWL_space
{
	size_t p, n, m, k;

	gsl_matrix *I, *D, *gamma_1, *gamma_2, *f_t, *g_t;
	gsl_vector *f_nrm2, *g_nrm2;

} FWL_space;

FWL_space * FWL_space_alloc(size_t p, size_t n, size_t m, size_t k, gsl_matrix *g_t);

void FWL_space_free(FWL_space *r);



/*
 *	If any column of R matrix is collinear with B, p-value: 1 and beta: inf is set.
 * */

int FWL_matrix(
		// workspace
		FWL_space *space,

		// 1: calculate all workspace variables
		// 0: use previous existing values for background variables
		const int calculate,

		// 1 do t-test, 0: no t-test
		const int t_test,

		// input
		const gsl_matrix *B_t,
		const gsl_matrix *R_t,
		const gsl_matrix *Y_t,

		// output
		gsl_matrix *beta_t,
		gsl_matrix *sderr_t,
		gsl_matrix *t_t,
		gsl_matrix *pv_t,
		gsl_matrix *RSS_t);


int FWL_vector(
		// workspace
		FWL_space *space,

		// 1: calculate all workspace variables
		// 0: use previous existing values for background variables
		const int calculate,

		// 1 do t-test, 0: no t-test
		const int t_test,

		// input
		const gsl_matrix *B_t, const gsl_matrix *R_t, const gsl_vector *Y,

		// output
		gsl_vector *beta, gsl_vector *sderr, gsl_vector *t, gsl_vector *pv, gsl_vector *RSS);


// workspace for FWL regression
typedef struct Forward_space
{
	gsl_matrix
		*X_t,	// combined matrix space
		*U_t	;	// unselected

	gsl_vector *Y, *product;

} Forward_space;

Forward_space * Forward_space_alloc(size_t p, size_t n, gsl_matrix *X_t);

void Forward_space_free(Forward_space *r);

// resize the forward selection workspace.
void Forward_space_resize(Forward_space *r, const size_t p, const size_t n);

int forward_selection(Forward_space *space, size_t index_array[], double metric_array[], size_t stop_cnt);


#endif /* LM_H_ */
