#include "lm.h"
#include "lin.h"
#include "gsl_util.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>

extern int verbose_flag;

void OLS_space_resize(OLS_space *r, const size_t p, const size_t n, const size_t k)
{
	if(p>=n)
	{
		fprintf(stderr, "Bad OLS: p %lu >= n %lu\n", p, n);
		exit(1);
	}

	r->p = p;
	r->n = n;
	r->k = k;

	resize_matrix(r->I, p,p, NULL);
	resize_matrix(r->T, p,k, NULL);
	resize_matrix(r->e_t, k,n, NULL);
}


// allocate the internal variable space
OLS_space * OLS_space_alloc(const size_t p, const size_t n, const size_t k)
{
	OLS_space *r = (OLS_space*)malloc(sizeof(OLS_space));

	if(p>=n)
	{
		fprintf(stderr, "Bad OLS: p %lu >= n %lu\n", p, n);
		exit(1);
	}

	r->n = n;
	r->p = p;
	r->k = k;

	r->I = gsl_matrix_alloc(p, p);
	r->T = gsl_matrix_alloc(p, k);
	r->e_t = gsl_matrix_alloc(k, n);

	// handle errors by non-exit procedure
	gsl_set_error_handler_off();

	return r;
}


void OLS_space_free(OLS_space *r)
{
	gsl_matrix_free(r->I);
	gsl_matrix_free(r->T);
	gsl_matrix_free(r->e_t);
	free(r);
}

// X_t:p*n, Y_t: k*n, beta_t:k*p
int OLS_matrix(
		OLS_space *space,

		const gsl_matrix *X_t,
		const gsl_matrix *Y_t,

		gsl_matrix *beta_t,
		gsl_matrix *sderr_t,
		gsl_matrix *t_t,
		gsl_matrix *pv_t,

		gsl_vector *sigma2_v)
{
	double sigma, t;

	// hook up workspace
	size_t p=space->p, n=space->n, k=space->k, degree=n-p, i, j;

	gsl_matrix *I = space->I, *T = space->T, *e_t = space->e_t;

	if(p>=n)
	{
		fprintf(stderr, "bad OLS: p %lu >= n %lu\n", p, n);
		return OLS_FAIL;
	}

	check_matrix_dimension(Y_t, k, n, "Y", 1);
	check_matrix_dimension(beta_t, k, p, "beta", 1);
	check_matrix_dimension(sderr_t, k, p, "sderr", 1);
	check_matrix_dimension(t_t, k, p, "t", 1);
	check_matrix_dimension(pv_t, k, p, "pv", 1);

	check_vector_dimension(sigma2_v, k, "RSS", 1);

	// I = X'X
	gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1, X_t, 0, I);

	if(gsl_linalg_cholesky_decomp(I) == GSL_EDOM)
	{
		fprintf(stderr, "Cholesky decomposition failed on X'X.\n");
		return OLS_COLINEAR;
	}

	// I = (X'X)^-1
	gsl_linalg_cholesky_invert(I);

	// T = X'Y
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, X_t, Y_t, 0, T);

	// beta = (X'X)^-1 X'Y
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, T, I, 0, beta_t);

	// e = X beta - Y
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, beta_t, X_t, 0, e_t);
	gsl_matrix_sub(e_t, Y_t);

	// calculate standard error
	gsl_vector_const_view diag = gsl_matrix_const_diagonal(I);

	for (i=0;i<k;i++)
	{
		gsl_vector_const_view e = gsl_matrix_const_row(e_t, i);

		gsl_blas_ddot(&e.vector, &e.vector, &sigma);
		sigma /= degree;

		gsl_vector_set(sigma2_v, i, sigma);

		gsl_vector_view se = gsl_matrix_row(sderr_t, i);

		for(j=0;j<p;j++)
		{
			t = sigma * gsl_vector_get(&diag.vector, j);

			if(fabs(t)<EPS)
			{
				if(verbose_flag)
					fprintf(stderr, "Warning: standard error %e smaller than %e\n", t, EPS);

				t = EPS;
			}else{
				t = sqrt(t);
			}

			gsl_vector_set(&se.vector, j, t);
		}
	}

	// t-test
	for(i=0;i<k;i++)
	{
		for(j=0;j<p;j++)
		{
			t = gsl_matrix_get(beta_t, i, j) / gsl_matrix_get(sderr_t, i, j);
			gsl_matrix_set(t_t, i, j, t);

			// two tailed test
			t = 2*gsl_cdf_tdist_Q(fabs(t), degree);
			gsl_matrix_set(pv_t, i, j, t);
		}
	}

	return OLS_SUCCESS;
}

int OLS_vector(
		OLS_space *space,

		const gsl_matrix *X_t,
		const gsl_vector *Y,

		gsl_vector *beta,
		gsl_vector *sderr,
		gsl_vector *t,
		gsl_vector *pv,

		double *sigma2)
{
	// make vectors to matrix views and call matrix version OLS

	size_t p = beta->size;

	check_vector_dimension(beta, p, "beta", 1);
	check_vector_dimension(sderr, p, "sderr", 1);
	check_vector_dimension(t, p, "t", 1);
	check_vector_dimension(pv, p, "pv", 1);

	gsl_matrix_const_view
		Y_t = gsl_matrix_const_view_vector(Y, 1, Y->size);

	gsl_matrix_view
		beta_t = gsl_matrix_view_vector(beta, 1, p),
		sderr_t = gsl_matrix_view_vector(sderr, 1, p),
		t_t = gsl_matrix_view_vector(t, 1, p),
		pv_t = gsl_matrix_view_vector(pv, 1, p);

	gsl_vector_view sigma_v = gsl_vector_view_array(sigma2, 1);

	return OLS_matrix(space, X_t, &Y_t.matrix,
			&beta_t.matrix, &sderr_t.matrix, &t_t.matrix, &pv_t.matrix,
			&sigma_v.vector);
}



// allocate the internal variable space
FWL_space * FWL_space_alloc(size_t p, size_t n, size_t m, size_t k, gsl_matrix *g_t)
{
	FWL_space *r = (FWL_space*)malloc(sizeof(FWL_space));

	if(p+1>=n)
	{
		fprintf(stderr, "Bad FWL OLS: p+1 %lu >= n %lu\n", p+1, n);
		exit(1);
	}

	r->p = p;
	r->n = n;
	r->m = m;
	r->k = k;

	r->I = gsl_matrix_alloc(p, p);
	r->D = gsl_matrix_alloc(p, n);
	r->gamma_1 = gsl_matrix_alloc(p,k);
	r->gamma_2 = gsl_matrix_alloc(p,m);
	r->f_t = gsl_matrix_alloc(k, n);

	// g_t dimension is m x n, the largest matrix when m is very huge. Share it with external matrix
	r->g_t = g_t;

	r->f_nrm2 = gsl_vector_alloc(k);
	r->g_nrm2 = gsl_vector_alloc(m);

	// handle errors by non-exit procedure
	gsl_set_error_handler_off ();

	return r;
}


void FWL_space_free(FWL_space *r)
{
	gsl_matrix_free(r->I);
	gsl_matrix_free(r->D);
	gsl_matrix_free(r->gamma_1);
	gsl_matrix_free(r->gamma_2);
	gsl_matrix_free(r->f_t);
	gsl_vector_free(r->f_nrm2);
	gsl_vector_free(r->g_nrm2);

	free(r);
}



// B_t: p*n, R_t: m*n, Y_t: k*n, beta_t, sderr_t, t_t, pv_t: k*m

int FWL_matrix( FWL_space *space,
		const int calculate,
		const int t_test,

		const gsl_matrix *B_t,
		const gsl_matrix *R_t,
		const gsl_matrix *Y_t,

		gsl_matrix *beta_t,
		gsl_matrix *sderr_t,
		gsl_matrix *t_t,
		gsl_matrix *pv_t,
		gsl_matrix *RSS_t)
{
	double se, t, denom;

	// hook up to workspace
	size_t p = space->p, n= space->n, m = space->m, k = space->k, degree=n-p-1, i, j, colinear=0;

	gsl_matrix *I = space->I, *D = space->D,
			*gamma_1 = space->gamma_1, *gamma_2 = space->gamma_2,
			*f_t = space->f_t, *g_t = space->g_t;

	gsl_vector *f_nrm2 = space->f_nrm2, *g_nrm2 = space->g_nrm2;

	// make sure dimension
	check_matrix_dimension(B_t, p, n, "B", 1);
	check_matrix_dimension(R_t, m, n, "R", 1);
	check_matrix_dimension(Y_t, k, n, "Y", 1);
	check_matrix_dimension(beta_t, k, m, "beta", 1);
	check_matrix_dimension(sderr_t, k, m, "sderr", 1);
	check_matrix_dimension(t_t, k, m, "t", 1);
	check_matrix_dimension(pv_t, k, m, "pv", 1);
	check_matrix_dimension(RSS_t, k, m, "pv", 1);

	// initialize
	gsl_matrix_set_zero(beta_t);
	gsl_matrix_set_zero(sderr_t);
	gsl_matrix_set_zero(t_t);
	gsl_matrix_set_zero(RSS_t);
	gsl_matrix_set_all(pv_t, 2);	// default p-value is bigger than 1

	if(calculate)
	{
		// I = B'B
		gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1, B_t, 0, I);

		if(gsl_linalg_cholesky_decomp(I) == GSL_EDOM)
		{
			fprintf(stderr, "Cholesky decomposition failed on B'B.\n");
			return OLS_COLINEAR;
		}

		// I = (B'B)^-1
		gsl_linalg_cholesky_invert(I);

		// D = I B'
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, I, B_t, 0, D);

		//gamma_2 = D R
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, D, R_t, 0, gamma_2);
	}

	// g = R - B gamma_2
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, gamma_2, B_t, 0, g_t);

	j = g_t->size1 * g_t->tda;
	for(i=0;i< j; i++) g_t->data[i] = R_t->data[i] - g_t->data[i];

	// |g|^2 = g'g
	for(i=0;i<m;i++)
	{
		gsl_vector_view v = gsl_matrix_row(g_t, i);
		gsl_blas_ddot(&v.vector, &v.vector, &t);
		gsl_vector_set(g_nrm2, i, t);

		if(fabs(t) < EPS && verbose_flag)
		{
			fprintf(stderr, "Warning in FWL: column %lu of X is collinear with background.\n", i);
		}
	}

	//gamma_1 = D Y
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, D, Y_t, 0, gamma_1);

	// f = Y - B gamma_1
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, gamma_1, B_t, 0, f_t);

	j = f_t->size1 * f_t->tda;
	for(i=0; i<j; i++) f_t->data[i] = Y_t->data[i] - f_t->data[i];

	// |f|^2 = f'f
	for(i=0;i<k;i++)
	{
		gsl_vector_view v = gsl_matrix_row(f_t, i);
		gsl_blas_ddot(&v.vector, &v.vector, &t);
		gsl_vector_set(f_nrm2, i, t);
	}

	// f'g
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, f_t, g_t, 0, beta_t);

	// beta = f'g / g'g
	for(i=0;i<k;i++)
	{
		gsl_vector_view v = gsl_matrix_row(beta_t, i);

		for(j=0;j<m;j++)
		{
			denom = gsl_vector_get(g_nrm2, j);

			// g is residual of R~B. If two small, means collinearity exists
			if(fabs(denom) < EPS)
			{
				colinear ++;
				t = DBL_MAX;
			}else{
				t = gsl_vector_get(&v.vector, j)/denom;
			}

			gsl_vector_set(&v.vector, j, t);
		}
	}

	// Univariate FWL shortcut: sderr^2 = (|f|^2/|g|^2 - beta^2 )/degree

	for(i=0;i<k;i++)
	{
		for(j=0;j<m;j++)
		{
			t = gsl_matrix_get(beta_t, i, j);
			if (t==DBL_MAX) continue;

			denom = gsl_vector_get(g_nrm2,j);

			// RSS first
			se = gsl_vector_get(f_nrm2,i)-t*t*denom;
			gsl_matrix_set(RSS_t, i,j, se);

			// OLS t-test
			if(t_test)
			{
				se /= (degree*denom);

				if(fabs(se)<EPS)
				{
					if(verbose_flag)
						fprintf(stderr, "Warning: standard error %e smaller than %e\n", se, EPS);

					se = EPS;
				}else{
					se = sqrt(se);
				}

				gsl_matrix_set(sderr_t, i,j, se);

				t /= se;
				gsl_matrix_set(t_t, i,j, t);

				t = 2*gsl_cdf_tdist_Q(fabs(t), degree);
				gsl_matrix_set(pv_t, i,j, t);
			}
		}
	}

	return colinear;
}



int FWL_vector( FWL_space *space,
		const int calculate,
		const int t_test,

		const gsl_matrix *B_t,
		const gsl_matrix *R_t,
		const gsl_vector *Y,

		gsl_vector *beta,
		gsl_vector *sderr,
		gsl_vector *t,
		gsl_vector *pv,
		gsl_vector *RSS)
{
	// make vectors to matrix views and call matrix version FWL

	size_t m = beta->size;

	check_vector_dimension(beta, m, "beta", 1);
	check_vector_dimension(sderr, m, "sderr", 1);
	check_vector_dimension(t, m, "t", 1);
	check_vector_dimension(pv, m, "pv", 1);
	check_vector_dimension(RSS, m, "RSS", 1);

	gsl_matrix_const_view
		Y_t = gsl_matrix_const_view_vector(Y, 1, Y->size);

	gsl_matrix_view
		beta_t = gsl_matrix_view_vector(beta, 1, m),
		sderr_t = gsl_matrix_view_vector(sderr, 1, m),
		t_t = gsl_matrix_view_vector(t, 1, m),
		pv_t = gsl_matrix_view_vector(pv, 1, m),
		RSS_t = gsl_matrix_view_vector(RSS, 1, m);

	return FWL_matrix(space, calculate, t_test, B_t, R_t, &Y_t.matrix,
			&beta_t.matrix, &sderr_t.matrix, &t_t.matrix, &pv_t.matrix, &RSS_t.matrix);
}

void Forward_space_resize(Forward_space *r, const size_t p, const size_t n)
{
	resize_matrix(r->X_t, p, n, NULL);
	resize_vector(r->Y, n, NULL);
	resize_vector(r->product, p, NULL);
}

Forward_space * Forward_space_alloc(size_t p, size_t n, gsl_matrix *X_t)
{
	Forward_space *r = (Forward_space*)malloc(sizeof(Forward_space));

	r->X_t = X_t;

	r->Y = gsl_vector_alloc(n);
	r->product = gsl_vector_alloc(p);

	// U_t is a pointer matrix to unselected elements
	r->U_t = (gsl_matrix*)malloc(sizeof(gsl_matrix));
	r->U_t->block = (gsl_block*)malloc(sizeof(gsl_block));

	return r;
}

void Forward_space_free(Forward_space *r)
{
	gsl_vector_free(r->Y);
	gsl_vector_free(r->product);

	free(r->U_t->block);
	free(r->U_t);

	free(r);
}


void forward_selection_QR(Forward_space *space, double *Y_last)
{
	gsl_matrix *X_t = space->X_t;
	size_t p = X_t->size1, n = X_t->size2;

	if(p>=n) return;

	QR_condense(X_t->data, space->Y->data, n, p, Y_last);

	resize_matrix(X_t, p, p, NULL);
	resize_vector(space->Y, p, NULL);
}


/*
 * Forward selection problem of Y ~ X. X dimension is n*p.
 * When p <= n, there is a more efficient implementation based on QR decomposition.
 * However, we assume p > n for some cases and use naive implementation here.
 * */

int forward_selection(Forward_space *space, size_t index_array[], double metric_array[], size_t stop_cnt)
{
	double beta, x2, beta_max, x2_max, Y_last = 0;

	// speed up by reducing n*p to p*p
	forward_selection_QR(space, &Y_last);
	Y_last = Y_last * Y_last;

	//CAUTION: After QR decomposition, p is unchanged. But n will reduce to p if n>p.

	gsl_matrix *U_t = space->U_t, *X_t = space->X_t;
	gsl_vector *Y = space->Y, *product = space->product;

	size_t i, j, inx_max, selected, unselected, collinear_jumped = 0, p = X_t->size1, n = X_t->size2;

	check_vector_dimension(Y, n, "Y", 1);

	resize_matrix(U_t, p, n, X_t->data);
	resize_vector(product, p, NULL);

	selected = 0;
	unselected = p;

	if (p < stop_cnt) stop_cnt = p;

	while(unselected > 0 && selected < stop_cnt)
	{
		// RSS = \sum(Yi-Xi*beta)^2 = \sum Yi^2 - (X'Y)^2/(\sum Xi^2)
		// minimum of RSS is maximum of (X'Y)^2/(\sum Xi^2)

		gsl_blas_dgemv(CblasNoTrans, 1, U_t, Y, 0, product);

		beta_max = 0;
		inx_max = 0;

		for(j=i=0 ; i<unselected ; i++)
		{
			gsl_vector_view xi = gsl_matrix_row(U_t, i);
			gsl_blas_ddot(&xi.vector, &xi.vector, &x2);

			// collinear with previous selected variables or start with zeros
			if(x2 < EPS) continue;

			beta = gsl_vector_get(product, i);
			beta = beta * beta/x2;

			// lag behind by jumped vectors
			if(j<i)
			{
				gsl_matrix_set_row(U_t, j, &xi.vector);
				index_array[j] = index_array[i];
			}

			if(beta > beta_max)
			{
				beta_max = beta;
				x2_max = x2;
				inx_max = j;
			}

			j ++ ;
		}

		//fprintf(stdout, "max index: %lu beta = %e\n", inx_max, beta_max);

		collinear_jumped += unselected - j;

		unselected = j;

		// nothing can be selected
		if(unselected == 0) break;


		gsl_blas_ddot(Y,Y,&beta);

		// adjust Y last element if QR if applied
		metric_array[selected] = beta + Y_last - beta_max;

		// swap selected variable to the beginning
		if(inx_max != 0)
		{
			gsl_matrix_swap_rows(U_t, 0, inx_max);

			i = index_array[0];
			index_array[0] = index_array[inx_max];
			index_array[inx_max] = i;
		}

		// Orthogonalize the rest to selected variable x0
		gsl_vector_view x0 = gsl_matrix_row(U_t, 0);

		// Y = Y - X0*beta
		gsl_blas_ddot(&x0.vector, Y, &beta);
		gsl_blas_daxpy(-beta/x2_max, &x0.vector, Y);


		// step forward
		selected++;
		index_array++;
		unselected--;

		resize_matrix(U_t, unselected, U_t->size2, U_t->data + n);
		resize_vector(product, unselected, NULL);

		// Orthogonalize X
		gsl_blas_dgemv(CblasNoTrans, 1, U_t, &x0.vector, 0, product);

		for(i=0;i<unselected;i++)
		{
			gsl_vector_view xi = gsl_matrix_row(U_t, i);

			beta = gsl_vector_get(product, i)/x2_max;

			// Xi = Xi - X0 * beta , then standardize by RSS
			gsl_blas_daxpy(-beta, &x0.vector, &xi.vector);
		}
	}

	if(collinear_jumped && verbose_flag)
	{
		fprintf(stderr, "Warning: %lu jumped by collinearity.\n", collinear_jumped);
	}

	// set X_t matrix size to selected
	resize_matrix(X_t, selected, X_t->size2, NULL);

	return OLS_SUCCESS;
}
