#ifndef GSL_UTIL_H_
#define GSL_UTIL_H_

#include <gsl/gsl_blas.h>

#define EPS 1e-20

// EPS correction for sqrt function on nearly 0 value
#define SQRT_EPS(x) (fabs(x)<EPS?0:sqrt(x))



// return 1 if v variation is 0 or size is smaller than 2
// if standardize, make norm2 = 1
int Z_normalize(gsl_vector *v, const int standardize);


// change the dimension of matrices and vectos
void resize_matrix(gsl_matrix *m, const size_t size1, const size_t size2, double *data);
void resize_vector(gsl_vector *v, const size_t size, double *data);

void print_vector(const gsl_vector *v, const char *title, FILE *fp);
void print_matrix(const gsl_matrix *X, const char *title, FILE *fp);


// check the dimension of matrices and vectos
void check_matrix_dimension(const gsl_matrix *X, const size_t n, const size_t p,
		const char *title, const int equal_stride);

void check_vector_dimension(const gsl_vector *X, const size_t n,
		const char *title, const size_t stride);

#endif /* GSL_UTIL_H_ */
