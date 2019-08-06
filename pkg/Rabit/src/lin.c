#include "lin.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

// external link to lapack
extern void dgeqrf_(int *M, int *N, double *A, int *LDA, double *TAU, double *WORK, int *LWORK, int *INFO);


void print_vector_product(const double A[], const double B[], const int n)
{
	double r=0;
	int i;

	for(i=0;i<n;i++) r += A[i] * B[i];

	fprintf(stdout, "%e\n", r);
}


void QR_condense(double X[], double Y[], int n, int p, double *Y_last)
{
	int i, lwork = -1, info = 0;
	double *dst, *src, *work = NULL, work_query;

	if(p >= n)
	{
		fprintf(stderr, "Please don't use QR decomposition for p>=n\n");
		exit(1);
	}

	// append Y at the last column of X
	memcpy(X + n*p, Y, n*sizeof(double));

	// QR decomposition
	i = p+1;

	// determine best work size first
	dgeqrf_(&n, &i, X, &n, Y, &work_query, &lwork, &info);

	if(info != 0)
	{
		fprintf(stderr, "Error: lapack DGEQRF memory error.\n");
		exit(1);
	}

	// allocated work space
	lwork = (int)work_query;
	work = (double*)malloc(lwork*sizeof(double));

	// start working
	dgeqrf_(&n, &i, X, &n, Y, work, &lwork, &info);
	free(work);

	// compact space
	for(i=0, dst=X, src=X; i<p; i++, dst+=p, src+= n)
	{
		if(dst != src) memcpy(dst, src, (i+1)*sizeof(double));

		memset(dst+i+1, 0, (p-i-1)*sizeof(double));
	}

	memcpy(Y, src, p*sizeof(double));

	*Y_last = src[p];
}
