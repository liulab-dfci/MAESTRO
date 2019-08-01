#include "gsl_util.h"

#include <math.h>

void print_matrix(const gsl_matrix *X, const char *title, FILE *fp)
{
	size_t i,j;

	fprintf(fp,"Matrix %s:\n", title);

	for(i=0;i<X->size1;i++)
	{
		for(j=0;j<X->size2;j++)
		{
			fprintf(fp,"%f\t\t",gsl_matrix_get(X,i,j));
		}
		fprintf(fp,"\n");
	}
}

void print_vector(const gsl_vector *v, const char *title, FILE *fp)
{
	size_t i;

	fprintf(fp, "Vector %s:", title);

	for(i=0;i<v->size;i++)
	{
		fprintf(fp,"\t\t%f",gsl_vector_get(v,i));
	}
	fprintf(fp,"\n");
}


int Z_normalize(gsl_vector *v, const int standardize)
{
	double aver=0, sd=0, t;
	size_t i, n=v->size;

	if(n<=1) return 1;

	for(i=0 ; i<n ; i++)
	{
		t = gsl_vector_get(v, i);
		aver += t;
		sd += t*t;
	}

	aver /= n;
	sd -= n*aver*aver;

	// use n-1 for sample variance
	if(standardize == 0) sd /= n-1;

	if(fabs(sd) < EPS) return 1;

	sd = sqrt(sd);

	for(i=0 ; i<n ; i++)
	{
		t = gsl_vector_get(v, i);
		gsl_vector_set(v, i, (t-aver)/sd);
	}

	return 0;
}



void resize_matrix(gsl_matrix *m, const size_t size1, const size_t size2, double *data)
{
	m->size1 = size1;
	m->tda = m->size2 = size2;

	m->block->size = m->size1*m->size2;

	if(data!=NULL) m->block->data = m->data = data;
}

void resize_vector(gsl_vector *v, const size_t size, double *data)
{
	if(size>0) v->block->size = v->size = size;

	if(data!=NULL) v->data = v->block->data = data;
}


void check_matrix_dimension(const gsl_matrix *X, const size_t n, const size_t p,
		const char *title, const int equal_stride)
{
	if(X->size1 != n || X->size2 != p)
	{
		fprintf(stderr, "Matrix %s dimension is %lu * %lu but not expected %lu * %lu.\n", title, X->size1, X->size2, n, p);
		exit(1);
	}

	if(equal_stride != 0 && X->size2 != X->tda)
	{
		fprintf(stderr, "Matrix %s tda %lu is not equal to column %lu.\n", title, X->tda, X->size2);
		exit(1);
	}
}

void check_vector_dimension(const gsl_vector *X, const size_t n, const char *title, const size_t stride)
{
	if(X->size != n)
	{
		fprintf(stderr, "Vector %s length is %lu but not expected %lu.\n", title, X->size, n);
		exit(1);
	}

	if(X->stride != stride)
	{
		fprintf(stderr, "Vector %s stride %lu is not equal to %lu.\n", title, X->stride, stride);
		exit(1);
	}
}



