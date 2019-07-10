/*
 * kfunc.h
 *
 *  Created on: May 1, 2015
 *      Author: nek3d
 */

#ifndef KFUNC_H_
#define KFUNC_H_

#include <math.h>
#include <stdlib.h>


long double _lbinom(long long n, long long k);
long double _hypergeo(long long n11, long long n1_, long long n_1, long long n);

typedef struct {
    long long n11, n1_, n_1, n;
    long double p;
} _hgacc_t;

// incremental version of hypergenometric distribution
long double _hypergeo_acc(long long n11, long long n1_, long long n_1, long long n, _hgacc_t *aux);
long double _kt_fisher_exact(long long n11, long long n12, long long n21, long long n22, long double *_left, long double *_right, long double *two);


#endif /* KFUNC_H_ */
