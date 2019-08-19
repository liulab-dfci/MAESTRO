#ifndef UTIL_H_
#define UTIL_H_

#include <gsl/gsl_cdf.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

#define	EPS	1e-20
#define	DISPLAY_BOOL(x) ((x)?("1 (yes)"):("0 (no)"))


// intersection of two sorted vectors. Assume vector is ordered.
template <class T>
void intersect(const vector<T> &s1, const vector<T> &s2, vector<T> &result)
{
	result.clear();

	typename vector<T>::const_iterator
		first1 = s1.begin(), last1 = s1.end(),
		first2 = s2.begin(), last2 = s2.end();

	while (first1 != last1 && first2 != last2)
	{
		if (*first1 < *first2)
			first1++;
		else if (*first2 < *first1)
			first2++;
		else {
			result.push_back(*first1);
			first1++;
			first2++;
		}
	}
}


// print vector elements
template <class T>
void print_vector(const vector<T> &v, const string &title="")
{
	typename vector<T>::const_iterator iter;

	cout << title;
	for(iter=v.begin(); iter!=v.end(); iter++) cout << '\t' << *iter;
	cout << endl;
}


// rank node structure for rank transform function
template <class T>
class rank_node
{
public:
	rank_node(T v, size_t i):v(v),i(i){}
	bool operator < (const rank_node &r) const {return v < r.v;}

	T v;
	size_t i;
};

// transform values in src to quantiles in dst. If dst == src, transform in place
template <class T>
void quantile_transform(double dst[], T src[], const size_t N)
{
	size_t i;
	vector<rank_node<T> > sortvec;

	for(i=0;i<N;i++) sortvec.push_back(rank_node<T>(src[i],i));

	sort(sortvec.begin(), sortvec.end());

	for(i=0;i<N;i++) dst[sortvec[i].i] = (double)(i+1)/(N+1);
}


// transform values in src to quantiles in dst. If dst == src, transform in place
template <class T>
void normal_transform(double dst[], T src[], const size_t N)
{
	size_t i;
	double aver=0, sd=0, v;

	// no point to do normal transformation
	if(N==1) return;

	for (i=0;i<N;i++)
	{
		v = (double)src[i];
		aver += v;
		sd += v*v;
	}

	aver /= N;
	sd = (sd - N*aver*aver)/(N-1);
	sd = (fabs(sd)<EPS?0:sqrt(sd));

	quantile_transform(dst, src, N);

	for(i=0;i<N;i++) dst[i] = gsl_cdf_gaussian_Pinv(dst[i], sd) + aver;
}

size_t triangle_index(size_t i, size_t j);

void Benjamini_Hochberg(double FDR[], double pvalue[], const size_t N);

#endif /* UTIL_H_ */
