#include "util.h"

void Benjamini_Hochberg(double FDR[], double pvalue[], const size_t N)
{
	size_t i;
	double qvalue = 1;
	vector<rank_node<double> > sortvec;

	for(i=0;i<N;i++) sortvec.push_back(rank_node<double>(pvalue[i],i));

	sort(sortvec.begin(), sortvec.end());

	for(i=0;i<N;i++)
	{
		qvalue = min<double>(qvalue, sortvec[N-1-i].v*N/(N-i));
		FDR[sortvec[N-1-i].i] = qvalue;
	}
}

size_t triangle_index(size_t i, size_t j)
{
	if(i==j){
		cerr << "Error: cannot use i==j in triangle index." << endl;
		exit(1);

	}else if(i<j) swap<size_t>(i,j);

	return ((i*(i-1))>>1) + j;
}
