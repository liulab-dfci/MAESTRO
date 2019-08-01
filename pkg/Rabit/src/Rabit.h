#ifndef RABIT_H_
#define RABIT_H_

extern "C"{
#include "lm.h"
#include "gsl_util.h"
}

#include "util.h"
#include <fstream>
#include <cstring>
#include <set>
#include <vector>
#include <algorithm>
using namespace std;


// result node for Result analysis
class Rabit_result
{
public:
	Rabit_result(
		const gsl_vector *beta, const gsl_vector *sderr, const gsl_vector *t, const gsl_vector *pv,
		const gsl_vector *t_FWL, const gsl_vector *FDR_FWL,
		const size_t *index_array, const double *metric_array, const vector<pair<string, string> > &selection_map,
		const size_t p, const size_t count, const size_t index);

	~Rabit_result();

	// print result coefficients
	void print(
		ofstream &fout, ofstream &fout_t, ofstream &fout_FDR,
		const string title, const vector<string> backgrounds, const vector<string> features);

	// print selection map
	void print_selection(ofstream &fout, const string title);

	gsl_vector *beta, *sderr, *t, *pv, *t_FWL, *FDR_FWL;
	size_t *index_array;
	double *metric_array;

	vector<pair<string, string> > selection_map;

	// p: number of background factors
	// count: number of selected features
	// index: of current response
	size_t p, count, m, index;
};



class Rabit
{
public:
	Rabit(const size_t n, const size_t p, const size_t k, const size_t m);
	~Rabit();

	void run(
		const gsl_matrix *X,
		const gsl_matrix *Y,
		gsl_matrix *B,
		const vector<gsl_matrix*> C_vec,
		const bool run_forward,
		const bool select_best,
		const double FDR_thres,
		const double exit_point,
		const vector<string> features,

		// result output vector for each Y column
		vector<Rabit_result*> &results);

	// copy current result
	Rabit_result *clone_result(const size_t index);

private:
	// result vectors
	gsl_vector *beta, *sderr, *t, *pv, *RSS, *t_FWL, *FDR_FWL;

	// feature matrix for Orthogonalization and QR decomposition
	gsl_matrix *X_t;

	// dimensions
	size_t p, n, m, k, count;

	// index in selection
	size_t *index_array;

	// metric in selection
	double *metric_array;

	// variables after FWL selection
	bool *pass_array;

	double sigma2;

	vector<pair<string, string> > selection_map;

	// internal workspace

	// first step variable filtering
	FWL_space *fwl_space;

	// forward selection workspace
	Forward_space *forward_space;

	// estimate final model based on forward selection variables
	OLS_space *ols_space;

	size_t convert_RSS_to_metric(const size_t size, double exit_point);

};



#endif /* RABIT_H_ */
