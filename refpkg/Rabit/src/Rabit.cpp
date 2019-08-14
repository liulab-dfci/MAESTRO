#include "Rabit.h"


Rabit::Rabit(const size_t n, const size_t p, const size_t k, const size_t m)
:p(p), n(n), m(m), k(k), count(0), sigma2(0)
{
	// keep 1 buffer for Y in QR decomposition
	X_t = gsl_matrix_alloc(m+1, n);

	forward_space = Forward_space_alloc(m + 1, n, X_t);
	fwl_space = FWL_space_alloc(p, n, m, 1, X_t);
	ols_space = OLS_space_alloc(min<size_t>(m+p, n-1), n, 1);

	// results vectors
	beta = gsl_vector_alloc(m + p);
	sderr = gsl_vector_alloc(m + p);
	t = gsl_vector_alloc(m + p);
	pv = gsl_vector_alloc(m + p);
	RSS = gsl_vector_alloc(m);

	// results right after FWL without forward selection
	t_FWL = gsl_vector_alloc(m);
	FDR_FWL = gsl_vector_alloc(m);

	// make sure the vector stride is one
	check_vector_dimension(pv, m+p, "pv", 1);

	// variable index in selection order
	index_array = new size_t[m];

	// metrics in selection process
	metric_array = new double[m];
	memset(metric_array, 0, m * sizeof(double));

	pass_array = new bool[m];
	memset(pass_array, 0, m * sizeof(bool));
}

Rabit::~Rabit()
{
	gsl_matrix_free(X_t);
	gsl_vector_free(beta);
	gsl_vector_free(sderr);
	gsl_vector_free(t);
	gsl_vector_free(pv);
	gsl_vector_free(RSS);
	gsl_vector_free(t_FWL);
	gsl_vector_free(FDR_FWL);

	delete[] index_array;
	delete[] metric_array;
	delete[] pass_array;

	FWL_space_free(fwl_space);
	Forward_space_free(forward_space);
	OLS_space_free(ols_space);
}


void Rabit::run(const gsl_matrix *X, const gsl_matrix *Y, gsl_matrix *B, const vector<gsl_matrix*> C_vec,
		const bool run_forward, const bool select_best, const double FDR_thres, const double exit_point,
		const vector<string> features, vector<Rabit_result*> &results)
{
	string value;
	gsl_matrix *C;
	double *arr;
	int status;

	set<string> included;

	vector<rank_node<double> > FDR_sortvec;
	vector<rank_node<double> >::iterator FDR_sortvec_iter;

	size_t i, j, calculate, inx,
		// for progress bar
		step = max<size_t>(k/100,1),

		// extra space in B for Y specific confounding factors
		append_size = 1 + C_vec.size();


	vector<gsl_matrix*>::const_iterator C_vec_iter;

	// forward selection space input matrix handle
	gsl_matrix *X_t = forward_space->X_t;


	for(i=0;i<k;i++)
	{
		// progress bar
		if(i%step==0){cout << 100*(i+1)/k << '%' << endl;}

		gsl_vector_const_view Yi = gsl_matrix_const_row(Y,i);

		// merge X Y background factors to B
		for(	C_vec_iter=C_vec.begin(), j=0; C_vec_iter != C_vec.end(); C_vec_iter++, j++)
		{
			C = *C_vec_iter;

			gsl_vector_const_view Ci = gsl_matrix_const_row(C, i);
			gsl_matrix_set_row(B, j + p-append_size, &Ci.vector);
		}

		// Frish Waugh Lovell selection
		// If only X background exists, many FWL variables will keep the same and be calculated once
		if(i && C_vec.empty())
			calculate = 0;
		else
			calculate = 1;


		beta->size = sderr->size = t->size = pv->size = m;

		resize_matrix(fwl_space->g_t, m, n, NULL);
		status = FWL_vector(fwl_space, calculate, 1, B, X, &Yi.vector, beta, sderr, t, pv, RSS);

		if(status == OLS_COLINEAR)
		{
			cerr << "Error: Background matrix is collinear." << endl;
			exit(1);
		}

		arr = pv->data;

		// convert P-values to FDRs
		Benjamini_Hochberg(arr, arr, m);

		/////////////////////////////////////////
		// Stepwise forward selection

		// recover forward space to maximum capacity
		Forward_space_resize(forward_space, m+1, n);	// 1 for copy Y is QR decomposition

		// copy down the results right after FWL, before forward selection
		gsl_vector_memcpy(t_FWL, t);

		// sort from small FDR to large FDR
		for(j=0;j<m;j++)
		{
			FDR_sortvec.push_back(rank_node<double>(arr[j], j));
			gsl_vector_set(FDR_FWL, j, arr[j]);
		}


		sort(FDR_sortvec.begin(), FDR_sortvec.end());


		// set start
		for(count = 0, FDR_sortvec_iter = FDR_sortvec.begin(), selection_map.clear();
				FDR_sortvec_iter != FDR_sortvec.end() ; FDR_sortvec_iter++)
		{
			inx = FDR_sortvec_iter->i;

			// filter previous included regulators
			if(select_best)
			{
				value = features[inx];
				value = value.substr(0, value.find('.'));

				if(included.find(value) != included.end()) continue;

				selection_map.push_back(pair<string,string>(value, features[inx]));
				included.insert(value);
			}

			if(FDR_sortvec_iter->v <= FDR_thres) pass_array[inx] = true;
		}

		for(j=0, count=0 ; j<m ; j++)
		{
			if(pass_array[j])
			{
				if(j>count)
				{
					gsl_vector_view v = gsl_matrix_row(X_t, j);
					gsl_matrix_set_row(X_t, count, &v.vector);
				}

				index_array[count] = j;
				count++;

				pass_array[j] = false;
			}
		}

		FDR_sortvec.clear();
		included.clear();

		// nothing to select from
		if(count == 0) continue;

		// compact X to filled size
		resize_matrix(X_t, count, n, NULL);

		// copy orthogonalized Y
		gsl_matrix_get_row(forward_space->Y, fwl_space->f_t, 0);

		if(run_forward)
		{
			forward_selection(forward_space, index_array, metric_array, n-p-1);

			count = X_t->size1;

			// nothing is selected
			if(count == 0) continue;


			// convert RSS to selection metric and select at minimum metric position
			count = convert_RSS_to_metric(count, exit_point) + 1;
		}


		// estimate selected linear model

		// copy everything to X
		resize_matrix(X_t, p + count, n, NULL);

		for(j=0 ; j<p ; j++)
		{
			gsl_vector_view v = gsl_matrix_row(B, j);
			gsl_matrix_set_row(X_t, j, &v.vector);
		}

		for(j=0 ; j<count ; j++)
		{
			inx = index_array[j];

			gsl_vector_const_view v = gsl_matrix_const_row(X, inx);
			gsl_matrix_set_row(X_t, p + j, &v.vector);
		}


		// resize OLS space for current regression
		beta->size = sderr->size = t->size = pv->size = p + count;
		OLS_space_resize(ols_space, p + count, n, 1);

		// FDR is p-value from now
		OLS_vector(ols_space, X_t, &Yi.vector, beta, sderr, t, pv, &sigma2);

		// output selected model
		results.push_back(clone_result(i));
	}
}



size_t Rabit::convert_RSS_to_metric(const size_t size, double exit_point)
{
	size_t i, inx_min = string::npos;
	double metric_min = DBL_MAX, t;

	if(p+size >= n)
	{
		cerr << "Error: Illegal number of features selected p + size >= n" << endl;
		exit(1);
	}

	sigma2 = metric_array[size-1]/(n-size-p);

	bool smallsigma2 = (fabs(sigma2) < EPS);

	//if(verbose_flag && smallsigma2)
	//	cerr << "Warning: Sigma2 is too small. Use different form for Mallow's Cp calculation." << endl;

	for(i=0 ; i<size ; i++)
	{
		// Mallow's Cp
		if(smallsigma2){
			t = (metric_array[i] + 2*(p+i+1)*sigma2)/n;

			if(i==0)
				exit_point = exit_point*t;
			else
				if(t < exit_point) break;

		}else
			t = metric_array[i]/sigma2 + 2*(p+i+1) - n;

		metric_array[i] = t;

		if(t < metric_min)
		{
			metric_min = t;
			inx_min = i;
		}
	}

	return inx_min;
}


Rabit_result *Rabit::clone_result(const size_t index)
{
	return new Rabit_result(beta, sderr, t, pv, t_FWL, FDR_FWL, index_array, metric_array, selection_map, p, count, index);
}


Rabit_result::Rabit_result(
	const gsl_vector *beta, const gsl_vector *sderr, const gsl_vector *t, const gsl_vector *pv,
	const gsl_vector *t_FWL, const gsl_vector *FDR_FWL,
	const size_t *index_array, const double *metric_array, const vector<pair<string, string> > &selection_map,
	const size_t p, const size_t count, const size_t index
):
	beta(gsl_vector_alloc(beta->size)),
	sderr(gsl_vector_alloc(sderr->size)),
	t(gsl_vector_alloc(t->size)),
	pv(gsl_vector_alloc(pv->size)),
	t_FWL(gsl_vector_alloc(t_FWL->size)),
	FDR_FWL(gsl_vector_alloc(FDR_FWL->size)),
	index_array(new size_t[count]),
	metric_array(new double[count]),
	selection_map(selection_map),
	p(p), count(count), m(t_FWL->size), index(index)
{
	gsl_vector_memcpy(this->beta, beta);
	gsl_vector_memcpy(this->sderr, sderr);
	gsl_vector_memcpy(this->t, t);
	gsl_vector_memcpy(this->pv, pv);
	gsl_vector_memcpy(this->t_FWL, t_FWL);
	gsl_vector_memcpy(this->FDR_FWL, FDR_FWL);

	memcpy(this->index_array, index_array, count * sizeof(size_t));
	memcpy(this->metric_array, metric_array, count * sizeof(double));
}

Rabit_result::~Rabit_result()
{
	gsl_vector_free(beta);
	gsl_vector_free(sderr);
	gsl_vector_free(t);
	gsl_vector_free(pv);
	gsl_vector_free(t_FWL);
	gsl_vector_free(FDR_FWL);

	delete[] index_array;
	delete[] metric_array;
}

void Rabit_result::print(ofstream &fout, ofstream &fout_t, ofstream &fout_FDR,
		const string title, const vector<string> backgrounds, const vector<string> features)
{
	size_t i, inx;

	fout << '>' << title << "\tCp\tEstimate\tStd. Error\tt-value\tPr(>|t|)\n";

	for(i=0 ; i<p ; i++)
	{
		fout << backgrounds[i] << '\t' << 0 << '\t'
			<< gsl_vector_get(beta, i) << '\t'
			<< gsl_vector_get(sderr, i) << '\t'
			<< gsl_vector_get(t, i) << '\t'
			<< gsl_vector_get(pv, i) << '\n';
	}

	for(i=0 ; i<count ; i++)
	{
		inx = index_array[i];

		fout << features[inx] << '\t' << metric_array[i] << '\t'
			<< gsl_vector_get(beta, p + i) << '\t'
			<< gsl_vector_get(sderr, p + i) << '\t'
			<< gsl_vector_get(t, p + i) << '\t'
			<< gsl_vector_get(pv, p + i) << '\n';
	}

	fout_t << title;
	fout_FDR << title;

	for(i=0 ; i<m ; i++)
	{
		fout_t << '\t' << gsl_vector_get(t_FWL, i);
		fout_FDR << '\t' << gsl_vector_get(FDR_FWL, i);
	}

	fout_t << '\n';
	fout_FDR << '\n';
}


void Rabit_result::print_selection(ofstream &fout, const string title)
{
	vector<pair<string, string> >::const_iterator iter;

	for(iter=selection_map.begin(); iter != selection_map.end(); iter++)
	{
		fout << title << '\t' << iter->first << '\t' << iter->second << '\n';
	}
}
