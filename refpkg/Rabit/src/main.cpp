#include "Matrixmap.h"
#include "Rabit.h"
#include <cstring>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;


// verbose output
int verbose_flag = 0;


void get_common_matrix_name(
		const vector<Matrix_map*> &matrix_vec,
		vector<string> &common_names, const bool row_flag)
{
	// align the row names of all matrices
	vector<string> A, B;
	vector<Matrix_map*>::const_iterator miter;

	for(miter=matrix_vec.begin(); miter!=matrix_vec.end(); miter++)
	{
		if(row_flag)
			B = (*miter)->rownames;
		else
			B = (*miter)->colnames;

		sort(B.begin(), B.end());

		if(miter == matrix_vec.begin())
			common_names = B;
		else
			intersect(A, B, common_names);

		A = common_names;
	}
}


gsl_matrix *matrix_map_to_gsl(
		const vector<string> &rownames,
		const vector<string> &colnames,
		const Matrix_map *mat,
		const bool transpose)
{
	size_t i,j,
		Nrow = rownames.size(), Ncol = colnames.size(),
		*arrinx = new size_t[Ncol];

	double *arr;
	gsl_matrix *result;

	if(transpose)
		result = gsl_matrix_alloc(Ncol, Nrow);
	else
		result = gsl_matrix_alloc(Nrow, Ncol);

	map<string, size_t>::const_iterator miter;

	for (i=0;i<Ncol;i++)
	{
		miter = mat->col_map.find(colnames[i]);

		if(miter==mat->col_map.end())
		{	// I assume this never happens
			cerr << "Fatal error: missing element \"" << colnames[i] << "\"." << endl;
			exit(1);
		}

		arrinx[i] = miter->second;
	}


	for(i=0; i<Nrow; i++)
	{
		miter = mat->row_map.find(rownames[i]);

		if(miter==mat->row_map.end())
		{	// I assume this never happens
			cerr << "Fatal error: missing element \"" << rownames[i] << "\"." << endl;
			exit(1);
		}

		arr = mat->mat[miter->second];

		for(j=0;j<Ncol;j++)
		{
			if(transpose)
				gsl_matrix_set(result, j, i, arr[arrinx[j]]);
			else
				gsl_matrix_set(result, i, j, arr[arrinx[j]]);
		}
	}

	return result;
}



void write_gsl_matrix(const gsl_matrix *m,
		const vector<string> &rownames, const vector<string> &colnames,
		ofstream &fout, const bool transpose)
{
	double v;
	size_t i,j, Nrow = rownames.size(), Ncol = colnames.size();

	for(i=0;i<Ncol;i++) fout << colnames[i] << (i==Ncol-1?'\n':'\t');

	for(i=0;i<Nrow;i++)
	{
		fout << rownames[i];

		for(j=0;j<Ncol;j++)
		{
			if (transpose)
				v = gsl_matrix_get(m,j,i);
			else
				v = gsl_matrix_get(m,i,j);

			fout << '\t' << v;
		}

		fout << '\n';
	}
}


void normal_transform_gsl_matrix(gsl_matrix *m)
{
	size_t i, j, n = m->size1, p = m->size2;

	double *arr = new double[p];

	for(i=0;i<n;i++)
	{
		gsl_vector_view r = gsl_matrix_row(m,i);

		for(j=0;j<p;j++) arr[j] = gsl_vector_get(&r.vector, j);

		normal_transform(arr, arr, p);

		for(j=0;j<p;j++) gsl_vector_set(&r.vector, j, arr[j]);
	}

	delete[] arr;
}



int main(int argc, char *argv[])
{
	gsl_matrix *X, *Y, *B, *C;
	size_t i, j, parseCnt = (argc-1)/2, n, p, k, m;

	string
		Xfile,	// X matrix file in regression
		Yfile,	// Y response matrix file
		Bfile,	// Background confounding factor matrix for X
		output, value, type, title;

	// Background confounding factor matrices for Y
	set<string> Cfile_vec;
	set<string>::iterator Cfile_vec_iter;

	double FDR_thres = 0.05, *arr, exit_point =0.1;

	// transform Y vectors to normal distribution
	bool normal_transform = false, select_best = false, run_forward = true, joint_output = true;

	////////////////////////////////////////////////////////////////////////////////////////
	// Part 0: parameter input and check
	if ( argc < 7 )
	{
		if(argc == 2 && (value = argv[1]) == "-help")
		{
			cout << "\nRABIT: Regression Analysis with Background InTegration\n" << endl;
			cout << "Usage: Rabit -x X -y Y -o output [OPTIONS]\n"<<endl;

			cout << "\tX: Matrix of variables in selection "<< endl;
			cout << "\tY: Matrix of response vectors "<< endl;

			cout << "\nOptions:" << endl;
			cout << "\t-b\tBackground factors of X. Default: empty "<< endl;
			cout << "\t-c\tBackground factors of Y. Default: empty "<< endl;
			cout << "\t-f\tFDR threshold. Range: (0,1]. Default: " << FDR_thres << '.' << endl;
			cout << "\t-t\tTransform Y to normal distribution.\n\t\tOptions: 1 (yes), 0 (no). Default: " << DISPLAY_BOOL(normal_transform) << '.' << endl;
			cout << "\t-e\tExit point if Y is almost predicted. Range: [0,1). Default: " << exit_point << '.' << endl;

			cout << "\t-s\tSelect one best variable in X for each category.\n\t\tOptions: 1 (yes) or 0 (no). Default: " << DISPLAY_BOOL(select_best) << '.' << endl;
			cout << "\t-r\tRun forward stepwise selection.\n\t\tOptions: 1 (yes) or 0 (no). Default: " << DISPLAY_BOOL(run_forward) << '.' << endl;
			cout << "\t-j\tJoint output.\n\t\tOptions: 1 (yes) or 0 (no). Default: " << DISPLAY_BOOL(joint_output) << '.' << endl;
			cout << "\t-v\tVerbose output such as warnings.\n\t\tOptions: 1 (yes) or 0 (no). Default: " << DISPLAY_BOOL(verbose_flag) << '.' << endl;

			cout << "\nReport bugs to Peng Jiang (peng.jiang.software@gmail.com)\n" <<endl;
			exit(0);
		}else{
			cerr << "Insufficient number of arguments, do \"Rabit -help\" for help."<<endl;
			exit(1);
		}
	}

	// read in all parameters
	for(i=0;i<parseCnt;i++)
	{
		type = argv[2*i+1];
		value = argv[2*i+2];

		if(type == "-x"){
			Xfile = value;

		}else if (type == "-y"){
			Yfile = value;

		}else if (type == "-b"){
			Bfile = value;

		}else if (type == "-c")
		{
			if(Cfile_vec.find(value) != Cfile_vec.end())
			{
				if (verbose_flag)
					cerr << "Warning: Ignore duplicated background file \"" << value << "\"." << endl;
			}else
				Cfile_vec.insert(value);

		}else if (type == "-f"){
			FDR_thres = string_to_number(value, "bad parameter to \"-f\"");

			if(FDR_thres <= 0 || FDR_thres > 1)
			{
				cerr << "Error: FDR threshold \"" << FDR_thres << "\" is out of range (0,1]." << endl;
				exit(1);
			}

		}else if (type == "-e"){
			exit_point = string_to_number(value, "bad parameter to \"-e\"");

			if(exit_point < 0 || exit_point >= 1)
			{
				cerr << "Error: exit_point \"" << exit_point << "\" is out of range [0,1)." << endl;
				exit(1);
			}

		}else if (type == "-t")
		{
			if(value != "0" && value != "1")
			{
				cerr << "Error: Invalid parameter input \"" << value << "\" for \"-t\"." << endl;
				cerr << "Please only use 1 (yes) or 0 (no)." << endl;
				exit(1);
			}

			normal_transform = (value[0] == '1');

		}else if (type == "-s")
		{
			if(value != "0" && value != "1")
			{
				cerr << "Error: Invalid parameter input \"" << value << "\" for \"-s\"." << endl;
				cerr << "Please only use 1 (yes) or 0 (no)." << endl;
				exit(1);
			}

			select_best = (value[0] == '1');

		}else if (type == "-r")
		{
			if(value != "0" && value != "1")
			{
				cerr << "Error: Invalid parameter input \"" << value << "\" for \"-r\"." << endl;
				cerr << "Please only use 1 (yes) or 0 (no)." << endl;
				exit(1);
			}

			run_forward = (value[0] == '1');

		}else if (type == "-j")
		{
			if(value != "0" && value != "1")
			{
				cerr << "Error: Invalid parameter input \"" << value << "\" for \"-j\"." << endl;
				cerr << "Please only use 1 (yes) or 0 (no)." << endl;
				exit(1);
			}

			joint_output = (value[0] == '1');

		}else if (type == "-v")
		{
			if(value != "0" && value != "1")
			{
				cerr << "Error: Invalid parameter input \"" << value << "\" for \"-v\"." << endl;
				cerr << "Please only use 1 (yes) or 0 (no)." << endl;
				exit(1);
			}

			verbose_flag = (value[0] == '1'?1:0);

		}else if (type == "-o"){
			output = value;

		}else if (type == "-help"){
			cerr << "Please don't use \"-help\" as parameter input." << endl;
			exit(1);

		}else{
			cerr << "Cannot recognize parameter \""<< type << "\"." << endl;
			exit(1);
		}
	}

	// required fields

	if(Xfile.empty())
	{
		cerr << "Error: Cannot find X file." << endl;
		exit(1);
	}

	if(Yfile.empty())
	{
		cerr << "Error: Cannot find Y file." << endl;
		exit(1);
	}

	if(output.empty())
	{
		cerr << "Error: Cannot find output file." << endl;
		exit(1);
	}

	vector<string>
		rownames,	// common row names across all matrices
		colnames		// common column names across Y and all its background matrices
		;

	// read in matrix maps from file
	Matrix_map *X_map = new Matrix_map, *Y_map = new Matrix_map, *B_map = new Matrix_map, *C_map;

	if (verbose_flag) cout << "Load data" << endl;

	X_map->read(Xfile);
	Y_map->read(Yfile);

	// extra space size for Y background matrix plus 1 for Intercept
	size_t append_size = 1 + Cfile_vec.size();


	// construct X background matrix B

	if(Bfile.empty())
	{
		// no X background, construct all 1 background as Intercept
		map<string, size_t>::iterator map_iter;

		// same row dimension with X
		B_map->Nrow = X_map->Nrow;
		B_map->row_map = X_map->row_map;
		B_map->rownames = X_map->rownames;

		// repeated column
		B_map->Ncol = append_size;
		for(i=0;i<append_size;i++) B_map->colnames.push_back("Intercept");

		for(i=0; i<B_map->Nrow; i++)
		{
			arr = new double[append_size];
			for(j=0;j<append_size;j++) arr[j] = 1;
			B_map->mat.push_back(arr);
		}

	}else{
		B_map->read(Bfile, append_size);
	}

	// fill in B map column names
	B_map->col_map["Intercept"] = B_map->Ncol-1;


	for(i=B_map->Ncol-append_size, Cfile_vec_iter = Cfile_vec.begin();
			Cfile_vec_iter!=Cfile_vec.end();	Cfile_vec_iter++, i++)
	{
		value = *Cfile_vec_iter;

		// search for title file separator
		j = value.find(':');
		if(j!=string::npos) value = value.substr(0, j);

		B_map->colnames[i] = value;
		B_map->col_map[value] = i;
	}

	// align all matrices for common names

	vector<Matrix_map*>
		matrix_vec,		// All matrices together for common row names
		matrix_vec_Y;	// Y and Y background matrices for common column names

	matrix_vec.push_back(X_map);
	matrix_vec.push_back(B_map);
	matrix_vec.push_back(Y_map);

	for(Cfile_vec_iter = Cfile_vec.begin(); Cfile_vec_iter!=Cfile_vec.end();	 Cfile_vec_iter++)
	{
		value = *Cfile_vec_iter;

		// search for title file separator
		i = value.find(':');
		if(i!=string::npos) value = value.substr(i+1);

		C_map = new Matrix_map;
		C_map->read(value);

		matrix_vec.push_back(C_map);
		matrix_vec_Y.push_back(C_map);
	}

	// push Y at last
	matrix_vec_Y.push_back(Y_map);

	get_common_matrix_name(matrix_vec, rownames, true);

	if(rownames.empty())
	{
		cerr << "Error: No common row names exist across all matrices." << endl;
		exit(1);
	}

	get_common_matrix_name(matrix_vec_Y, colnames, false);

	if(colnames.empty())
	{
		cerr << "Error: No common column names exist across all matrices." << endl;
		exit(1);
	}

	// converting all matrix maps to gsl matrices

	vector<string> colnames_X = X_map->colnames, colnames_B = B_map->colnames;

	// convert all matrices to transposed for row major order of each vector

	X = matrix_map_to_gsl(rownames, colnames_X, X_map, true);
	delete X_map;

	B = matrix_map_to_gsl(rownames, colnames_B, B_map, true);
	delete B_map;

	vector<gsl_matrix*> C_vec;

	for(vector<Matrix_map *>::iterator viter = matrix_vec_Y.begin(); viter != matrix_vec_Y.end(); viter++)
	{
		C_map = *viter;

		C = matrix_map_to_gsl(rownames, colnames, C_map, true);
		delete C_map;

		C_vec.push_back(C);
	}

	// pop Y from last
	Y = C_vec.back();
	C_vec.pop_back();

	// transform Y vectors to normal distribution to improve t-test accuracy
	if(normal_transform) normal_transform_gsl_matrix(Y);


	n = rownames.size();
	p = colnames_B.size();
	k = colnames.size();
	m = colnames_X.size();

	// start using RABIT here
	Rabit *rabit = new Rabit(n,p,k,m);

	if(verbose_flag)
	{
		cout << "Regression units:\t" << n << endl;
		cout << "Background factors:\t" << p-1 << endl;
		cout << "Candidate features:\t" << m << endl;
		cout << "Response variables:\t" << k << endl;
	}

	ofstream
		fout_t((output + ".t").c_str()),
		fout_FDR((output + ".FDR").c_str()),
		fout(output.c_str()),
		fout_selected;

	if(select_best) fout_selected.open((output + ".selected").c_str());

	// algorithm start
	if (verbose_flag) cout << "Start algorithm" << endl;

	Rabit_result *rptr;
	vector<Rabit_result*> results;
	vector<Rabit_result*>::iterator result_iter;

	rabit->run(X, Y, B, C_vec, run_forward, select_best, FDR_thres, exit_point, colnames_X, results);

	for(i=0;i<m;i++)
	{
		fout_t << colnames_X[i] << (i==m-1?'\n':'\t');
		fout_FDR << colnames_X[i] << (i==m-1?'\n':'\t');
	}

	for(result_iter = results.begin(); result_iter != results.end(); result_iter++)
	{
		rptr = *result_iter;

		rptr->print(fout, fout_t, fout_FDR, colnames[rptr->index], colnames_B, colnames_X);

		if (select_best)
			rptr->print_selection(fout_selected, colnames[rptr->index]);

		delete rptr;
	}

	fout.close();
	fout_t.close();
	fout_FDR.close();
	if(select_best) fout_selected.close();

	gsl_matrix_free(X);
	gsl_matrix_free(Y);
	gsl_matrix_free(B);

	for(vector<gsl_matrix*>::iterator C_vec_iter = C_vec.begin(); C_vec_iter != C_vec.end(); C_vec_iter++) delete *C_vec_iter;

	delete rabit;
	return 0;
}
