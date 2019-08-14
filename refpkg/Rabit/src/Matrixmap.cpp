#include "Matrixmap.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
using namespace std;


double string_to_number(const string &s, const string &errmsg)
{
	double r;
	char *endstr;

	r = strtod(s.c_str(), &endstr);

	if (*endstr!='\0' || endstr == s.c_str())
	{
		cerr << "Error: " << errmsg << ". \"" << s << "\" is not a number." << endl;
		exit(1);
	}

	return r;
}


Matrix_map::Matrix_map()
:Nrow(0), Ncol(0) {}


Matrix_map::~Matrix_map()
{
	for (vector<double*>::iterator iter = mat.begin(); iter!=mat.begin(); iter++) delete[] *iter;
}


void Matrix_map::read(const string &file, const size_t append_count)
{
	// temporary variables for parsing
	size_t i, j;
	string line;
	double *arr;

	vector<string> parse_vec;
	map<string, size_t>::iterator miter;

	// read in matrix to temporary space
	ifstream fin(file.c_str());

	if(fin.fail())
	{
		cerr << "Error: Cannot open \"" << file << "\"." << endl;
		exit(1);
	}

	getline(fin, line, '\n');

	istringstream iss(line);

	// load column names
	for(Ncol=0; getline(iss,line,'\t'); Ncol++)
	{
		if(col_map.find(line) != col_map.end())
		{
			cerr << "Error: Duplicated column ID \"" << line << "\"." << endl;
			exit(1);
		}

		colnames.push_back(line);
		col_map[line] = Ncol;
	}

	iss.clear();

	if(Ncol==0)
	{
		cerr << "Error: No columns in matrix \"" << file << "\"." << endl;
		exit(1);
	}

	// adjust column numbers for extra space
	for(i=0; i<append_count; i++, Ncol++)
	{
		colnames.push_back("Intercept");
	}

	// read in numbers
	for(i=0; getline(fin,line,'\n'); i++, Nrow++)
	{
		iss.str(line);
		getline(iss,line,'\t');

		if(iss.fail())
		{
			cerr << "Error: Empty line " << i << " in \"" << file << "\"." << endl;
			exit(1);
		}

		if(row_map.find(line) != row_map.end())
		{
			cerr << "Error: Duplicated row ID \"" << line << "\"." << endl;
			exit(1);
		}

		rownames.push_back(line);
		row_map[line] = i;

		arr = new double[Ncol];
		mat.push_back(arr);

		while(getline(iss, line, '\t')) parse_vec.push_back(line);
		iss.clear();

		for(j=0;j<append_count;j++) parse_vec.push_back("1");

		if(Ncol != parse_vec.size())
		{
			cerr << "Error: row ID \"" << rownames[i] << "\" has different number of elements " << parse_vec.size()
				<< " than column number " << Ncol << '.' << endl;

			exit(1);
		}

		// convert line strings to numbers
		for (j=0;j<Ncol;j++) arr[j] = string_to_number(parse_vec[j], "row ID \"" + rownames[i] + "\"");

		parse_vec.clear();
	}

	fin.close();
}

void Matrix_map::print() const
{
	size_t i, j;
	map<string, size_t>::const_iterator miter;

	double *arr;

	for (i=0;i<Ncol;i++) cout << '\t' << colnames[i];
	cout << '\n';

	for(i=0;i<Nrow;i++)
	{
		cout << rownames[i];

		miter = row_map.find(rownames[i]);
		arr = mat[miter->second];

		for (j=0; j < Ncol; j++)
		{
			cout << '\t' << arr[j];
		}

		cout << '\n';
	}

	cout << endl;
}
