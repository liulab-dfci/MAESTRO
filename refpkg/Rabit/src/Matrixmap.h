#ifndef MATRIXMAP_H_
#define MATRIXMAP_H_

#include <string>
#include <vector>
#include <map>
using namespace std;

double string_to_number(const string &s, const string &errmsg="");

class Matrix_map {
public:
	Matrix_map();
	~Matrix_map();

	map<string, size_t> row_map, col_map;
	vector<string> rownames, colnames;

	size_t Nrow, Ncol;

	vector<double*> mat;

	// read matrix from file, append 1 to extra space
	void read(const string &file, const size_t append_count = 0);

	// print matrix
	void print() const;
};

#endif /* MATRIXMAP_H_ */
