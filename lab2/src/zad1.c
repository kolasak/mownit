#include "matrix.h"


int main(void) {
	
	const int N = 6;
	
	//Av = b
	double A[6*6]; // coeficients matrix
	
	int v[6];		// variables' order matrix
	double r[6];	// results' matrix
	
	fill_matricies(A, v, r, N);
	
	solve_matrix(A, v, r, N);

	show_matricies(A, v, r, 1, N);
	
	return 0;
}
