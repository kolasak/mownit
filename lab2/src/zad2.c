#include "matrix.h"


int main(void) {
	
	const int N = 5;
	
	double A[] = {  2.3,  19.0,  11.0,  -1.7,   5.4,
				   11.0,  12.5,  18.0,   6.0,  -5.0,
				    8.1,   8.1,   5.8,  12.7,  14.1,
				   -9.0,  -9.2,   5.0,   8.0,  11.6,
				    4.0,   9.2,  -4.0,  14.0,  19.0 };
						  
	printf("\nA");
	show_matrix(A, N);
	
	lu_decomp(A, N);
	double L[N*N];
	double U[N*N];
	lu_decomp_split(A, L, U, N);
	
	
	printf("\nL");
	show_matrix(L, N);
	printf("\nU");
	show_matrix(U, N);
	
	return 0;
}
