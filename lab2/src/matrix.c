#include "matrix.h"

/*ZADANIE 1*/
void solve_matrix(double A[], int v[], double r[], const int N) {
	
	
	int X, Y;
	for(int k=0; k<N; k++) {
#ifdef DEBUG
		show_matricies(A, v, r, 1, N);
#endif

		/*******************************/
		//PIVOTING
		/*******************************/
		find_pivot(A, k, &X, &Y, N);
		
		// swapping rows
		// must swap rows in resuls matrix too
		if (k != Y) {
			for(int x=k; x<N; x++)			// start from k - no need to swap in columns where value is already 0.0
				swap(&A[x + N*k], &A[x + N*Y]);
			swap(&r[k], &r[Y]);
		}
		
		// swapping columns
		// must swap rows in variables matrix too
		if (k != X) {
			for(int y=0; y<N; y++)		// swap whole column
				swap(&A[k + N*y], &A[X + N*y]);
			swap_int(&v[k], &v[X]);
		}
		// at this point biggest value in matrix A is at A[k][k]
		
		
		
		/*******************************/
		//SCALING
		/*******************************/
		// if this will throw divBy0Fault that will mean that system of equations was indeterminate
		double pivot_scl = 1/A[k + N*k];
		
		// scaling row
		// must scale row in resuls matrix too
		for(int x=k+1; x<N; x++)
			A[x + N*k] *= pivot_scl;
		r[k] *= pivot_scl;
		A[k + N*k] = 1.0;
		// at this point minor's pivot = 1.0 and all other values in row are <= 1.0

		
		/*******************************/
		//ZEROING LOWER TRIANGLE
		/*******************************/
		// zero every cell bellow pivot
		// and for every cell scale it's corresponding row
		for(int y=k+1; y<N; y++) { 		// for remaining rows 
			pivot_scl = A[k + N*y];
			
			for(int x=k+1; x<N; x++)	// scale corresponding row
				A[x + N*y] -= pivot_scl*A[x + N*k];
			r[y] -= pivot_scl*r[k];
			
			A[k + N*y] = 0.0;				// zero it by hand(never trust floating point operations)
		}
		// at this point every value bellow pivot = 0.0
	}
	// at this point matrix A is upper triangular one
	
	
	
	/*******************************/
	//ZEROING UPPER TRIANGLE
	/*******************************/
	// start from bottom-right and go up. Repeat for remaining columns
	// must scale result matrix as well
	for(int x=N-1; x>0; x--)
		for(int y=x-1; y>=0; y--) {
			r[y] -= r[x]*A[x + N*y];
			A[x + N*y] = 0.0;			
		}
	// at this point A is identity matrix, v matrix stores 
	// order of variables with its corresponding values in matrix r
	
	return;
}

// inplace decomposition
void lu_decomp(double A[], const int N) {
	
	for(int k=0; k<N-1; k++) {
		// normalizacja kolumny
		for(int y=k+1; y<N; y++)
			A[k + N*y] /= A[k + N*k];
		
		// modyfikwanie podmacierzy
		for(int y=k+1; y<N; y++)
			for(int x=k+1; x<N; x++)
				A[x + N*y] -= A[x + N*k]*A[k + N*y];
	}
	return;
}




void lu_decomp_split(double A[], double L[], double U[], const int N) {
	// zerowanie macierzy
	for(int y=0; y<N; y++) {
		for(int x=0; x<N; x++)
			L[x + N*y] = U[x + N*y] = 0.0;
		L[y + N*y] = 1.0;
	}
	
	for(int y=0; y<N; y++)
		for(int x=0; x<y; x++)
			L[x + N*y] = A[x + N*y];

	for(int y=0; y<N; y++)
		for(int x=y; x<N-y; x++)
			U[x + N*y] = A[x + N*y];
			
	return;
}




void fill_matricies(double A[], int x[], double y[], const int N) {
	
	srand(time(NULL));
	
	//set coeficients
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
			A[i + N*j] = (double)(rand()%1000)/(rand()%100+1) - (double)(rand()%1000)/(rand()%100+1);
	
	//set initial viaraible's order
	for(int i=0; i<N; i++)
		x[i] = i;
	
	//set results
	for(int i=0; i<N; i++)
		y[i] = (double)(rand()%1000)/(rand()%100+1) - (double)(rand()%1000)/(rand()%100+1);
	
	
	return;
}

void fill_matrix(double A[], const int N) {
	
	srand(time(NULL));
	
	//set coeficients
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
			A[i + N*j] = (double)(rand()%1000)/(rand()%100+1) - (double)(rand()%1000)/(rand()%100+1);
	
	return;
}


void identity_matrix(double A[], const int N) {
	for(int y=0; y<N; y++) {
		for(int x=0; x<N; x++)
			A[x + N*y] = 0;
		A[y + N*y] = 1;
	}
	
	return;
}




void find_pivot(double matrix[], int k, int *X, int *Y, const int N) {
	
	double max = 0.0;
	
	for(int y=k; y<N; y++) 
		for(int x=k; x<N; x++)
			if (max < abs(matrix[x + N*y])) {
				max = matrix[x + N*y];
				*X = x;
				*Y = y;
			}
	return;
}


void show_matricies(double A[], int v[], double r[], int show_A, const int N) {
	
	putchar('\n');
	
	if(show_A)
		for (int y=0; y<N; y++) {
			putchar('|');
			for(int x=0; x<N; x++)
				printf("%4.1f ", A[x + N*y]);
			printf("| \t|%3d| \t|%f|\n", v[y], r[y]);
		}
	else
		for(int y=0; y<N; y++) 
			printf("|%3d| \t|%f|\n", v[y], r[y]);
	
	return;
}


void
show_matrix(double A[], const int N) {
	
	printf("\n{");

	for (int y=0; y<N; y++) {
		printf("\n  {");
		for(int x=0; x<N; x++)
			printf(" % 4.2f, ", A[x + N*y]);
		printf("\b\b },");
	}
	
	printf("\n}\n");
	return;
}








