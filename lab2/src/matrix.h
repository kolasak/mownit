#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// main functions
void solve_matrix(double A[], int v[], double r[], const int N);
void lu_decomp(double A[], const int N);
void find_pivot(double matrix[], int k, int *X, int *Y, const int N);

// matrix initializatoin functions
void fill_matricies(double A[], int v[], double r[], const int N);
void fill_matrix(double A[], const int N);
void identity_matrix(double A[], const int N);

// visualization fucntions
void show_matricies(double A[], int v[], double r[], int show_A, const int N);
void show_matrix(double A[], const int N);
void lu_decomp_split(double A[], double L[], double U[], const int N);


static inline void swap(double *a, double *b) {
	double tmp = *a;
	*a = *b;
	*b = tmp;
	
	return;
}

static inline void swap_int(int *a, int *b) {
	int tmp = *a;
	*a = *b;
	*b = tmp;
	
	return;
}
