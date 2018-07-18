#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#ifdef MKL_ILP64
#include <mkl.h>
#else
#include <cblas.h>
#include <cblas_f77.h>
#endif

#define DEBUG false

int add( double* arr, double* res, int n) {
	int i;
	res[0] = cblas_dasum ( n, arr, 1);
/*
	res[0] = 0;
	for ( i = 0; i < n; i++ ) {
		res[0] += arr[i];
	}
*/
	return 0;
}

void printvec(double* arr, int n) {
	int i;
	printf("Vector:");
	for ( i = 0; i < n; i++ ) {
		printf (" %g", arr[i]);
	}
	printf("\n");
}

int main() {

	int i, isucc;

	int n;
	n = 5000000;

	clock_t begin = clock();

	double *A = (double *) calloc( n * sizeof( double ), 64 );
	if ( A == NULL ) {goto cleanup_A;}

	double *S = (double *) calloc( 1 * sizeof( double ), 64 );
	if ( S == NULL ) {goto cleanup_S;}

	for ( i = 0; i < n; i++ ) {
		//A[i] = (double) i * 2.0;
		A[i] = (double) pow(-1, i);
	}

	if (DEBUG) {printvec(A, n);}

	isucc = add(A, S, n);
	printf ("Sum: %g\n", S[0]);

cleanup:
cleanup_S:
	free(S);
cleanup_A:
	free(A);

	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf ("Time taken: %f sec.\n", time_spent);

	return 0;
}
