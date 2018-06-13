/*
 * =====================================================================================
 *
 *       Filename:  trial.c
 *
 *    Description:  Trial
 *
 *        Version:  1.0
 *        Created:  12.06.2018 19:57:52
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Saikat Banerjee (banskt), bnrj.saikat@gmail.com
 *   Organization:  Max Planck Institute for Biophysical Chemistry
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include <mkl_lapacke.h>

void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );

int main() {

	printf("Ok, so start \n");
	int m, n, lda, ldu, ldvt, info;
    m = 30000;
    n = 338;
    lda = n;
    ldu = m;
    ldvt = n;

    double * S;
    double * A;
    double * U;
    double * VT;
    S  = (double *) mkl_malloc( (unsigned long) n * sizeof( double ), 64 );
    A  = (double *) mkl_malloc( (unsigned long) m * n * sizeof( double ), 64 );
    U  = (double *) mkl_malloc( (unsigned long) m * m * sizeof( double ), 64 );
    VT = (double *) mkl_malloc( (unsigned long) n * n * sizeof( double ), 64 );


    double * superb;
    superb = (double *) mkl_malloc( (unsigned long) (n - 1) * sizeof( double ), 64 );

	for ( int i = 0; i < m * n; i++) {
		A[i] = rand();
	}


    printf("Calculating SVD.\n");
    info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'A', 'A', m, n, A, lda, S, U, ldu, VT, ldvt, superb );
    printf( "LAPACKE_dgesvd (row-major, high-level) Example Program Results\n" );
    print_matrix( "Singular values", 1, n, S, 1 );

    //double S[n];//, A[m * n], VT[n * n];

    printf("Done.\n");

}

void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
    int i, j;
    printf( "\n %s\n", desc );
    for( i = 0; i < m; i++ ) {
        for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
        printf( "\n" );
    }
}
