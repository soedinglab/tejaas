/*
 * =====================================================================================
 *
 *       Filename:  lapacke_dgesvd.c
 *
 *    Description:  Test
 *
 *        Version:  1.0
 *        Created:  12.06.2018 19:17:55
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Saikat Banerjee (banskt), bnrj.saikat@gmail.com
 *   Organization:  Max Planck Institute for Biophysical Chemistry
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <stdio.h>
#include "mkl_lapacke.h"

#define min(a,b) ((a)>(b)?(b):(a))

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );

/* Parameters */
/*
#define M 5000
#define N 338
#define LDA N
#define LDU M
#define LDVT N
*/

/* Main program */
int main() {
        /* Locals */
        printf("Ok, so start \n");
        MKL_INT m, n, lda, ldu, ldvt, info;
        m = 10000;
        n = 338;
        lda = n;
        ldu = m;
        ldvt = n;


        double superb[min(m,n)-1];
        /* Local arrays */
        printf("Dimension: %lld\n", m * n);
        //double* S;
        //double* U;
        //double* VT;
        //double* A;

        printf("Assigning arrays.\n");
        printf ("Size required:  %f Gb\n", (double)((unsigned long) m * m * sizeof(double)) / (1024 * 1024 * 1024));
        //A  = (double *)mkl_malloc( (unsigned long) m * n * sizeof( double ), 64 );
        //U  = (double *)mkl_malloc( (unsigned long) m * m * sizeof( double ), 64 );
        //VT = (double *)mkl_malloc( (unsigned long) n * n * sizeof( double ), 64 );
        //S  = (double *)mkl_malloc( n * sizeof( double ), 64 );
        //double s[n], u[ldu*m], vt[ldvt * n];
        //
        double S[n], A[m * n], VT[n * n];
        /*
        if ( A == NULL || S == NULL || U == NULL || VT == NULL ) {
                printf( "C Error: Can't allocate memory. Aborting... \n");
		mkl_free(A);
		mkl_free(U);
		mkl_free(VT);
		mkl_free(S);
		exit(0);
	}
	*/

        printf("Created arrays.\n");
        //double a[lda * m];
        for ( int i = 0; i < m * n; i++) {
		A[i] = rand();
	}
        /*
         = {
            8.79,  9.93,  9.83, 5.45,  3.16,
            6.11,  6.91,  5.04, -0.27,  7.98,
           -9.15, -7.93,  4.86, 4.85,  3.01,
            9.57,  1.64,  8.83, 0.74,  5.80,
           -3.49,  4.02,  9.80, 10.00,  4.27,
            9.84,  0.15, -8.99, -6.02, -5.31
        }; */
        /* Executable statements */
        //printf( "LAPACKE_dgesvd (row-major, high-level) Example Program Results\n" );
        /* Compute SVD */
        //info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'A', 'A', m, n, a, lda,
        //                s, u, ldu, vt, ldvt, superb );
        /* Check for convergence */
        /*if( info > 0 ) {
                printf( "The algorithm computing SVD failed to converge.\n" );
                exit( 1 );
        }*/
        /* Print singular values */
        //print_matrix( "Singular values", 1, n, s, 1 );
        /* Print left singular vectors */
        //print_matrix( "Left singular vectors (stored columnwise)", m, n, u, ldu );
        /* Print right singular vectors */
        //print_matrix( "Right singular vectors (stored rowwise)", n, n, vt, ldvt );
        /*
        mkl_free(A);
        mkl_free(U);
        mkl_free(VT);
        mkl_free(S);
	*/
        exit( 0 );
} /* End of LAPACKE_dgesvd Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
        MKL_INT i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
                printf( "\n" );
        }
}
