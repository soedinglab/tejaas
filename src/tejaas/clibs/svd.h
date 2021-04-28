/*
 * =====================================================================================
 *
 *       Filename:  svd.h
 *
 *    Description:  Calculation of SVD using Intel MKL 
 *
 *        Version:  1.0
 *        Created:  12.06.2018 14:41:44
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Saikat Banerjee (banskt), bnrj.saikat@gmail.com
 *   Organization:  Max Planck Institute for Biophysical Chemistry
 *
 * =====================================================================================
 */

#include <stdbool.h>
#include <stdio.h>

#ifdef MKL_ILP64
#include <mkl.h>
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
/*
	static int
LAPACKE_dgesvd ( int matrix_layout, char jobu, char jobvt,
			int m, int n, double* a,
			int lda, double* s, double* u, int ldu,
			double* vt, int ldvt, double* superb ) {
	int info = 0;
	int lwork = -1;
	double* work = NULL;
	double work_query;

	// Query optimal work array
	dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
	lwork = (int) work_query;

	
	return info;
}
*/
#endif

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  dsvd
 *  Description:  SVD of rectangular double precision matrix
 *                A is matrix of size m x n
 *                A = USV' [U size m x k, S size k x k, V size n x k] ... k = min(m, n)
 *                Returns U, S and VT
 * =====================================================================================
 */
	bool
dsvd ( double *A, int m, int n, double *S, double *U, double *VT )
{
	bool success = false;
	int info;
	int k;
	if (m < n) k = m;
	else k = n;
	
	double * superb;
	superb = (double *) calloc( (unsigned long) (k - 1) * sizeof( double ), 64 );
	if (superb == NULL) {success = false; goto cleanup_superb;}

	info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'S', 'S', m, n, A, n, S, U, k, VT, n, superb);
	
	if ( info != 0 ) {
		success = false;
		printf ("Error in calculating SVD with info = %d\n", info);
		goto cleanup;
	} else {
		success = true;
	}

cleanup:
cleanup_superb:
	free(superb);

	return success;
}		/* -----  end of function dsvd  ----- */
