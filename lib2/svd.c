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
#include <mkl.h>
#include <mkl_lapacke.h>
#include "utils.h"

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  dsvd
 *  Description:  SVD of rectangular double precision matrix
 *                A is matrix of size m x n
 *                A = USV' [U size m x k, S size k x k, V size n x k] ... k = min(m, n)
 *                Returns U and S
 * =====================================================================================
 */
	bool
dsvd ( double *A, int m, int n, double *S, double *U )
{
	bool success = false;
	int info;
	int k = min(m, n);
	
	double * superb;
	double * VT;
	superb = (double *) mkl_malloc( (unsigned long) (k - 1) * sizeof( double ), 64 );
	if (superb == NULL) {success = false; goto cleanup_gxT;}

	VT     = (double *) mkl_malloc( (unsigned long) k * n   * sizeof( double ), 64 );
    if (VT == NULL) {success = false; goto cleanup_gxT;}
	
	info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'S', 'S', m, n, A, n, S, U, k, VT, n, superb);
	
	if ( info != 0 ) {
		success = false;
		printf ("Error in calculating SVD with info = %g\n", info);
		goto cleanup;
	} else {
		success = true;
	}

cleanup:
cleanup_VT:
	mkl_free(VT);
cleanup_superb:
	mkl_free(superb);

	return success;
}		/* -----  end of function dsvd  ----- */
