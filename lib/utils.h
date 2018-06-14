/*
 * =====================================================================================
 *
 *       Filename:  utils.h
 *
 *    Description:  Helpers for common operations
 *
 *        Version:  1.0
 *        Created:  12.06.2018 18:59:50
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Saikat Banerjee (banskt), bnrj.saikat@gmail.com
 *   Organization:  Max Planck Institute for Biophysical Chemistry
 *
 * =====================================================================================
 */

#include <stdbool.h>

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  min
 *  Description:  Returns minimum of two integers  
 * =====================================================================================
 */
	int
min ( int i, int j )
{
	if(i < j) return i;
	else return j;
}               /* -----  end of function min  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  transpose
 *  Description:  Transpose of a matrix row-major matrix A (m x n)
 * =====================================================================================
 */
	bool
transpose ( double *A, int m, int n, double *AT )
{
	/*
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			AT[ j * m + i ] = A[ i * n + j ];
		}
	}
	*/
	mkl_domatcopy ('R', 'T', m, n, 1.0, A, n, AT, m);
	return true;
}		/* -----  end of function transpose  ----- */
