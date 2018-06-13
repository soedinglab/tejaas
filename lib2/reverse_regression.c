/*
 * =====================================================================================
 *
 *       Filename:  reverse_regression.c
 *
 *    Description:  Q-score statistics and null model 
 *                  for reverse regression of trans-eQTLs
 *
 *        Version:  1.0
 *        Created:  13.06.2018 11:08:24
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Saikat Banerjee (banskt), bnrj.saikat@gmail.com
 *   Organization:  Max Planck Institute for Biophysical Chemistry
 *
 * =====================================================================================
 */

#include <stdbool.h>                            /* datatype bool */
#include <mkl.h>                                /* mkl_malloc */
#include <math.h>                               /* sqrt */

#include "svd.h"
#include "utils.h"                              /* min, transpose */
#include "chisqr.h"

#define MAF_NULL 1
#define PERM_NULL 2

void genotype_variance ( double* GT, int nsnp, int nsample, double* SX2, int null );
bool getSmod( double* S, double* Smod, double* sx2, double* sb2, int nsnp, int nS );
bool getW(double* U, double* Smod, double* W, int nsample, int nS, int ioff);
double vecT_smat_vec ( int n, double* v, double* A, double* D );
double permuted_null ( int N, double* X, double* W, double Q, double* muQ, double* sigQ, int nsample );

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  qscore
 *  Description:  Interfaced with python
 * =====================================================================================
 */
	bool
qscore ( double* GT, double* GX, double* SB2, int ngene, int nsnp, int nsample, int null, double* MAF, double* Q, double* P, double* MUQ, double* SIGQ )
{
	bool success;
	int info;
	int nS;                                     /* number of components of SVD */
	double dQ;

	success = false;
	nS = min(ngene, nsample);

	double *GXT = (double *) mkl_malloc( (unsigned long) ngene   * nsample * sizeof( double ), 64 );
	if (GXT == NULL) {success = false; goto cleanup_GXT;}

	double *S   = (double *) mkl_malloc(                           nS      * sizeof( double ), 64 );
	if (S == NULL) {success = false; goto cleanup_S;}

	double *U   = (double *) mkl_malloc( (unsigned long) nsample * nS      * sizeof( double ), 64 );
	if (U == NULL) {success = false; goto cleanup_U;}

	double *SM  = (double *) mkl_malloc( (unsigned long) nsnp    * nS      * sizeof( double ), 64 );
	if (SM == NULL) {success = false; goto cleanup_SM;}

	double *W   = (double *) mkl_malloc( (unsigned long) nS      * nS      * sizeof( double ), 64 );
	if (W == NULL) {success = false; goto cleanup_W;}

	double *X   = (double *) mkl_malloc(                           nsample * sizeof( double ), 64 );
	if (X == NULL) {success = false; goto cleanup_X;}

	double *SX2 = (double *) mkl_malloc(                           nsnp    * sizeof( double ), 64 );
	if (SX2 == NULL) {success = false; goto cleanup_SX2;}

	double *D1  = (double *) mkl_malloc(                           nS      * sizeof( double ), 64 );
	if (D1 == NULL) {success = false; goto cleanup_D1;}

	

	success = transpose(GX, ngene, nsample, GXT);
	if ( success == false ) goto cleanup;

	success = dsvd (GXT, nsample, ngene, S, U);
	if ( success == false ) goto cleanup;

	genotype_variance(GT, nsnp, nsample, SX2, null);
	printf("%g\n", SX2[2]);

	success = getSmod (S, SM, SX2, SB2, nsnp, nS);
	if ( success == false ) goto cleanup;

	if (null == MAF_NULL) {                            /* fixed W for all SNPs with MAF null, because sigmax2 = 1 */
		success = getW (U, SM, W, nsample, nS, 0);
		if ( success == false ) goto cleanup;
	}

	for ( int i = 0; i < nsnp; i++ ) {

		if (null == PERM_NULL) {                        /* gets W for every SNP if using permutation null, because sigmax2 is not fixed */
			success = getW (U, SM, W, nsample, nS, i*nS);
			if ( success == false ) goto cleanup;
		}

		for ( int j = 0; j < nsample; j++) {
			X[j] = GT[ i*nsample + j ];
		}
		Q[i] = vecT_smat_vec ( nS, X, W, D1 );

		if ( null == PERM_NULL ) {
			P[i] = permuted_null ( nS, X, W, Q[i], &MUQ[i], &SIGQ[i], nsample );
		}
		else if ( null == MAF_NULL ) {
			//success = maf_null ();
			P[i] = 1;
		}

	}

cleanup:
cleanup_D1:
	mkl_free(D1);
cleanup_SX2:
	mkl_free(SX2);
cleanup_X:
	mkl_free(X);
cleanup_W:
	mkl_free(W);
cleanup_SM:
	mkl_free(SM);
cleanup_U:
	mkl_free(U);
cleanup_S:
	mkl_free(S);
cleanup_GXT:
	mkl_free(GXT);

	return success;

}		/* -----  end of function qscore  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  second_moment
 * =====================================================================================
 */
	double
second_moment ( double* A, int n )
{
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += A[i] * A[i];
	}
	return sum / n;
}		/* -----  end of function second_moment  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fourth_moment
 * =====================================================================================
 */
    double
fourth_moment ( double* A, int n )
{
    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum += A[i] * A[i] * A[i] * A[i];
    }
    return sum / n;
}       /* -----  end of function fourth_moment  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  permuted_null
 *  Description:  Calculate p-values for permuted null model
 * =====================================================================================
 */
	double
permuted_null ( int N, double* X, double* W, double Q, double* muQ, double* sigQ, int nsample )
{
	bool success;
    double pval, sig2;

	double q11 = 0, q2 = 0, q31 = 0, q4 = 0, q22 = 0, q211 = 0, v31 = 0, v22 = 0, v211 = 0, v1111 = 0;

	double *q31_left  = (double *) mkl_malloc( N * sizeof( double ), 64 );
	if (q31_left == NULL)  {success = false; goto cleanup_q31_left;}

	double *q31_right = (double *) mkl_malloc( N * sizeof( double ), 64 );
	if (q31_right == NULL) {success = false; goto cleanup_q31_right;}

    double gtmu2  = second_moment(X, nsample);
    double gtmu4  = fourth_moment(X, nsample);


	for( int i=0; i < N; i++ ){
		q2 += W[ i*N + i ];
		q4 += W[ i*N + i ] * W[ i*N + i ];

		double q211_inner = 0;
		for( int j=0; j < N; j++ ){
			q11 += W[ i*N + j ];
			q22 += W[ i*N + j ] * W[ i*N + j ];
			q211_inner = W[ i*N + j ];
		}

		q31_right[i] = q211_inner;
		q31_left [i] = W[ i*N + i ];
		q211 += q211_inner * q211_inner;
	}

	q31 = cblas_ddot (N, q31_left, 1, q31_right, 1);

	v31     = - gtmu4 / (N - 1);
	v22     = (N * gtmu2 * gtmu2 / (N - 1)) + v31;
	v211    = - (v31 + v22) / (N - 2);
	v1111   = - 3 * v211 / (N - 3);

	*muQ   = (gtmu2 / (N - 1)) * (N * q2 - q11);
	sig2 = v1111 * (q11*q11 - 2*q2*q11 - 4*q211 + 8*q31 + 2*q22 + q2*q2 - 6*q4)
				+ 2 * v211 * (q2 * q11 + 2*q211 - 6*q31 - 2*q22 - q2*q2 + 6*q4)
				+ v22 * (q2*q2 + 2*q22 - 3*q4)
				+ 4 * v31 *(q31 - q4)
				+ gtmu4 * q4;

	sig2  =  sig2 - (*muQ) * (*muQ);
	*sigQ = sqrt(sig2);
	double mscale = (sig2 / *muQ ) / 2.0;
	double df = *muQ / mscale;
	double Qscaled = Q / mscale;
	pval = chisqr((int)df, (double)Qscaled);

cleanup:
cleanup_q31_right:
	mkl_free(q31_right);
cleanup_q31_left:
	mkl_free(q31_left);

	return pval;
}		/* -----  end of function permuted_null  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getSmod
 *  Description:  Return Smod from S and sigma_beta2 and sigma_x2
 * =====================================================================================
 */
	bool
getSmod ( double* S, double* Smod, double* sx2, double* sb2, int nsnp, int nS )
{
	bool success;
	double S2i;
	for ( int i = 0; i < nS; i++ ) {
		S2i = S[i] * S[i];
		for ( int j = 0; j < nsnp; j++ ) {
			Smod[ j * nS + i ] = S2i / (S2i + (sx2[j] / sb2[j]));
		}
	}

	return true;
}		/* -----  end of function getSmod  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getW
 *  Description:  Return the W
 * =====================================================================================
 */
	bool
getW ( double* U, double* Smod, double* W, int nsample, int nS, int ioff )
{
	for ( int i = 0; i < nsample; i++ ) {
		for ( int j = 0; j < nS; j++ ) {
			W[ i*nS + j ] = 0;
			for ( int k = 0; k < nS; k++ ) {
				W[ i*nS + j ] += U[ i*nS + k ] * U[ j*nS + k ] * Smod[ ioff + k ];
			}
		}
	}
	return true;
}		/* -----  end of function getW  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  genotype_variance
 *  Description:  
 * =====================================================================================
 */
	void
genotype_variance ( double* GT, int nsnp, int nsample, double* SX2, int null )
{
	double mean, var;
	for ( int i = 0; i < nsnp; i++ ) {
		var = 1;
		if (null == PERM_NULL) {
			mean = 0;
			for ( int j = 0; j < nsample; j++ ) {
				mean += GT[ i*nsample + j ];
			}
			mean = mean / nsample;

			var = 0;
			for ( int j = 0; j < nsample; j++ ) {
				var += (GT[ i*nsample + j] - mean) * (GT[ i*nsample + j ] - mean);
			}
			var = var / nsample;
		}
		SX2[i] = var;
	}
	return;
}		/* -----  end of function genotype_variance  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  vecT_smat_vec
 *  Description:  calculate v'Av where v is a vector of length n, 
 *                and A is a n-by-n symmetric matrix
 *                D is a dummy vector of length n
 * =====================================================================================
 */
    double
vecT_smat_vec ( int n, double* v, double* A, double* D )
{
    double one = 1.0;
    double zero = 0.0;
    double res;

    cblas_dsymv(CblasRowMajor, CblasLower, n, one, A, n, v, 1, zero, D, 1);
    res = cblas_ddot(n, v, 1, D, 1);

    return res;
}		/* -----  end of function vecT_smat_vec  ----- */
