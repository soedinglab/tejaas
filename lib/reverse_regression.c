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
 *         			Raktim Mitra (timkartar), timkartar7879@gmail.com
 *   Organization:  Max Planck Institute for Biophysical Chemistry
 *
 * =====================================================================================
 */

#include <stdbool.h>                            /* datatype bool */
#include <math.h>                               /* sqrt, fabs */
#include <time.h>

#include "svd.h"
#include "utils.h"                              /* min, transpose */

#ifdef MKL_ILP64
#include <mkl.h>
void my_cdfnorm( double* X, double* P) {
	vdCdfNorm( 1, X, P );
}
#else
#include <cblas.h>
#include <cblas_f77.h>
#include "dcdflib/src/dcdflib.c"
#include "dcdflib/src/ipmpar.c"
void my_cdfnorm( double* X, double* P) {
	int which[1] = {1}; // iwhich = 1 : Calculate P and Q from X,MEAN and SD
	double Q[1] = {0.0};
	double MEAN[1] = {0.0};
	double SD[1] = {1.0};
	int status[1] = {0};
	double bound[1];
    cdfnor(which, P, Q, X, MEAN, SD, status, bound);
}
#endif

#define NO_NULL 0
#define MAF_NULL 1
#define PERM_NULL 2

bool genotype_variance ( double* GT, int nsnp, int nsample, double* SX2, int null );
bool getSmod( double* S, double* Smod, double* sx2, double* sb2, int nsnp, int nS );
bool getW(double* U, double* Smod, double* W, int nsample, int nS, int ioff, double sx2);
double vecT_smat_vec ( int n, double* v, double* A, double* D );
double permuted_null ( int N, double* X, double* W, double Q, double* muQ, double* sigQ, int nsample );
void getWnullmaf ( double* W, double* SM, double* muQmaf, double* sig2Qmaf, double* sumW2nn, int N, int ioff );
//double qnull_maf ( double Q, double mu, double sigma );
double gt4maf ( double f );
double cdf_norm ( double x, double mu, double sigma );


bool A_vecV(double* A, double* v, double* B, int ngene, int nsample, int ioff);
bool matmulAB( double* A, double* B, double* C, int ngene, int nS, int nsample) ;
bool diagA_B ( double* A, double* B, double* C, int nS, int nsample, int ioff);
bool getLinv_ST ( double* S, double* Linv_ST, double* sx2, double* sb2, int nsnp, int nS );



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
	double muQmaf, sig2Qmaf, sumW2nn;

	success = false;
	nS = min(ngene, nsample);

	double *GXT = (double *) calloc( (unsigned long) ngene   * nsample * sizeof( double ), 64 );
	if (GXT == NULL) {success = false; goto cleanup_GXT;}

	double *S   = (double *) calloc(                           nS      * sizeof( double ), 64 );
	if (S == NULL) {success = false; goto cleanup_S;}

	double *U   = (double *) calloc( (unsigned long) nsample * nS      * sizeof( double ), 64 );
	if (U == NULL) {success = false; goto cleanup_U;}

	double *VT  = (double *) calloc( (unsigned long )  nS * ngene      * sizeof( double ), 64 );
	if (VT == NULL) {success = false; goto cleanup_VT;}

	double *SM  = (double *) calloc( (unsigned long)    nsnp * nS      * sizeof( double ), 64 );
	if (SM == NULL) {success = false; goto cleanup_SM;}

	double *W   = (double *) calloc( (unsigned long)      nS * nS      * sizeof( double ), 64 );
	if (W == NULL) {success = false; goto cleanup_W;}

	double *X   = (double *) calloc(                           nsample * sizeof( double ), 64 );
	if (X == NULL) {success = false; goto cleanup_X;}

	double *SX2 = (double *) calloc(                           nsnp    * sizeof( double ), 64 );
	if (SX2 == NULL) {success = false; goto cleanup_SX2;}

	double *D1  = (double *) calloc(                           nS      * sizeof( double ), 64 );
	if (D1 == NULL) {success = false; goto cleanup_D1;}

	

	success = transpose(GX, ngene, nsample, GXT);
	if ( success == false ) goto cleanup;

	clock_t start, end;
	double cpu_time_used;
	start = clock();
	success = dsvd (GXT, nsample, ngene, S, U, VT);
	if ( success == false ) goto cleanup;

	success = genotype_variance(GT, nsnp, nsample, SX2, null);

	success = getSmod (S, SM, SX2, SB2, nsnp, nS);
	if ( success == false ) goto cleanup;

	if (null == MAF_NULL) {                            /* fixed W for all SNPs with MAF null, because sigmax2 = 1 */
		success = getW (U, SM, W, nsample, nS, 0, 1);
		if ( success == false ) goto cleanup;
		getWnullmaf ( W, SM, &muQmaf, &sig2Qmaf, &sumW2nn, nS, 0);
	}

	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	// printf("SVD calculation took %f seconds \n", cpu_time_used);

	for ( int i = 0; i < nsnp; i++ ) {

		if (null == NO_NULL ) {
			// duplicated from PERM_NULL, what should we do by default?
			success = getW (U, SM, W, nsample, nS, i*nS, SX2[i]);
			if ( success == false ) goto cleanup;
			MUQ[i] = -1.0;
			SIGQ[i] = -1.0;
			P[i] = 1.0;
		}

		if (null == PERM_NULL) {                        /* gets W for every SNP if using permutation null, because sigmax2 is not fixed */
			// Here probably add some check when G < N
		    // [X] and also multiply by 1/SX2[i] because it is not fixed here
		    // DONE!
			success = getW (U, SM, W, nsample, nS, i*nS, SX2[i]);
			if ( success == false ) goto cleanup;
		}

		for ( int j = 0; j < nsample; j++) {
			X[j] = GT[ i*nsample + j ];
		}

		Q[i] = vecT_smat_vec ( nS, X, W, D1 );

		if ( null == PERM_NULL ) {
			P[i] = permuted_null ( nS, X, W, Q[i], &MUQ[i], &SIGQ[i], nsample );
			// printf("Pval: %f \n", P[i]);
		}
		else if ( null == MAF_NULL ) {
			MUQ[i] = muQmaf;
			SIGQ[i] = sqrt( sig2Qmaf + (gt4maf(MAF[i]) - 3) * sumW2nn );
			P[i] = cdf_norm (Q[i], MUQ[i], SIGQ[i] ); //qnull_maf(Q[i], MUQ[i], SIGQ[i]);
		}



	}

	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	// printf("RR calculation took %f seconds \n", cpu_time_used);

cleanup:
cleanup_D1:
	//printf("cleaning D1\n");
	free(D1);
cleanup_SX2:
	//printf("cleaning SX2\n");
	free(SX2);
cleanup_X:
	//printf("cleaning X\n");
	free(X);
cleanup_W:
	//printf("cleaning W\n");
	free(W);
cleanup_SM:
	//printf("cleaning SM\n");
	free(SM);
cleanup_VT:
	//printf("cleaning VT\n");
	free(VT);
cleanup_U:
	//printf("cleaning U\n");
	free(U);
cleanup_S:
	//printf("cleaning S\n");
	free(S);
cleanup_GXT:
	//printf("cleaning GXT\n");
	free(GXT);

	return success;

}		/* -----  end of function qscore  ----- */


bool betas(double* GT, double* GX, double* SB2, int ngene, int nsnp, int nsample, double* B)
{
	bool success;
	int info;
	int nS;                                     /* number of components of SVD */

	success = false;
	nS = min(ngene, nsample);

	double *GXT = (double *) calloc( (unsigned long)  ngene * nsample  * sizeof( double ), 64 );
	if (GXT == NULL) {success = false; goto cleanup_GXT;}

	double *preB = (double *) calloc( (unsigned long)  ngene * nsample  * sizeof( double ), 64 );
	if (preB == NULL) {success = false; goto cleanup_preB;}

	double *Bi = (double *) calloc( (unsigned long)           ngene    * sizeof( double ), 64 );
	if (Bi == NULL) {success = false; goto cleanup_Bi;}

	double *S   = (double *) calloc(                           nS      * sizeof( double ), 64 );
	if (S == NULL) {success = false; goto cleanup_S;}

	double *U   = (double *) calloc( (unsigned long)   nsample * nS    * sizeof( double ), 64 );
	if (U == NULL) {success = false; goto cleanup_U;}

	double *UT   = (double *) calloc( (unsigned long)    nS * nsample  * sizeof( double ), 64 );
	if (UT == NULL) {success = false; goto cleanup_UT;}

	double *VT  = (double *) calloc( (unsigned long )   nS * ngene     * sizeof( double ), 64 );
	if (VT == NULL) {success = false; goto cleanup_VT;}

	double *V  = (double *) calloc( (unsigned long )    ngene * nS     * sizeof( double ), 64 );
	if (V == NULL) {success = false; goto cleanup_V;}

	// contains nsnps vector of size nS (diagonal elements of the S matrix for each snp)
	double *Linv_ST  = (double *) calloc( (unsigned long)  nS * nsnp      * sizeof( double ), 64 );
	if (Linv_ST == NULL) {success = false; goto cleanup_Linv_ST;}

	// this could be extended to contain the matrix for all snps? maybe to much memory (multiply by *nsnp)
	double *Linv_STUT  = (double *) calloc( (unsigned long) nS * nsample  * sizeof( double ), 64 );
	if (Linv_STUT == NULL) {success = false; goto cleanup_Linv_STUT;}

	double *X   = (double *) calloc(                           nsample * sizeof( double ), 64 );
	if (X == NULL) {success = false; goto cleanup_X;}

	double *SX2 = (double *) calloc(                           nsnp    * sizeof( double ), 64 );
	if (SX2 == NULL) {success = false; goto cleanup_SX2;}

	success = transpose(GX, ngene, nsample, GXT);
	if ( success == false ) goto cleanup;

	clock_t start, end;
	double cpu_time_used;
	start = clock();
	success = dsvd (GXT, nsample, ngene, S, U, VT);
	if ( success == false ) goto cleanup;

	success = transpose(U, nsample, nS, UT);
	if ( success == false ) goto cleanup;

	success = transpose(VT, nS, ngene, V);
	if ( success == false ) goto cleanup;

	success = genotype_variance(GT, nsnp, nsample, SX2, PERM_NULL);
	getLinv_ST (S, Linv_ST, SX2, SB2, nsnp, nS);

	for ( int i = 0; i < nsnp; i++ ) {
		for ( int j = 0; j < nsample; j++) {
			X[j] = GT[ i*nsample + j ];
		}

		// Do Linv_ST * UT
		diagA_B ( Linv_ST, UT, Linv_STUT, nS, nsample, i*nS);

		// Now do V * Linv_STUT
		matmulAB( V, Linv_STUT, preB, ngene, nS, nsample);

		// Finally, multipliy preB * X[i]
		// preB is ngene x nsample
		A_vecV(preB, X, B, ngene, nsample, i*ngene);
	}

	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	// printf("Betas calculation took %f seconds \n", cpu_time_used);
	success = true;


cleanup:
cleanup_SX2:
	free(SX2);
cleanup_X:
	free(X);
cleanup_Linv_STUT:
	free(Linv_STUT);
cleanup_Linv_ST:
	free(Linv_ST);
cleanup_V:
	free(V);
cleanup_VT:
	free(VT);
cleanup_UT:
	free(UT);
cleanup_U:
	free(U);
cleanup_S:
	free(S);
cleanup_Bi:
	free(Bi);
cleanup_preB:
	free(preB);
cleanup_GXT:
	free(GXT);

	return success;

}

// A is ngene x nsample
// v is nsample x 1
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  A_vecV
 *  Description:  Multiplies matrix A of size m x n
 *                with vector V of length n x 1
 *                and outputs vector B of size m x 1
 * =====================================================================================
 */
bool A_vecV(double* A, double* v, double* B, int m, int n, int ioff)
{
	bool success;
	double alpha = 1.0;
	double beta  = 0.0;
/*	int m = ngene;
	int n = nsample;*/
	double * y;
	y = (double *) calloc( (unsigned long) (m - 1) * sizeof( double ), 64 );
	if (y == NULL) {success = false; goto cleanup_y;}

	for (int i=0; i < m; i++) {
		y[i] = 0.0;
	}

	cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, alpha, A, n, v, 1, beta, y, 1 );

	for (int i=0; i<m; i++)
	{
		B[ioff + i] = y[i];
		//printf("B[%d]=%f - ", ioff +i, B[ioff+i]);
	}
	success = true;

	cleanup_y:
		free(y);
	return success;
}


// V is          ngene x nS
// Linv_STUT is  nS x nsample

// VT is k x G
// V  is G x k
// LinvSTUT is k x N 
// bool matmulAB( double* V, double* Linv_STUT, double* preB, int ngene, int nS, int nsample) 
// returns G x N matrix

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  matmulAB
 *  Description:  Multiplies matrix A of size m x k
 *                with matrix B of size k x n and returns 
 *                matrix C of size m x n
 * =====================================================================================
 */

bool matmulAB( double* A, double* B, double* C, int m, int k, int n) 
{
	double alpha = 1.0;
	double beta  = 0.0;
/*	int m = ngene;
	int n = nsample;
	int k = nS;*/

	for (int i = 0; i < (m*n); i++) {
        C[i] = 0.0;
    }

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                m, n, k, alpha, A, k, B, n, beta, C, n);

    return true;

}



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  diagA_B
 *  Description:  Multiplies a diagonal matrix A of size m x m (nS x nS) which is actually a vector of length m
 				  with dense matrix B of size m x n ( nS x nsample )
 				  and outputs C of size (m x n) (nS x nsample)
 * =====================================================================================
 */

bool diagA_B ( double* A, double* B, double* C, int nS, int nsample, int ioff)
{
	for (int i = 0; i < nS; i++) {
			for (int j = 0; j < nsample; j++) {
				C[ i*nsample + j ] = A[ioff + i] * B[i*nsample + j];
			}
		}
	return true;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getLinv_ST
 *  Description:  Return L_inv * ST from S and sigma_beta2 and sigma_x2
 *                (STS + sigma_x2/sigma_b2)^-1 * ST
 * =====================================================================================
 */
	bool
getLinv_ST ( double* S, double* Linv_ST, double* sx2, double* sb2, int nsnp, int nS )
{
	double S2i;
	double Si;
	for ( int i = 0; i < nS; i++ ) {
		Si  = S[i];
		S2i = Si * Si;
		for ( int j = 0; j < nsnp; j++ ) {
			Linv_ST[ j * nS + i ] =  Si / (S2i + (sx2[j] / sb2[j]));
		}
	}

	return true;
}		/* -----  end of function getLinv_ST  ----- */



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
 *         Name:  gt4maf
 *  Description:  
 * =====================================================================================
 */
	double
gt4maf ( double f )
{
	double f0 = (1 - f) * (1 - f);
	double f1 = 2.0 * f * (1 - f);
	double f2 = f * f;
	double mu = 2 * f;
	double sig = sqrt(f1);
	double x4 = f0 *  pow((-mu/sig), 4) + f1 * pow(((1 - mu)/sig), 4) + f2 * pow(((2 - mu)/sig), 4);

	return x4;
}		/* -----  end of function gt4maf  ----- */


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

	double *q31_left  = (double *) calloc( N * sizeof( double ), 64 );
	if (q31_left == NULL)  {success = false; goto cleanup_q31_left;}

	double *q31_right = (double *) calloc( N * sizeof( double ), 64 );
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
			q211_inner += W[ i*N + j ];
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

	// Due to some numerical error with large sig_b, the variance can be negative, small fix
	if (sig2 < 0) {
		printf("Warning: negative sigQ: %f", sig2);
		sig2 = -sig2;
	}

	*sigQ = sqrt(sig2);  // sometimes this turns out -nan
	pval = cdf_norm (Q, *muQ, *sigQ);


cleanup:
cleanup_q31_right:
	free(q31_right);
cleanup_q31_left:
	free(q31_left);

	return pval;
}		/* -----  end of function permuted_null  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  qnull_maf
 *  Description:  
 * =====================================================================================
	double
qnull_maf ( double Q, double mu, double sigma )
{
	double sig2 = sigma * sigma;
	double mscale = (sig2 / mu) / 2.0;
	double df = mu / mscale;
	double Qscale = Q / mscale;
	double pval = chisqr( (int)df, (double)Qscale );
	printf ("M: %g, S: %g, Q: %g, Qs: %g, df: %g, p: %g\n", mu, sigma, Q, Qscale, df, pval);
	return pval;
}
 */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  cdf_norm
 *  Description:  
 * =====================================================================================
 */
	double
cdf_norm ( double x, double mu, double sigma )
{

	double pval;

	double *QS = (double *) calloc( 1 * sizeof( double ), 64);
	if (QS == NULL) {goto cleanup_QS;}

	double *P = (double *) calloc( 1 * sizeof( double ), 64);
	if (P == NULL) {goto cleanup_P;}

	/** QS[0] = fabs((x - mu) / sigma); 
	P[0] = 0;
	my_cdfnorm( QS, P );
	pval = 2.0 * (1 - P[0]);
        **/ // we don't need 2-sided p-value
        QS[0] = (x - mu) / sigma; 
        P[0] = 0;
        my_cdfnorm( QS, P );
        pval = (1 - P[0]);
	
	//printf ("M: %g, S: %g, Q: %g, p: %g\n", mu, sigma, x, pval);

cleanup:
cleanup_QS:
	free(QS);
cleanup_P:
	free(P);

	return pval;
}		/* -----  end of function qnull_maf  ----- */

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
 *                success = getW (U, SM, W, nsample, nS, i*nS);
 * =====================================================================================
 */
	bool
getW ( double* U, double* Smod, double* W, int nsample, int nS, int ioff, double sx2 )
{
	for ( int i = 0; i < nsample; i++ ) {
		for ( int j = 0; j < nS; j++ ) {
			W[ i*nS + j ] = 0;
			for ( int k = 0; k < nS; k++ ) {
				W[ i*nS + j ] += U[ i*nS + k ] * U[ j*nS + k ] * Smod[ ioff + k ];
			}
			W[ i*nS + j ] /= sx2;
		}
	}

	return true;
}		/* -----  end of function getW  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getWnullmaf
 *  Description:  
 * =====================================================================================
 */
	void
getWnullmaf ( double* W, double* Smod, double* muQmaf, double* sig2Qmaf, double* sumW2nn, int N, int ioff )
{
	double mu = 0;
	double sig2 = 0;
	double sum = 0;
	
	for ( int i = 0; i < N; i++ ) {
		mu += Smod [ioff + i];
		sig2 += 2.0 * Smod[ioff + i] * Smod[ioff + i];
		sum += W[ i*N + i ] * W[ i*N + i ];
	}

	*muQmaf = mu;	
	*sig2Qmaf = sig2;
	*sumW2nn = sum;

	return;
}		/* -----  end of function getWnullmaf  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  genotype_variance
 *  Description:  
 * =====================================================================================
 */
	bool
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
	return true;
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
