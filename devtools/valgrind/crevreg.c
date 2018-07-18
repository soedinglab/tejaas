#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include "reverse_regression.h"

#define DEBUG true

double randn (double mu, double sigma);
double randb (double f);
double randu (double a, double b);
void   printGX (double *GX, int ngene, int nsample, int gprint, int sprint);
void   printGT (double *GT, double* MAF, int nsnp,  int nsample, int gprint, int sprint);
void   printvec(char* name, double* arr, int n);

bool main() {
	bool success;
	int ngene;
	int nsample;
	int nsnp;
	int i, j;
	int isucc;

	clock_t begin = clock();
	char mstring[30];
	
	success = false;
	ngene = 23973;
	nsample = 338;
	nsnp = 10;
	
	double *GX    = (double *) calloc( (unsigned long) ngene   * nsample * sizeof( double ), 64 );
	if (GX == NULL) {success = false; goto cleanup_GX;}
	
	double *GT    = (double *) calloc( (unsigned long) nsnp    * nsample * sizeof( double ), 64 );
	if (GT == NULL) {success = false; goto cleanup_GT;}
	
	double *MAF    = (double *) calloc( (unsigned long) nsnp             * sizeof( double ), 64 );
	if (MAF == NULL) {success = false; goto cleanup_MAF;}

	double *S      = (double *) calloc( (unsigned long) nsample          * sizeof( double ), 64 );
	if (S == NULL) {success = false; goto cleanup_S;}

	double *U      = (double *) calloc( (unsigned long) ngene  * nsample * sizeof( double ), 64 );
	if ( U == NULL) {success = false; goto cleanup_U;}

	double *SB2    = (double *) calloc( (unsigned long) nsnp             * sizeof( double ), 64 );
	if (SB2 == NULL) {success = false; goto cleanup_SB2;}

	double *Q      = (double *) calloc( (unsigned long) nsnp             * sizeof( double ), 64 );
	if (Q == NULL) {success = false; goto cleanup_Q;}

	double *P      = (double *) calloc( (unsigned long) nsnp             * sizeof( double ), 64 );
	if (P == NULL) {success = false; goto cleanup_P;}

    double *MUQ    = (double *) calloc( (unsigned long) nsnp             * sizeof( double ), 64 );
	if (MUQ == NULL) {success = false; goto cleanup_MUQ;}

    double *SIGQ   = (double *) calloc( (unsigned long) nsnp             * sizeof( double ), 64 );
	if (SIGQ == NULL) {success = false; goto cleanup_SIGQ;}

	for ( i = 0; i < ngene * nsample; i++ ) {
		GX[i] = randn(0, 1);
		U[i] = 0.0;
	}

	for ( i = 0; i < nsnp; i++ ) {
		MAF[i] = randu (0.1, 0.2);
		SB2[i] = 0.001;
		Q[i] = 0.0;
		P[i] = 0.0;
		MUQ[i] = 0.0;
		SIGQ[i] = 0.0;
	}

	for ( i = 0; i < nsnp; i++ ) {
		for ( j = 0; j < nsample; j++ ) {
			GT[ i*nsample + j ] = randb (MAF[i]);
		}
	}

	for ( i = 0; i < nsample; i++ ) {
		S[i] = 0.0;
	}

	if (DEBUG) {
		printGX (GX, ngene, nsample, 5, 10);
		printGT (GT, MAF, nsnp,  nsample, 5, 30);
	}

	// Check the SVD
	printf ("Checking SVD ...\n");
	success = dsvd ( GX, ngene, nsample, S, U );
	if (DEBUG) { strcpy(mstring, "First 10 elements of S: "); printvec (mstring, S, 10); }

	// Check Q-Score permuted null
	printf ("Checking Q-Score Perm Null ...\n");
	success = qscore (GT, GX, SB2, ngene, nsnp, nsample, 2, MAF, Q, P, MUQ, SIGQ);
	if (DEBUG) { strcpy(mstring, "Q of first 10 SNPs: "); printvec (mstring, Q, 10); }

	// Check Q-Score MAF
	printf ("Checking Q-Score Perm Null ...\n");
	success = qscore (GT, GX, SB2, ngene, nsnp, nsample, 1, MAF, Q, P, MUQ, SIGQ);
	if (DEBUG) { strcpy(mstring, "Q of first 10 SNPs: "); printvec (mstring, Q, 10); }

cleanup:
cleanup_SIGQ:
	free(SIGQ);
cleanup_MUQ:
	free(MUQ);
cleanup_P:
	free(P);
cleanup_Q:
	free(Q);
cleanup_SB2:
	free(SB2);
cleanup_U:
	free(U);
cleanup_S:
	free(S);
cleanup_MAF:
	free(MAF);
cleanup_GT:
	free(GT);
cleanup_GX:
	free(GX);

	if (DEBUG) printf ("Memory free :)\n");

	clock_t end = clock();
        double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf ("Time taken: %f sec.\n", time_spent);

	return success;
}

double randn (double mu, double sigma) {
	double U1, U2, W, mult;
	static double X1, X2;
	static int call = 0;

	if (call == 1) {
		call = !call;
		return (mu + sigma * (double) X2);
	}
	do {
		U1 = randu(-1, 1);
		U2 = randu(-1, 1);
		W = pow (U1, 2) + pow (U2, 2);
	}
	while (W >= 1 || W == 0);
	mult = sqrt ((-2 * log (W)) / W);
	X1 = U1 * mult;
	X2 = U2 * mult;
	call = !call;
	return (mu + sigma * (double) X1);
}

double randb (double f) {
	int i;
	double r, x;
	double f0, f1, f2;
	
	f2 = f * f;
	f1 = 2 * f * (1 - f);
	f0 = (1 - f) * (1 - f);
	r = randu (0, 1);
	if (r < f0) {
		x = 0;
	} else if ( r < (f0 + f1) ) {
		x = 1;
	} else {
		x = 2;
	}
	return x;
}

double randu (double a, double b) {
	double r, x;
	r = ((double) rand() / RAND_MAX);
	x = (r * (b - a)) + a;
	return x;
}

void printGX (double *GX, int ngene, int nsample, int gprint, int sprint) {
	int i, j;
	printf ("Gene Expression (%d x %d)\n", ngene, nsample);
	printf ("======================================================\n");
	for ( i = 0; i < gprint; i++ ) {
		printf ("Gene %d:", i+1);
		for ( j = 0; j < sprint; j++ ) {
			printf(" %6.3f", GX[ i*nsample + j ]);
		}
		printf("\n");
	}
	printf("\n");
	return;
}

void printGT (double *GT, double *MAF, int nsnp, int nsample, int gprint, int sprint) {
	int i, j;
	printf ("Genotype (%d x %d)\n", nsnp, nsample);
	printf ("======================================================\n");
	for ( i = 0; i < gprint; i++ ) {
		printf ("SNP %d (%4.2f):", i+1, MAF[i]);
		for ( j = 0; j < sprint; j++ ) {
			printf(" %g", GT[ i*nsample + j ]);
		}
		printf("\n");
	}
	printf("\n");
	return;
}

void printvec(char* name, double* arr, int n) {
	int i;
	printf(name);
	for ( i = 0; i < n; i++ ) {
		printf (" %g", arr[i]);
	}
	printf("\n");
	printf("\n");
}
