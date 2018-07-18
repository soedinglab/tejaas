#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "linear_regression.h"

#define DEBUG true

double randn (double mu, double sigma);
double randb (double f);
double randu (double a, double b);
void   printGX (double *GX, int ngene, int nsample, int gprint, int sprint);
void   printGT (double *GT, double* MAF, int nsnp,  int nsample, int gprint, int sprint);
void   printFSTAT (double *A, int nsnp, int ngene, int sprint, int gprint);

bool main() {
	bool success;
	int ngene;
	int nsample;
	int nsnp;
	int i, j;
	int isucc;
	
	success = false;
	ngene = 219;
	nsample = 80;
	nsnp = 19;
	
	double *GX    = (double *) calloc( (unsigned long) ngene   * nsample * sizeof( double ), 64 );
	if (GX == NULL) {success = false; goto cleanup_GX;}
	
	double *GT    = (double *) calloc( (unsigned long) nsnp    * nsample * sizeof( double ), 64 );
	if (GT == NULL) {success = false; goto cleanup_GT;}
	
	double *MAF    = (double *) calloc( (unsigned long) nsnp             * sizeof( double ), 64 );
	if (MAF == NULL) {success = false; goto cleanup_MAF;}

	double *GTS    = (double *) calloc( (unsigned long) nsample          * sizeof( double ), 64 );
	if (GTS == NULL) {success = false; goto cleanup_GTS;}

	double *FSTAT = (double *) calloc( (unsigned long) nsnp    * ngene   * sizeof( double ), 64 );
	if ( FSTAT == NULL) {success = false; goto cleanup_FSTAT;}

	for ( i = 0; i < ngene * nsample; i++ ) {
		GX[i] = randn(0, 1);
	}

	for ( i = 0; i < nsnp; i++ ) {
		MAF[i] = randu (0.1, 0.2);
	}

	for ( i = 0; i < nsnp; i++ ) {
		for ( j = 0; j < nsample; j++ ) {
			GT[ i*nsample + j ] = randb (MAF[i]);
		}
	}

	if (DEBUG) {
		printGX (GX, ngene, nsample, 5, 10);
		printGT (GT, MAF, nsnp,  nsample, 5, 30);
	}

	for ( i = 0; i < nsnp * ngene; i++ ) {
		FSTAT[i] = 0.0;
	}

	isucc = fit(GT, GX, nsnp, ngene, nsample, FSTAT);

	if (DEBUG) printFSTAT (FSTAT, nsnp, ngene, 10, 20);

cleanup:
cleanup_FSTAT:
	free(FSTAT);
cleanup_GTS:
	free(GTS);
cleanup_MAF:
	free(MAF);
cleanup_GT:
	free(GT);
cleanup_GX:
	free(GX);

	if (DEBUG) printf ("Memory free :)\n");

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
	printf ("Gene Expression (%d genes, %d samples)\n", ngene, nsample);
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
	printf ("Genotype (%d SNPs, %d samples)\n", nsnp, nsample);
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

void printFSTAT (double *A, int nsnp, int ngene, int sprint, int gprint) {
	int i, j;
	printf ("F-statistic matrix from C (%d SNPs, %d genes):\n", nsnp, ngene);
	printf ("======================================================\n");
	for ( i = 0; i < sprint; i++ ) {
		printf ("SNP %2d:", i+1);
		for ( j = 0; j < gprint; j++ ) {
			printf (" %6.3f", A[ i*ngene + j ]);
		}
		printf("\n");
	}
	printf("\n");
	return;
}
