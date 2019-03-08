#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double linear_regression(double* X, double* Y, int N) {
    int     i;
    double  sum_x   = 0;
    double  sum_y   = 0;
    double  sum_xy  = 0;
    double  sum_xsq = 0;
    double  sum_ysq = 0;
    double  NCxx, NCxy;
    double  b0, b1;
    double  sum_sq_total, sum_sq_residual, sum_sq_explained;
    double  mean_sq_residual, mean_sq_explained;
    double  fstat;
    int     dfres, dfexp;

    dfres = N - 2; // residual degrees of freedom
    dfexp = 2 - 1; // model    degrees of freedom

    for(i = 0; i < N; i = i + 1) {
        sum_x   += X[i];
        sum_y   += Y[i];
        sum_xy  += X[i]*Y[i];
        sum_xsq += X[i]*X[i];
        sum_ysq += Y[i]*Y[i];
    }

    NCxx = sum_xsq - (sum_x * sum_x) / N;
    NCxy = sum_xy  - (sum_x * sum_y) / N;
    b1 = NCxy / NCxx;
    b0 = (sum_y - (b1 * sum_x)) / N;
    sum_sq_total      = sum_ysq - sum_y*sum_y / N;
    sum_sq_residual   = sum_ysq + b1*b1*sum_xsq + b0*b0*N - 2.0*b0*sum_y - 2.0*b1*sum_xy + 2.0*b0*b1*sum_x;
    sum_sq_explained  = sum_sq_total - sum_sq_residual;
    mean_sq_residual  = sum_sq_residual  / dfres;
    mean_sq_explained = sum_sq_explained / dfexp;
    fstat = mean_sq_explained / mean_sq_residual;
    return fstat;
}


double linear_regression_zstat(double* X, double* Y, int N) {
    int     i;
    double  sum_x  = 0;
    double  sum_y  = 0;
    double  sum_xy = 0;
    double  sum_xx = 0;
    double  sum_yy = 0;
    double  NCxx, NCxy;
    double  b0, b1;
    double  sum_sq_yres, sum_sq_xres;
    double  sesq, se;
    double  zstat;

    for(i = 0; i < N; i = i + 1) {
        sum_x  += X[i];
        sum_y  += Y[i];
        sum_xy += X[i]*Y[i];
        sum_xx += X[i]*X[i];
        sum_yy += Y[i]*Y[i];
    }

    NCxx = sum_xx - (sum_x * sum_x) / N;
    NCxy = sum_xy - (sum_x * sum_y) / N;
    b1 = NCxy / NCxx;
    b0 = (sum_y - (b1 * sum_x)) / N;
    sum_sq_yres = sum_yy - sum_y*sum_y / N;
    sum_sq_xres = sum_xx - sum_x*sum_x / N;
    sesq = (sum_sq_yres - b1 * b1 * sum_sq_xres) / (sum_sq_xres * (N - 2));
    se = sqrt(sesq);
    zstat = (b1 - b0) / se;
    return zstat;
}


int fit(double* genotype, double* expression, int nsnps, int ngene, int nsample, double* farr) {

    int i, j, k;
    int pos;
    double fval;
    //double x[nsample];
    //double y[nsample];
    double *x;
    double *y;

    x = (double*) malloc(nsample * sizeof(double));
    y = (double*) malloc(nsample * sizeof(double));
    if ( x == NULL || y == NULL) {
        printf ( "\ndynamic memory allocation failed\n" );
        free ( x );
        free ( y );
        exit ( 0 );
    }



    /* Perform operations using ndonor elements from offset  */
    fval = 0;

    for (i = 0 ; i < nsnps; i++) {
        // Get the x for regression (ndonor elements)
        for (j = 0; j < nsample; j++){
            x[j] = genotype[i * nsample + j];
        }

        // For G genes, calculate F-statistic
        for (k = 0; k < ngene; k++) {
            // Get the y for regression (ndonor elements)
            for (j = 0; j < nsample; j++){
                y[j] = expression[k * nsample + j];
            }
            /* Calculate whatever you need to with X and Y.
               X is the genotype
               Y is the expression
             */
            fval = linear_regression_zstat(x, y, nsample);
            pos = (i * ngene) + k;
            farr[pos] = fval;
        }
    }

    free ( x );
    free ( y );
    x = NULL;
    y = NULL;
    return 0;
}
