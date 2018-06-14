#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mkl.h>
#include <mkl_lapacke.h>
#include "chisqr.h"

double ** transpose(double ** mat, int m, int n);
double ** memalloc(double ** mat, int m, int n);
void     mat_square(double ** mat, int m, int n);
void	 arr_square(double * arr, int n);
double ** matmult(double ** a, double ** b, int m, int n, int p);
double    mean(double * values, int n);
double    var(double * values, int n);
double *  mat_var( double ** mat, int m, int n);
double *  mat_mean( double ** mat, int m, int n);
double *  arr_add(double * a, double *b, int n);
double *  arr_add_val(double * a, double val, int n);
double *  arr_div(double * a, double *b, int n);
double    arr_dot(double * a, double *b, int n);
double *  arr_mult(double * a, double b, int n);
double **  mat_sum(double ** a, double ** b, int m, int n);
void 	 free_mat(double ** ptr, int n);
int      min(int a, int b);
int      max(int a, int b);
double 	 fourth_moment(double * a, int n);
double 	 second_moment(double * a, int n);
void     null_maf(double * s2term, double * Rscore, double * maf, double ** Q, int nsamples, int nsnps, double* pvals, double * mu, double * sigma);
void     null_perm(double * gt, double * Rscore, double ** Q, int nsamples, double* pvals, double * mu, double * sigma);
double *  toRowMajor(double ** A, int m, int n);


int dsvd ( double *A, int m, int n, double *S, double *U )
{
    int info;
    int k = min(m, n);
    double superb[k - 1];
    double Vt[k * n];
    info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'S', 'S', m, n, A, n, S, U, k, Vt, n, superb);
    return info;
}		/*  -----  end of function dsvd  ----- */ 


void rscore_maf(double ** gt, double ** gx, double sigmabeta2, double * maf, int ngenes, int nsnps, int nsamples, double * Rscore, double * pvals, double * mu, double * sigma ){
	double * Urowmajor = (double *) malloc (nsamples * ngenes * sizeof(double)); 
	double * S  	  = (double *) malloc (nsamples * sizeof(double));
        double ** gxT = transpose(gx, ngenes, nsamples);
        double *  gxTrowmajor = toRowMajor(gxT, nsamples, ngenes);
        dsvd(gxTrowmajor, nsamples, ngenes, S, Urowmajor);
        double ** U = (double**) calloc (nsamples, sizeof(double *));
        for(int i=0; i < nsamples; i++){
                *(U + i) = Urowmajor + i*nsamples;
        }
	double * s2 	= S;
	arr_square(s2, nsamples);
	double * s2mod 	= (double *) malloc (nsamples * sizeof(double));
	memcpy(s2mod, s2, sizeof(double)*nsamples);
	double sigmax2 = 1;
	double fact 	=  sigmax2 / sigmabeta2;
	double * s2term_deno = arr_add_val(s2, fact, nsamples);
	double * s2term 	= arr_div(s2, s2term_deno, nsamples);
	double ** Q;
	Q = memalloc(Q, nsamples, nsamples);
	for(int j =0; j<nsamples;j++){
                        for(int k=0; j<nsamples;j++){
                                for(int l=0; l< nsamples; l++){
                                        Q[j][k] += U[j][l] * U[j][l]*s2term[l];
                                }
                        }
                }


	double ** gt_U = matmult(gt, transpose(U, nsamples, nsamples), nsnps, nsamples, nsamples);
	mat_square(gt_U, nsnps, nsamples);
	for(int i=0; i< nsnps; i++){
		for(int j=0;j<nsamples;j++){
			Rscore[i] += gt_U[i][j] * s2term[j];
		}
	}
	
	null_maf(s2term, Rscore, maf, Q, nsamples, nsnps, pvals, mu, sigma);
	
	return ;
}

void  null_maf(double * s2term, double * Rscore, double * maf, double ** Q, int nsamples, int nsnps, double* pvals, double * mu, double * sigma){
	*mu = 0;
	*sigma = 0;
	for(int i=0; i<nsamples;i++){
		*mu += s2term[i];
		*sigma += 2*s2term[i]*s2term[i];
	}
	double x4 = fourth_moment(maf, nsnps);
	double q_term = 0;
	for(int i=0;i<nsamples;i++){
		q_term += Q[i][i]*Q[i][i];
	}
	double corr_maf = (x4 - 3)*q_term;
	*sigma += corr_maf;
	double mscale = (*sigma / *mu)/2;
	double df = *mu / mscale;
	double Rscaled[nsnps];
	for(int i=0;i<nsnps;i++){
		Rscaled[i] = Rscore[i]/mscale;
		*(pvals + i) = chisqr((int)df, (double)Rscaled[i]);
	}
	
	
}

void rscore_perm(double ** gt, double ** gx, double * sigmabeta2, int ngenes, int nsnps, int nsamples , double * Rscore, double* pvals, double * mu, double * sigma){
        double * Urowmajor = (double *) malloc (nsamples * ngenes * sizeof(double));
        double * S         = (double *) malloc (nsamples * sizeof(double));
        double ** gxT = transpose(gx, ngenes, nsamples);
        double *  gxTrowmajor = toRowMajor(gxT, nsamples, ngenes);
        dsvd(gxTrowmajor, nsamples, ngenes, S, Urowmajor);
        double ** U = (double**) calloc (nsamples, sizeof(double *));
        for(int i=0; i < nsamples; i++){
                *(U + i) = Urowmajor + i*nsamples;
        }

	double * s2      = S;
        arr_square(s2, nsamples);
        double * s2mod   = (double *) malloc (nsamples * sizeof(double));
        memcpy(s2mod, s2, sizeof(double)*nsamples);
        double * sigmax2        = mat_var(gt, nsnps, nsamples);
        //double * Rscore = (double *)calloc(nsnps , sizeof(double));
	double ** Ut;
        Ut = transpose(U,nsamples, nsamples);
        double ** gt_U = matmult(gt, Ut, nsnps, nsamples, nsamples);
        mat_square(gt_U, nsnps, nsamples);
	for (int i=0; i<nsnps; i++){
		double fact      =  sigmax2[i] / sigmabeta2[i];
        	double * s2term_deno = arr_add_val(s2, fact, nsamples);
        	double * s2term  = arr_div(s2, s2term_deno, nsamples);
		double ** Q;
        	Q = memalloc(Q, nsamples, nsamples);


		for(int j =0; j<nsamples;j++){
			for(int k=0; j<nsamples;j++){
				for(int l=0; l< nsamples; l++){
					Q[j][k] += U[j][l]*U[j][l]*s2term[l];
				}
			}
		}
		for(int j=0;j<nsamples;j++){
			Rscore[i] += gt_U[i][j] * s2term[j];
                }
		null_perm(*(gt+i),(Rscore + i), Q, nsamples, pvals, (mu + i), (sigma + i));
		free_mat(Q, nsamples);
		free(s2term);
		free(s2term_deno);
	}
	free_mat(gt_U, nsnps);
	return ;
}
void  null_perm(double* gt, double * Rscore, double ** Q, int nsamples, double* pvals, double * mu, double * sigma){
	int N 		= nsamples;
	double mu2 	= second_moment(gt, nsamples);
	double mu4	= fourth_moment(gt, nsamples);
	double q11 = 0, q2 = 0, q31 = 0, q4 = 0, q22 = 0, q211 = 0, v31 = 0, v22 = 0, v211 = 0, v1111 = 0;	
	double * q31_left = (double * )calloc(nsamples, sizeof(double));
	double * q31_right = (double * )calloc(nsamples, sizeof(double));
	for(int i=0; i < nsamples; i++){
		q2 += Q[i][i];
		q4 += Q[i][i]*Q[i][i];
		double q211_inner = 0;
		for(int j=0;j< nsamples;j++){
			q11 += Q[i][j];
			q22 += Q[i][j]*Q[i][j];
			q211_inner = Q[i][j];
		}
		q31_right[i] = q211_inner;
		q31_left [i] = Q[i][i];
		q211 += q211_inner*q211_inner;
	}
	q31 	= arr_dot(q31_left,q31_right,nsamples);
	
	v31 	= -(mu4/(N-1));
	v22 	= (N*mu2*mu2 - mu4)/(N-1);
	v211 	= -(v31 + v22)/(N-2);
	v1111	= -3 * v211/(N-3);
	
	*sigma 	= v1111*(q11*q11 - 2*q2*q11 - 4*q211 + 8*q31 + 2*q22 + q2*q2 - 6*q4) + 2*v211*(q2*q11 + 2*q211 - 6*q31 - 2*q22 - q2*q2 + 6*q4) + v22*(q2*q2 + 2*q22 - 3*q4) + 4*v31*(q31 - q4) + mu4*q4; 
	
	*mu = (mu2/(N-1))*(q2 - q11);
	*sigma	=  *sigma - (*mu)*(*mu);
	double mscale = (*sigma / *mu)/2;
        double df = *mu / mscale;
        double Rscaled;
	Rscaled = *Rscore/mscale;
        *pvals = chisqr((int)df, (double)Rscaled);
	
}

int driver(double * geno, double * expr, double * sigmabeta2, int null, double * maf,  int ngenes, int nsnps, int nsamples, double * Rscore, double *pvals, double* mu, double* sigma){
	double ** gt = (double**)calloc(nsnps,sizeof(double *));
	double  ** gx = (double**)calloc(ngenes,sizeof(double *));
	for(int i=0;i< nsnps;i++){
		*(gt + i) = geno + i*nsamples;
	}	
	for(int i=0;i<ngenes;i++){
		*(gx + i) = expr + i*nsamples;
	}
	if(null == 1){
		rscore_maf(gt, gx, sigmabeta2[0], maf, ngenes, nsnps, nsamples, Rscore, pvals,mu,sigma);
	}
	else if (null == 2){
		rscore_perm(gt, gx, sigmabeta2, ngenes, nsnps, nsamples, Rscore, pvals, mu, sigma);
	}
	return 0;
}

int main(){//double ** gt, double ** gx, double sigmabeta2, double * maf, int ngenes, int nsnps, int nsamples, double * Rscore, double * pvals; double * mu, double * sigma
	int nsnps, ngenes, nsamples;
	ngenes = 5; //2300;
	nsnps = 10; //50;
	nsamples = 3; //338;
	double ** gt;
	double ** gx;
	double * sigmabeta2 = (double *)calloc(nsnps, sizeof(double));
	double * maf = (double *)calloc(nsnps, sizeof(double));
	gx = memalloc(gx, ngenes, nsamples);   //i = 5, n = 3, g = 5
        gt = memalloc(gt, nsnps, nsamples);
	for(int i=0; i < ngenes; i++){
		for(int j=0;j<nsamples;j++){
                gx[i][j] = i + j + 1;
		}
        }
        for(int i=0; i < nsnps; i++){
		for(int j=0; j<nsamples;j++){
                gt[i][j] = i + j + 1;
		}
                sigmabeta2[i] = 0.0016;
		maf[i] = 0.3;
        }
	clock_t start, end;
        double cpu_time_used;
	start = clock();
	double *Rscore, *mu, *sigma;
	double *pvals;
	
	Rscore = (double *)calloc(nsnps , sizeof(double));
        pvals = (double *)calloc(nsnps , sizeof(double));
        mu = (double *)calloc(nsnps, sizeof(double));
        sigma = (double *)calloc(nsnps, sizeof(double));
	
	rscore_perm(gt, gx, sigmabeta2, ngenes, nsnps, nsamples, Rscore, pvals, mu, sigma);
	//rscore_maf(gt, gx, 0.005, maf, ngenes, nsnps, nsamples, Rscore, pvals, mu, sigma);
	end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	for(int i=0; i<nsnps;i++){
		printf("%lf ", Rscore[i]);
	}
	
	printf("mu %f\n sigma %f\n %lf", *mu, *sigma, cpu_time_used);
	return 0;
}


void free_mat(double ** mat, int n){
	for(int i=0; i<n;i++){
		free(*(mat + i));
	}
	free(mat);
}

double ** mat_sum(double** a, double **b, int m, int n){
	double ** sum;
	sum = memalloc(sum, m, n);
	for(int i = 0; i< m; i++){
		for(int j=0;j< n; j++){
			*(*(sum + i) + j) = a[i][j] + b[i][j];
		}
	}
	return sum;
}

int min(int a, int b){
	if(a < b) return a;
	else return b;
}

int max(int a, int b){
	if(a > b) return a;
	else return b;
}

double ** transpose(double ** mat, int m, int n){
	double ** tr = (double **)malloc(n * sizeof(double *));
	for(int i = 0; i < n ; i ++){
		*(tr + i) = (double *)malloc(m * sizeof(double));
		for (int j = 0; j < m; j ++){
			*(*(tr + i) + j) = mat[j][i];
		}
	}

	return tr;
}

double ** memalloc(double ** mat, int m, int n){
	mat = (double **)calloc(m, sizeof(double *));                                              
        for(int i = 0; i < m ; i ++){                                                                     
                *(mat + i) = (double *)calloc(n, sizeof(double));                                                                                                                                    
        }
	return mat;
}

void mat_square(double ** mat, int m, int n){
	for(int i=0; i < m; i++){
		for(int j=0; j < n; j++){
			*(*(mat + i) + j) = mat[i][j] * mat[i][j];			
		}
      	}
	return;
}

void arr_square(double * arr, int n){
	for( int j =0; j < n; j++) {
		*(arr + j) = arr[j] * arr[j];
	}
	return;
}

double arr_dot(double *a, double *b, int n){
	double dot = 0;
        for (int i = 0; i< n; i++){
                dot += a[i] * b[i];
        }
        return dot;
}

double ** matmult(double ** a, double ** b, int m, int n, int p){
	double **result = (double **)malloc (m * sizeof (double *));
	// if (!result) throw error

	register int i=0, j=0, k=0;
	for (i = 0; i < m; i++)
	{
        /* calloc initializes all to '0' */
        result[i] = (double *)calloc (p, sizeof (double));
        // if (!result[i]) throw error
	}

	for (i = 0; i < m; i++)
	{
        	for (j = 0; j < p; j++)
        	{
            		for (k = 0; k < n; k++)
            		{
                	result [i][j] += a [i][k] * b [k][j];
            	}
        }
    }

    return result;

}

/**
 * calculate the mean of the provided list values, containing n values.
 */
double mean(double * values, int n)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
	{
        	sum += values[i];
    	}
   	 return sum / n;
}

double second_moment(double * values, int n)
{
        double sum = 0;
        for (int i = 0; i < n; i++)
        {
                sum += values[i]*values[i];
        }
         return sum / n;
}

double fourth_moment(double * values, int n)
{
        double sum = 0;
        for (int i = 0; i < n; i++)
        {
                double val2 = values[i]*values[i];
		sum += val2*val2;
        }
         return sum / n;
}



/**
 * calculate the variance of the provided list values, containing n values.
 */

double var(double * values, int n)
{
    	double valuesMean = mean(values, n);
    	double sum = 0;
    	for (int i = 0; i < n; i++)
    	{
        	sum += (values[i] - valuesMean) * (values[i] - valuesMean);
    	}
    	return sum / n;
}

double* mat_var(double ** mat, int m ,int n){
	double * vars = (double *)malloc(m * sizeof(double));
	for (int i = 0; i< m; i++){
		vars[i] = var(*(mat + i), n);
	}
	return vars;
} 

double* mat_mean(double ** mat, int m , int n){
        double * means = (double *)malloc(m * sizeof(double));
        for (int i = 0; i< m; i++){
                means[i] = mean(*(mat + i), n);
        }
        return means;
}

double* arr_add(double * a, double * b, int n){
	double * sum = (double *)malloc(n * sizeof(double));
        for (int i = 0; i< n; i++){
                sum[i] = a[i] + b[i];
        }
        return sum;

}

double* arr_add_val(double * a, double val, int n){
        double * sum = (double *)malloc(n * sizeof(double));
        for (int i = 0; i< n; i++){
                sum[i] = a[i] + val;
        }
        return sum;

}


double* arr_div(double * a, double * b, int n){
        double * div = (double *)malloc(n * sizeof(double));
        for (int i = 0; i< n; i++){
                div[i] = a[i] / b[i];
        }
        return div;

}

double* arr_mult(double * a, double  b, int n){
        double * mul = (double *)malloc(n * sizeof(double));
        for (int i = 0; i< n; i++){
                mul[i] = a[i] * b;
        }
        return mul;

}

double * toRowMajor(double ** A, int m, int n) {
        double * AT = (double *)calloc(m * n, sizeof(double));
        for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                        AT[i * n + j] = A[i][j];
                }
        }
        return AT;
}
