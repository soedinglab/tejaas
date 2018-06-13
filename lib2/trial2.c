/*
 * =====================================================================================
 *
 *       Filename:  trial2.c
 *
 *    Description:  Trial
 *
 *        Version:  1.0
 *        Created:  13.06.2018 16:45:28
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Saikat Banerjee (banskt), bnrj.saikat@gmail.com
 *   Organization:  Max Planck Institute for Biophysical Chemistry
 *
 * =====================================================================================
 */

#include <stdio.h>

void qscore(double* a) {
	*a = (*a) * (*a);
	return;
}

int main() {
	double Q[3];
	Q[0] = 7;
	Q[1] = 9;
	Q[2] = 3;
	printf ("Value of Q is %g\n", Q[1]);
	qscore(&Q[1]);
	printf ("New value of Q is %g\n", Q[1]);
	
}
