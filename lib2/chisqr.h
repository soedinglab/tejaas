
/*
    Implementation of the Chi-Square Distribution &
	Incomplete Gamma Function in C

	Written By Jacob Wells
	July 31, 2012
    Based on the formulas found here:

    Wikipedia - Incomplete Gamma Function -> Evaluation formulae -> Connection with Kummer's confluent hypergeometric function
    http://en.wikipedia.org/wiki/Regularized_Gamma_function#Connection_with_Kummer.27s_confluent_hypergeometric_function

    Wikipedia - Chi-squared Distribution -> Cumulative distribution function
    http://en.wikipedia.org/wiki/Chi-squared_distribution#Cumulative_distribution_function

    These functions are placed in the Public Domain, and may be used by anyone, anywhere, for any reason, absolutely free of charge.

*/


#include <stdio.h>
#include <math.h>
//#include "chisqr.h"
#include "gamma.h"

static double igf(double S, double Z);
static long double log_igf(long double S, long double Z);
static long double KM(long double S, long double Z);


double chisqr(int Dof, double Cv)
{
    //printf("Dof:  %i\n", Dof);
    //printf("Cv:  %f\n", Cv);
    if(Cv < 0 || Dof < 1)
    {
        return 0.0;
    }
	double K = ((double)Dof) * 0.5;
	double X = Cv * 0.5;
	if(Dof == 2)
	{
		return exp(-1.0 * X);
	}
	long double PValue, Gam;
    long double ln_PV;
    ln_PV = log_igf(K, X);

    Gam = approx_gamma(K);
    //Gam = lgammal(K);
    //Gam = log_gamma(K);

    ln_PV -= Gam;
    PValue = 1.0 - expl(ln_PV);

	return (double)PValue;

}


/*
	Returns the Natural Logarithm of the Incomplete Gamma Function
	
	I converted the ChiSqr to work with Logarithms, and only calculate 
	the finised Value right at the end.  This allows us much more accurate
	calculations.  One result of this is that I had to increase the Number 
	of Iterations from 200 to 1000.  Feel free to play around with this if 
	you like, but this is the only way I've gotten to work.  
	Also, to make the code easier to work it, I separated out the main loop.  
*/

static long double log_igf(long double S, long double Z)
{
	if(Z < 0.0)
	{
		return 0.0;
	}
	long double Sc, K;
	Sc = (logl(Z) * S) - Z - logl(S);

    K = KM(S, Z);

    return logl(K) + Sc;
}

static long double KM(long double S, long double Z)
{
	long double Sum = 1.0;
	long double Nom = 1.0;
	long double Denom = 1.0;

	for(int I = 0; I < 1000; I++) // Loops for 1000 iterations
	{
		Nom *= Z;
		S++;
		Denom *= S;
		Sum += (Nom / Denom);
	}

    return Sum;
}


/*
	Incomplete Gamma Function

	No longer need as I'm now using the log_igf(), but I'll leave this here anyway.
*/

static double igf(double S, double Z)
{
	if(Z < 0.0)
	{
		return 0.0;
	}
	long double Sc = (1.0 / S);
	Sc *= powl(Z, S);
	Sc *= expl(-Z);

	long double Sum = 1.0;
	long double Nom = 1.0;
	long double Denom = 1.0;

	for(int I = 0; I < 200; I++) // 200
	{
		Nom *= Z;
		S++;
		Denom *= S;
		Sum += (Nom / Denom);
	}

	return Sum * Sc;
}





