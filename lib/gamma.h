

/*
    Implementation of the Gamma function using Spouge's Approximation in C.

    Written By Jacob
	7/31/2012
    Public Domain

    This code may be used by anyone for any reason
    with no restrictions absolutely free of cost.
*/

#include <math.h>
//#include "gamma.h"


#define ACCURACY 15 // 15
/*
    A is the level of accuracy you wish to calculate.
    Spouge's Approximation is slightly tricky, as you
    can only reach the desired level of precision, if
    you have EXTRA precision available so that it can
    build up to the desired level.

    If you're using double (64 bit wide datatype), you
    will need to set A to 11, as well as remember to
    change the math functions to the regular
    (i.e. pow() instead of powl())

    double A = 11
    long double A = 15

   !! IF YOU GO OVER OR UNDER THESE VALUES YOU WILL !!!
              !!! LOSE PRECISION !!!
*/


double gamma(double N)
{
    /*
        The constant SQRT2PI is defined as sqrt(2.0 * PI);
        For speed the constant is already defined in decimal
        form.  However, if you wish to ensure that you achieve
        maximum precision on your own machine, you can calculate
        it yourself using (sqrt(atan(1.0) * 8.0))
    */

	//const long double SQRT2PI = sqrtl(atanl(1.0) * 8.0);
    const long double SQRT2PI = 2.5066282746310005024157652848110452530069867406099383;

    long double Z = (long double)N;
    long double Sc = powl((Z + ACCURACY), (Z + 0.5));
	Sc *= expl(-1.0 * (Z + ACCURACY));
    Sc /= Z;

	long double F = 1.0;
	long double Ck;
    long double Sum = SQRT2PI;
    int K;


	for(K = 1; K < ACCURACY; K++)
	{
	    Z++;
		Ck = powl(ACCURACY - K, K - 0.5);
		Ck *= expl(ACCURACY  - K);
		Ck /= F;

		Sum += (Ck / Z);

		F *= (-1.0 * K);
	}

	return (double)(Sum * Sc);
}

long double log_gamma(double N)
{
    /*
        The constant SQRT2PI is defined as sqrt(2.0 * PI);
        For speed the constant is already defined in decimal
        form.  However, if you wish to ensure that you achieve
        maximum precision on your own machine, you can calculate
        it yourself using (sqrt(atan(1.0) * 8.0))
    */

	//const long double SQRT2PI = sqrtl(atanl(1.0) * 8.0);
    const long double SQRT2PI = 2.5066282746310005024157652848110452530069867406099383;

    long double Z = (long double)N;
    long double Sc;

    Sc = (logl(Z + ACCURACY) * (Z + 0.5)) - (Z + ACCURACY ) - logl(Z);

	long double F = 1.0;
	long double Ck;
	long double Sum = SQRT2PI;
	int K;


	for(K = 1; K < ACCURACY; K++)
	{
	    Z++;
		Ck = powl(ACCURACY - K, K - 0.5);
		Ck *= expl(ACCURACY - K);
		Ck /= F;

		Sum += (Ck / Z);

		F *= (-1.0 * K);
	}

	return logl(Sum) + Sc;
}

double approx_gamma(double Z)
{
    const double RECIP_E = 0.36787944117144232159552377016147;  // RECIP_E = (E^-1) = (1.0 / E)
    const double TWOPI = 6.283185307179586476925286766559;  // TWOPI = 2.0 * PI

    double D = 1.0 / (10.0 * Z);
    D = 1.0 / ((12 * Z) - D);
    D = (D + Z) * RECIP_E;
    D = pow(D, Z);
    D *= sqrt(TWOPI / Z);

    return D;
}

double approx_log_gamma(double N)
{
    const double LOGPIHALF = 0.24857493634706692717563414414545; // LOGPIHALF = (log10(PI) / 2.0)

    double D;

    D = 1.0 + (2.0 * N);
    D *= 4.0 * N;
    D += 1.0;
    D *= N;
    D = log10(D) * (1.0 / 6.0);
    D += N + (LOGPIHALF);
    D = (N * log(N)) - D;
    return D;

}

/*
// Slightly faster
double approx_gamma(double Z)
{
    const double RECIP_E = 0.36787944117144232159552377016147;  // RECIP_E = (E^-1) = (1.0 / E)
    const double TWOPI = 6.283185307179586476925286766559;  // TWOPI = 2.0 * PI
    const double RECIP_Z = (1.0 / Z);

    double D = (0.1 * RECIP_Z);
    D = 1.0 / ((12 * Z) - D);
    D = (D + Z) * RECIP_E;
    D = pow(D, Z);
    D *= sqrt(TWOPI * RECIP_Z);

    return D;
}
*/

