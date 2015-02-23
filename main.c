// Data Source: http://www.nndc.bnl.gov/sigma/getInterpreted.jsp?evalid=15324&mf=2&mt=151
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include "Faddeeva.h"

#define N 10
#define Tm 12.0
#define Tm2 144.0

typedef struct{
	double Eo;
	double Tn;
	double Tg;
} Resonance;

// This one works!
double complex Abrarov_W( double complex Z )
{
	// Precomputed parts for speeding things up
	// (N = 10, Tm = 12.0)
	double complex prefactor = 8.124330e+01 * I;
	double an[N] = {
		2.758402e-01,
		2.245740e-01,
		1.594149e-01,
		9.866577e-02,
		5.324414e-02,
		2.505215e-02,
		1.027747e-02,
		3.676164e-03,
		1.146494e-03,
		3.117570e-04
	};
	double neg_1n[N] = {
		-1.0,
		1.0,
		-1.0,
		1.0,
		-1.0,
		1.0,
		-1.0,
		1.0,
		-1.0,
		1.0
	};

	double denominator_left[N] = {
		9.869604e+00,
		3.947842e+01,
		8.882644e+01,
		1.579137e+02,
		2.467401e+02,
		3.553058e+02,
		4.836106e+02,
		6.316547e+02,
		7.994380e+02,
		9.869604e+02
	};

	double complex W = I * ( 1 - cexp(I*Tm*Z) ) / (Tm * Z );
	double complex sum = 0;
	for( int n = 1; n <= N; n++ )
	{
		int idx = n-1;
		complex double top = neg_1n[idx] * cexp(I*Tm*Z) - 1.0;
		complex double bot = denominator_left[idx] - Tm2*Z*Z;
		sum += an[idx] * (top/bot);
	}
	W += prefactor * Z  * sum;
	return W;
}

double complex slow_Abrarov_W( double complex Z )
{
	double complex W = I * ( 1 - cexp(I*Tm*Z) ) / (Tm * Z );
	double complex sum = 0;
	for( int n = 1; n <= N; n++ )
	{
		double an = 2.0 * sqrt(M_PI) / Tm * exp( - n*n * M_PI * M_PI / (Tm * Tm));
		complex double top = pow(-1,n) * cexp(I*Tm*Z) - 1.0;
		complex double bot = n*n * M_PI*M_PI - Tm*Tm*Z*Z;
		sum += an * (top/bot);
	}
	W += I * Tm*Tm * Z / sqrt(M_PI) * sum;
	return W;
}

int main(void)
{
	int temperature_dependent = 1;

	int n_gridpoints = 1000;

	Resonance R[3];

	// First Resonance
	R[0].Eo = 6.673491e+0;
	R[0].Tn = 1.475792e-3;
	R[0].Tg = 2.300000e-2;

	// Second Resonance
	R[1].Eo = 2.087152e+1;
	R[1].Tn = 1.009376e-2;
	R[1].Tg = 2.286379e-2;

	// Third Resonance
	R[2].Eo = 3.668212e+1;
	R[2].Tn = 3.354568e-2;
	R[2].Tg = 2.300225e-2;

	// 4 pi a^2
	double sigma_pot = 11.29; // barns

	double k = 8.6173324e-5;
	double temp = 1000;

	double * E = (double *) calloc( 4 * n_gridpoints, sizeof(double));
	double * sigma_f = E + n_gridpoints; 
	double * sigma_n = E + 2 * n_gridpoints;
	double * sigma_t = E + 3 * n_gridpoints;

	// Initialize E to log scale (10^-2 <-> 10^2)
	for( int i = 0; i < n_gridpoints; i++ )
	{
		double delta = 4.0 / n_gridpoints;
		E[i] = pow(10.0,-2.0 + delta*i);
	}

	// Calculate sigmas
	for( int i = 0; i < n_gridpoints; i++ )
	{
		double r = 2603911.0 / E[i] * 239.0 / 238.0;
		double q = 2.0 * sqrt(r * sigma_pot);

		// Accumulate Contributions from Each Resonance
		for( int j = 0; j < 3; j++ )
		{
			double T = R[j].Tn + R[j].Tg;
			double x = 2.0 * (E[i] - R[j].Eo) / T;
			double psi, chi;
			if( temperature_dependent )
			{
				double xi = T * sqrt(238.0 / (4.0 * k * temp * R[j].Eo));
				double complex faddeeva_output = xi * Abrarov_W( (x + I) * xi);
				psi = sqrt(M_PI) * creal(faddeeva_output); 
				chi = sqrt(M_PI) * cimag(faddeeva_output);
			}
			else
			{
				psi = 1.0 / (1.0 + x*x);
				chi = x / (1.0 + x*x);
			}
			sigma_f[i] += R[j].Tn * R[j].Tg / (T*T) * sqrt(R[j].Eo / E[i]) * r * psi;
			sigma_n[i] += R[j].Tn * R[j].Tg / (T*T) * ( r * psi + q * chi ) + sigma_pot; 
		}

		// Calculate Total XS
		sigma_t[i] += sigma_f[i] + sigma_n[i];

	}


	// Save Data to File
	FILE * fp = fopen("data.dat", "w");
	for( int i = 0; i < n_gridpoints; i++ )
	{
		fprintf(fp, "%lf\t%lf\t%lf\t%lf\n",
				E[i],
				sigma_f[i],
				sigma_n[i],
				sigma_t[i]);
	}
	fclose(fp);

	free(E);
	return 0;
}
