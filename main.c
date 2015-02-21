// Data Source: http://www.nndc.bnl.gov/sigma/getInterpreted.jsp?evalid=15324&mf=2&mt=151
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include "Faddeeva.h"

typedef struct{
	double Eo;
	double Tn;
	double Tg;
} Resonance;

int main(void)
{
	int temperature_dependent = 0;

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
				double complex faddeeva_output = xi * Faddeeva_w( (x + I) * xi, 0.0);
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
