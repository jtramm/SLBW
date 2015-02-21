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

	double * E = (double *) calloc( 4 * n_gridpoints, sizeof(double));
	double * sigma_f = E + n_gridpoints; 
	double * sigma_n = E + 2 * n_gridpoints;
	double * sigma_t = E + 3 * n_gridpoints;

	// Initialize E to log scale (10^-2 <-> 10^2)
	for( int i = 0; i < n_gridpoints; i++ )
	{
		double delta = 4.f / n_gridpoints;
		E[i] = pow(10.f,-2.f + delta*i);
	}

	// Calculate sigmas
	for( int i = 0; i < n_gridpoints; i++ )
	{
		double r = 2603911 / E[i] * 239 / 238;
		double q = 2.f * sqrt(r * sigma_pot);

		// Accumulate Contributions from Each Resonance
		for( int j = 0; j < 3; j++ )
		{
			double T = R[j].Tn + R[j].Tg;
			double x = 2.f * (E[i] - R[j].Eo) / T;
			if( temperature_dependent )
			{
				double k = 8.6173324e-5;
				double temp = 1000;
				double xi = T * sqrt(238.f / (4.f * k * temp * R[j].Eo));
				double psi = sqrt(M_PI) * creal( xi * Faddeeva_w( (x + I) * xi ), 0.f ); 
			}
			else
			{
				double psi = 1 / (1 + x*x);
				double chi = x / (1 + x*x);
				sigma_f[i] += R[j].Tn * R[j].Tg / (T*T) * sqrt(R[j].Eo / E[i]) * r * psi;
				sigma_n[i] += R[j].Tn * R[j].Tg / (T*T) * ( r * psi + q * chi ) + sigma_pot; 
			}
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
