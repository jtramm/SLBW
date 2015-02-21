// Data Source: http://www.nndc.bnl.gov/sigma/getInterpreted.jsp?evalid=15324&mf=2&mt=151
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

typedef struct{
	float Eo;
	float Tn;
	float Tg;
} Resonance;

int main(void)
{
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

	// sigma_pot = 4 * pi * a^2 (where a is scattering radius)
	float sigma_pot = 11.29; // barns

	float * E = (float *) calloc( 4 * n_gridpoints, sizeof(float));
	float * sigma_f = E + n_gridpoints; 
	float * sigma_n = E + 2 * n_gridpoints;
	float * sigma_t = E + 3 * n_gridpoints;

	// Initialize E to log scale (10^-2 <-> 10^2)
	for( int i = 0; i < n_gridpoints; i++ )
	{
		float delta = 4.f / n_gridpoints;
		E[i] = powf(10.f,-2.f + delta*i);
	}

	// Calculate sigmas
	for( int i = 0; i < n_gridpoints; i++ )
	{
		float r = 2603911 / E[i] * 239 / 238;
		float q = 2.f * sqrtf(r * sigma_pot);

		// Accumulate Contributions from Each Resonance
		for( int j = 0; j < 3; j++ )
		{
			float T = R[j].Tn + R[j].Tg;
			float x = 2.f * (E[i] - R[j].Eo) / T;
			float psi = 1 / (1 + x*x);
			float chi = x / (1 + x*x);
			sigma_f[i] += R[j].Tn * R[j].Tg / (T*T) * sqrtf(R[j].Eo / E[i]) * r * psi;
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
