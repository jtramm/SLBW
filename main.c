// 22.211 PSET 1 - Problem 1 - SLBW Resonance Data Processing
// John Tramm

// ENDF Data Source: http://www.nndc.bnl.gov/sigma/getInterpreted.jsp?evalid=15324&mf=2&mt=151
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

typedef struct{
	double E;
	double sigma_g;
	double sigma_n;
	double sigma_t;
} XS;

Resonance * res_read(int * n_resonances);
double complex FNF( double complex Z );

XS calculate_XS( double E, Resonance * R, int nr )
{
	double sigma_pot = 11.2934;
	double k = 8.6173324e-5;
	double temp = 0.00001;
	double A = 238.05078826;
	XS xs = {0};
	xs.E = E;

	for( int j = 0; j < nr; j++ )
	{
		double r = 2603911.0 / R[j].Eo * (A+1) / A;
		double q = 2.0 * sqrt(r * sigma_pot);
		double T = R[j].Tn + R[j].Tg;
		double x = 2.0 * (E - R[j].Eo) / T;
		double xi = T * sqrt(A / (4.0 * k * temp * R[j].Eo));
		double complex faddeeva_in = x + I;
		faddeeva_in *= xi;
		//double complex faddeeva_out = xi * FNF( faddeeva_in);
		double complex faddeeva_out = xi * Faddeeva_w( faddeeva_in, 0.0);
		double psi = sqrt(M_PI) * creal(faddeeva_out); 
		double chi = sqrt(M_PI) * cimag(faddeeva_out);
		xs.sigma_g += R[j].Tn * R[j].Tg / (T*T) * sqrt(R[j].Eo / E) * r * psi;
		xs.sigma_n += R[j].Tn * R[j].Tn / (T*T) * ( r * psi + q * T/R[j].Tn * chi ); 
	}
	xs.sigma_n += sigma_pot;
	xs.sigma_t = xs.sigma_g + xs.sigma_n;

	return xs;
}


int main(void)
{
	int gp = 10000;
	int nr;
	Resonance * R = res_read(&nr);
	XS * xs = (XS *) malloc( gp * sizeof(XS));
	double * E = (double *) malloc( gp * sizeof(XS));

	for( int i = 0; i < gp; i++ )
	{
		double delta = 4.0 / gp;
		E[i] = pow(10.0, -2.0 + delta*i);
	}

	for( int i = 0; i < gp; i++ )
		xs[i] = calculate_XS(E[i], R, nr);

	FILE * fp = fopen("data.dat", "w");
	for( int i = 0; i < gp; i++ )
	{
		fprintf(fp, "%e\t%e\t%e\t%e\n",
				xs[i].E,
				xs[i].sigma_g,
				xs[i].sigma_n,
				xs[i].sigma_t);
	}
	fclose(fp);

	free(E);
	return 0;
}

Resonance * res_read(int * n_resonances)
{
	printf("Reading from \"resonances.dat\" file...\n");
	FILE * fp = fopen("resonances.dat", "r");
	if( fp == NULL )
	{
		printf("ERROR! No resonances file found!\n");
		exit(1);
	}

	int lines = 0;
	char c;
	while ( (c = getc(fp)) != EOF)
	{
		if (c == '\n')
			lines++;
	}
	rewind(fp);

	Resonance * R = (Resonance *) malloc( lines * sizeof(Resonance));

	for( int i = 0; i < lines; i++ )
	{
		float dummy;
		float * dum = &dummy;
		int num = fscanf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
				&R[i].Eo,
				dum,
				&R[i].Tn,
				&R[i].Tg,
				dum,
				dum);
	}

	fclose(fp);
	*n_resonances = lines;

	return R;
}

// My own implementation
// See https://github.com/jtramm/FNF for details
double complex FNF( double complex Z )
{
	// Abrarov 
	if( cabs(Z) < 6.0 )
	{
		// Precomputed parts for speeding things up
		// (N = 10, Tm = 12.0)
		double complex prefactor = 8.124330e+01 * I;
		double an[10] = {
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
		double neg_1n[10] = {
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

		double denominator_left[10] = {
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

		double complex W = I * ( 1 - cexp(I*12.*Z) ) / (12. * Z );
		double complex sum = 0;
		for( int n = 0; n < 10; n++ )
		{
			complex double top = neg_1n[n] * cexp(I*12.*Z) - 1.;
			complex double bot = denominator_left[n] - 144.*Z*Z;
			sum += an[n] * (top/bot);
		}
		W += prefactor * Z  * sum;
		return W;
	}

	// Pre-computed parameters
	double a = 0.512424224754768462984202823134979415014943561548661637413182;
	double b = 0.275255128608410950901357962647054304017026259671664935783653;
	double c = 0.051765358792987823963876628425793170829107067780337219430904;
	double d = 2.724744871391589049098642037352945695982973740328335064216346;

	// Three Term Asymptotic Expansion
	double complex W = I * Z * (a/(Z*Z - b) + c/(Z*Z - d));

	return W;
}
