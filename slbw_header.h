
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
XS calculate_XS( double E, double temp, Resonance * R, int nr );
void res_out( XS * xs, int gp );
void graph_driver(void);
void RI_driver(void);
void NR_WR_Driver(void);
void find_RI( double e1, double e2, int gp, double temp );
void find_NR_RI( double e1, double e2, int gp, double temp, double s_b );
void find_WR_RI( double e1, double e2, int gp, double temp, double s_b );
double phi_RI( double E, double temp, Resonance * R, int nr );
double integral_RI( double E, double temp, Resonance * R, int nr );
double phi_RI_Narrow( double E, double temp, Resonance * R, int nr, double s_b );
double integral_RI_Narrow( double E, double temp, Resonance * R, int nr, double s_b );
double phi_RI_Wide( double E, double temp, Resonance * R, int nr, double s_b );
double integral_RI_Wide( double E, double temp, Resonance * R, int nr, double s_b );
