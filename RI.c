#include "slbw_header.h"


void find_WR_RI( double e1, double e2, int gp, double temp, double s_b )
{
	int nr;
	Resonance * R = res_read(&nr);
	nr = 3;
	double s_p = 11.2934;

	double range = e2 - e1;
	double del = range / gp;
	double RI = 0;

	for( int i = 0; i < gp; i++ )
	{
		double low = e1 + del*i;
		double high = e1 + del*(i+1);
		double mid = (low+high)/2;

		XS A = calculate_XS(low, temp, R, nr);
		XS B = calculate_XS(mid, temp, R, nr);
		XS C = calculate_XS(high, temp, R, nr);

		double s_g = ((high-low)/6.0 * (A.sigma_g + 4.0*B.sigma_g + C.sigma_g));

		double D = s_g / (s_g + s_b );
		RI += D / mid;
	}
	RI *= s_b;

	double xs = RI / log(e2/e1);
	/*
	printf("T = %-4.0lfK  s_b = %-4.0lf  Range = %-2.0lf-%-2.0lf [eV] "
			"RIWR = %-8.3lf[b]  xs = %-8.3lf[b]\n",
			temp, s_b, e1, e2, RI, xs);
			*/
	printf("%-10.0lf%-10.0lf%7.0lf-%-7.0lf"
			"%-10s%-10.3lf%-10.3lf\n",
			temp, s_b, e1, e2,"Wide", RI, xs);
}

void find_NR_RI( double e1, double e2, int gp, double temp, double s_b )
{
	int nr;
	Resonance * R = res_read(&nr);
	nr = 3;
	double s_p = 11.2934;

	double range = e2 - e1;
	double del = range / gp;
	double RI = 0;

	for( int i = 0; i < gp; i++ )
	{
		double low = e1 + del*i;
		double high = e1 + del*(i+1);
		double mid = (low+high)/2;

		XS A = calculate_XS(low, temp, R, nr);
		XS B = calculate_XS(mid, temp, R, nr);
		XS C = calculate_XS(high, temp, R, nr);

		double s_g = ((high-low)/6.0 * (A.sigma_g + 4.0*B.sigma_g + C.sigma_g));
		double s_n = ((high-low)/6.0 * (A.sigma_n + 4.0*B.sigma_n + C.sigma_n));
		double s_t = ((high-low)/6.0 * (A.sigma_t + 4.0*B.sigma_t + C.sigma_t));

		double D = s_g / (s_t + s_b );
		RI += D / mid;
	}
	RI *= (s_p + s_b);

	double xs = RI / log(e2/e1);
	/*
	printf("T = %-4.0lfK  s_b = %-4.0lf  Range = %-2.0lf-%-2.0lf [eV] "
			"RINR = %-8.3lf[b]  xs = %-8.3lf[b]\n",
			temp, s_b, e1, e2, RI, xs);
			*/
	printf("%-10.0lf%-10.0lf%7.0lf-%-7.0lf"
			"%-10s%-10.3lf%-10.3lf\n",
			temp, s_b, e1, e2,"Narrow", RI, xs);
}

void find_RI( double e1, double e2, int gp, double temp )
{
	int nr;
	Resonance * R = res_read(&nr);
	nr = 3;

	double range = e2 - e1;
	double del = range / gp;
	double RI = 0;

	// Simpson's
	for( int i = 0; i < gp; i++ )
	{
		double low = e1 + del*i;
		double high = e1 + del*(i+1);
		double mid = (low+high)/2;

		double A = integral_RI(low, temp, R, nr);
		double B = integral_RI(mid, temp, R, nr);
		double C = integral_RI(high, temp, R, nr);

		RI += (high-low)/6.0 * (A + 4.0*B + C);
	}

	double xs = RI / log(e2/e1);
	/*
	printf("T = %-6.1lfK  Range = (%4.1lf - %-4.1lf) eV  "
			"RI = %-8.3lf[b]  xs = %-8.3lf[b]\n",
			temp, e1, e2, RI, xs);
			*/
	printf("%-10.0lf%-10s%7.0lf-%-7.0lf"
			"%-10s%-10.3lf%-10.3lf\n",
			temp, "-", e1, e2,"InfD", RI, xs);
}

double phi_RI( double E, double temp, Resonances * R, int nr )
{
	return 1.0/E;
}

double integral_RI( double E, double temp, Resonances * R, int nr )
{
	XS xs = calculate_XS(E, temp, R, nr);
	return xs.sigma_g * phi_RI(E, temp, R, nr);
}

double phi_RI_Narrow( double E, double temp, XS xs, Resonances * R, int nr, double s_b )
{
	double s_p = 11.2934;
	return 1.0 / E * ( (s_p + s_b) / (xs.sigma_t + s_b) );
}

double integral_RI_Narrow( double E, double temp, Resonances * R, int nr, double s_b )
{
	XS xs = calculate_XS(E, temp, R, nr);
	return xs.sigma_g * phi_RI_Narrow(E, temp, xs, R, nr);
}

double phi_RI_Wide( double E, double temp, XS xs, Resonances * R, int nr, double s_b )
{
	XS xs = calculate_XS(E, temp, R, nr);
	return 1.0 / E * ( s_b / (xs.sigma_g + s_b) );
}

double integral_RI_Wide( double E, double temp, Resonances * R, int nr, double s_b )
{
	XS xs = calculate_XS(E, temp, R, nr);
	return xs.sigma_g * phi_RI_Wide(E, temp, xs, R, nr);
}
