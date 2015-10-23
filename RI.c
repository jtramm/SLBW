#include "slbw_header.h"

void find_specific_NR_RI( double e1, double e2, int gp, double temp, double s_b )
{
	int nr;
	Resonance * R = (Resonance *) malloc(sizeof(Resonance));
	R->Eo = 1050.0;
	R->Tn = 0.023;
	R->Tg = 0.095;
	nr = 1;
	double s_p = 11.4;

	double range = e2 - e1;
	double del = range / gp;
	double RI = 0;

	for( int i = 0; i < gp; i++ )
	{
		double low = e1 + del*i;
		double high = e1 + del*(i+1);
		double mid = (low+high)/2;

		double A = integral_RI_Narrow(low, temp, R, nr, s_b);
		double B = integral_RI_Narrow(mid, temp, R, nr, s_b);
		double C = integral_RI_Narrow(high, temp, R, nr, s_b);

		RI += (high-low)/6.0 * (A + 4.0*B + C);
	}
	double phi = 0;
	for( int i = 0; i < gp; i++ )
	{
		double low = e1 + del*i;
		double high = e1 + del*(i+1);
		double mid = (low+high)/2;

		double A = phi_RI_Narrow(low, temp, R, nr, s_b);
		double B = phi_RI_Narrow(mid, temp, R, nr, s_b);
		double C = phi_RI_Narrow(high, temp, R, nr, s_b);

		phi += (high-low)/6.0 * (A + 4.0*B + C);
	}

	double xs = RI / phi;
	/*
	printf("T = %-4.0lfK  s_b = %-4.0lf  Range = %-2.0lf-%-2.0lf [eV] "
			"RINR = %-8.3lf[b]  xs = %-8.3lf[b]\n",
			temp, s_b, e1, e2, RI, xs);
			*/
	printf("%-10.0lf%-10.0lf%7.0lf-%-7.0lf"
			"%-10s%-10.3lf%-10.6lf\n",
			temp, s_b, e1, e2,"Narrow", RI, xs);
}


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

		double A = integral_RI_Wide(low, temp, R, nr, s_b);
		double B = integral_RI_Wide(mid, temp, R, nr, s_b);
		double C = integral_RI_Wide(high, temp, R, nr, s_b);

		RI += (high-low)/6.0 * (A + 4.0*B + C);
	}
	
	double phi = 0;
	for( int i = 0; i < gp; i++ )
	{
		double low = e1 + del*i;
		double high = e1 + del*(i+1);
		double mid = (low+high)/2;

		double A = phi_RI_Wide(low, temp, R, nr, s_b);
		double B = phi_RI_Wide(mid, temp, R, nr, s_b);
		double C = phi_RI_Wide(high, temp, R, nr, s_b);

		phi += (high-low)/6.0 * (A + 4.0*B + C);
	}

	double xs = RI / phi;
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

		double A = integral_RI_Narrow(low, temp, R, nr, s_b);
		double B = integral_RI_Narrow(mid, temp, R, nr, s_b);
		double C = integral_RI_Narrow(high, temp, R, nr, s_b);

		RI += (high-low)/6.0 * (A + 4.0*B + C);
	}
	double phi = 0;
	for( int i = 0; i < gp; i++ )
	{
		double low = e1 + del*i;
		double high = e1 + del*(i+1);
		double mid = (low+high)/2;

		double A = phi_RI_Narrow(low, temp, R, nr, s_b);
		double B = phi_RI_Narrow(mid, temp, R, nr, s_b);
		double C = phi_RI_Narrow(high, temp, R, nr, s_b);

		phi += (high-low)/6.0 * (A + 4.0*B + C);
	}

	double xs = RI / phi;
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

double phi_RI( double E, double temp, Resonance * R, int nr )
{
	return 1.0/E;
}

double integral_RI( double E, double temp, Resonance * R, int nr )
{
	XS xs = calculate_XS(E, temp, R, nr);
	return xs.sigma_g * phi_RI(E, temp, R, nr);
}

double phi_RI_Narrow( double E, double temp, Resonance * R, int nr, double s_b )
{
	double s_p = 11.2934;
	XS xs = calculate_XS(E, temp, R, nr);
	return 1.0 / E * ( (s_p + s_b) / (xs.sigma_t + s_b) );
}

double integral_RI_Narrow( double E, double temp, Resonance * R, int nr, double s_b )
{
	XS xs = calculate_XS(E, temp, R, nr);
	return xs.sigma_g * phi_RI_Narrow(E, temp, R, nr, s_b);
}

double phi_RI_Wide( double E, double temp, Resonance * R, int nr, double s_b )
{
	XS xs = calculate_XS(E, temp, R, nr);
	return 1.0 / E * ( s_b / (xs.sigma_g + s_b) );
}

double integral_RI_Wide( double E, double temp, Resonance * R, int nr, double s_b )
{
	XS xs = calculate_XS(E, temp, R, nr);
	return xs.sigma_g * phi_RI_Wide(E, temp, R, nr, s_b);
}

