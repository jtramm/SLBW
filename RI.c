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

	for( int i = 0; i < gp; i++ )
	{
		double low = e1 + del*i;
		double high = e1 + del*(i+1);
		double mid = (low+high)/2;

		XS A = calculate_XS(low, temp, R, nr);
		XS B = calculate_XS(mid, temp, R, nr);
		XS C = calculate_XS(high, temp, R, nr);

		RI += ((high-low)/6.0 * (A.sigma_g + 4.0*B.sigma_g + C.sigma_g)) / mid;
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
