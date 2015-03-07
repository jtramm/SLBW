#include "slbw_header.h"

int main(void)
{

	RI_driver();
	NR_WR_Driver();

	return 0;
}

void NR_WR_Driver(void)
{
	int gp = 10000;
	double low[3] = {6, 10, 25};
	double high[3] = {10, 25, 50};
	double s_b[3] = {2000, 200, 20};
	double temp = 300;
	printf("%-10s%-10s%-15s%-10s%-10s%-10s\n",
			"T [K]",
			"s_b [b]",
			"Range [eV]",
			"Type",
			"RI [b]",
			"XS [b]"
			);
	for( int i = 0; i < 3; i++ )
		for( int j = 0; j < 3; j++ )
		{
			find_NR_RI( low[j], high[j], gp, temp, s_b[i] );
			find_WR_RI( low[j], high[j], gp, temp, s_b[i] );
		}
}

void RI_driver(void)
{
	int gp = 1000;
	double low, high, temp, RI, xs, s_b;
	temp = 300.0;
	low = 6.0; high = 10.0;
	printf("%-10s%-10s%-15s%-10s%-10s%-10s\n",
			"T [K]",
			"s_b [b]",
			"Range [eV]",
			"Type",
			"RI [b]",
			"XS [b]"
			);
	find_RI( low, high, gp, temp );
	low = 10.0; high = 25.0;
	find_RI( low, high, gp, temp );
	low = 25.0; high = 50.0;
	find_RI( low, high, gp, temp );
	temp = 1000.0;
	low = 6.0; high = 10.0;
	find_RI( low, high, gp, temp );
	low = 10.0; high = 25.0;
	find_RI( low, high, gp, temp );
	low = 25.0; high = 50.0;
	find_RI( low, high, gp, temp );
}


void graph_driver(void)
{
	int gp = 10000;
	double temp = 0.01;
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
		xs[i] = calculate_XS(E[i], temp, R, nr);

	res_out(xs, gp);

	free(E);
	free(xs);
	free(R);
}
