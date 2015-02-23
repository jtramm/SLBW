#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>

#define N 10
#define Tm 12.0

int main(void)
{
	double complex prefactor = I * Tm * Tm / sqrt(M_PI);
	double an[N+1] = {0};
	double neg_1n[N+1] = {0};
	double denominator_left[N+1] = {0};
	for( int n = 1; n <= N; n++ )
	{
		an[n] = 2.0 * sqrt(M_PI) / Tm * exp( -n*n*M_PI*M_PI / (Tm*Tm));
		neg_1n[n] = pow(-1.0,n);	
		denominator_left[n] = n*n*M_PI*M_PI;
	}
	
	printf("prefactor = %e + %e i\n", creal(prefactor), cimag(prefactor));
	printf("an\n");
	for( int n = 1; n <= N; n++ )
		printf("%e\n", an[n]);
	printf("neg_1n\n");
	for( int n = 1; n <= N; n++ )
		printf("%e\n", neg_1n[n]);
	printf("denominator_left\n");
	for( int n = 1; n <= N; n++ )
		printf("%e\n", denominator_left[n]);
	return 0;
}
