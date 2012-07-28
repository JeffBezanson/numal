#include <stdio.h>

real_t b(int n, real_t x[])
{
	return x[2];
}

real_t fxj(int n, int k, real_t x[])
{
	return ((k == 1) ? x[2] : 10.0*(1.0-x[1]*x[1])*x[2]-x[1]);
}

void main ()
{
	void rk4na(real_t [], real_t [], real_t (*)(int, real_t[]),
					real_t (*)(int, int, real_t[]), real_t [], real_t [],
					int, int, int, int);
	int j,first;
	real_t x0,e[8],xa[3],x[3],d[6];

	for (j=0; j<=5; j++) e[j]=0.1e-4;
	e[6]=e[7]=1.0e-4;
	printf("VAN DER POL\n\nEPS = %e\n\n"
			"The values of x[0],x[1],x[2],p:\n",e[0]);
	x0=xa[0]=xa[2]=0.0;
	xa[1]=2.0;
	printf("  %8.5f  %8.5f  %8.5f  %8.5f\n",xa[0],xa[1],xa[2],x0);
	first=1;
	for (j=1; j<=4; j++) {
		rk4na(x,xa,b,fxj,e,d,first,2,0,1);
		x0=x[0]-x0;
		printf("  %8.5f  %8.5f  %8.5f  %8.5f\n",x[0],x[1],x[2],x0);
		first=0;
		x0=x[0];
	}
}

