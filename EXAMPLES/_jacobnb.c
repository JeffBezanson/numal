#include <stdio.h>

int func(int n, int l, int u, real_t x[], real_t f[])
{
	int i;

	for (i=l; i<=((u == 5) ? 4 : u); i++) {
		f[i]=(3.0-2.0*x[i])*x[i]+1.0-2.0*x[i+1];
		if (i != 1) f[i] -= x[i-1];
	}
	if (u == 5) f[5]=4.0-x[4]-x[5]*2.0;
}

real_t di(int i)
{
	return ((i == 5) ? 1.0 : 1.0e-6);
}

void main ()
{
	void jacobnbndf(int, int, int, real_t [], real_t [],
						real_t [], real_t (*)(int),
						int (*)(int, int, int, real_t[], real_t[]));
	int i;
	real_t jac[14],x[6],f[6];

	for (i=1; i<=5; i++) x[i] = -1.0;
	func(5,1,5,x,f);
	jacobnbndf(5,1,1,x,f,jac,di,func);
	printf("The calculated tridiagonal jacobian is:\n"
		"   %4.0f  %4.0f\n   %4.0f  %4.0f  %4.0f\n"
		"         %4.0f  %4.0f  %4.0f\n"
		"               %4.0f  %4.0f  %4.0f\n"
		"                     %4.0f  %4.0f\n",
		jac[1],jac[2],jac[3],jac[4],jac[5],jac[6],jac[7],jac[8],
		jac[9],jac[10],jac[11],jac[12],jac[13]);
}

