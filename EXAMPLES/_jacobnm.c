#include <stdio.h>

void f1(int n, int m, real_t x[], real_t f[])
{
	f[1]=x[1]*x[1]*x[1]+x[2];
	f[2]=x[2]*10.0+x[2]*x[1]*x[1];
	f[3]=x[1]*x[2];
}

real_t di(int i)
{
	return ((i == 2) ? 1.0 : 1.0e-5);
}

void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void jacobnmf(int, int, real_t [], real_t [], real_t **,
					real_t (*)(int), void (*)(int, int, real_t[], real_t[]));
	int i;
	real_t **jac,x[3],f[4];

	jac=allocate_real_matrix(1,3,1,2);
	x[1]=2.0;	x[2]=1.0;
	f1(3,2,x,f);
	jacobnmf(3,2,x,f,jac,di,f1);
	printf("The calculated jacobian is:\n"
		"   %7.1f  %7.1f\n   %7.1f  %7.1f\n   %7.1f  %7.1f\n",
		jac[1][1],jac[1][2],jac[2][1],jac[2][2],jac[3][1],jac[3][2]);
	free_real_matrix(jac,1,3,1);
}

