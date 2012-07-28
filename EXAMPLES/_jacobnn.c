#include <stdio.h>

void f1(int n, real_t x[], real_t f[])
{
	f[1]=x[1]*x[1]*x[1]+x[2];
	f[2]=x[2]*10.0;
}

real_t di(int i)
{
	return ((i == 1) ? 1.0e-5 : 1.0);
}

void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void jacobnnf(int, real_t [], real_t [], real_t **,
					real_t (*)(int), void (*)(int, real_t[], real_t[]));
	int i;
	real_t **jac,x[3],f[3];

	jac=allocate_real_matrix(1,2,1,2);
	x[1]=2.0;	x[2]=1.0;
	f1(2,x,f);
	jacobnnf(2,x,f,jac,di,f1);
	printf("The calculated jacobian is:\n %6.1f %6.1f\n %6.1f %6.1f\n",
			jac[1][1],jac[1][2],jac[2][1],jac[2][2]);
	free_real_matrix(jac,1,2,1);
}

