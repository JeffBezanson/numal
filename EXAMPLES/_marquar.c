#include <math.h>
#include <stdio.h>

real_t x[7],y[7];

int expfunct(int m, int n,real_t par[], real_t rv[])
{
	int i;

	for (i=1; i<=m; i++) {
		if (par[3]*x[i] > 680.0) return 0;
		rv[i]=par[1]+par[2]*exp(par[3]*x[i])-y[i];
	}
	return 1;
}

void jacobian(int m, int n, real_t par[], real_t rv[], real_t **jac)
{
	int i;
	real_t ex;

	for (i=1; i<=m; i++) {
		jac[i][1]=1.0;
		jac[i][2]=ex=exp(par[3]*x[i]);
		jac[i][3]=x[i]*par[2]*ex;
	}
}

void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void marquardt(int, int, real_t [], real_t [], real_t **,
						int (*)(int, int, real_t[], real_t[]),
						void (*)(int, int, real_t[], real_t[], real_t **),
						real_t [], real_t []);
	real_t in[7],out[8],rv[7],par[4],**jjinv;

	jjinv=allocate_real_matrix(1,3,1,3);
	in[0]=1.0e-6;  in[3]=1.0e-4;  in[4]=1.0e-1;  in[5]=75.0;
	in[6]=1.0e-2;
	x[1] = -5.0;  x[2] = -3.0;  x[3] = -1.0;  x[4]=1.0;
	x[5]=3.0;  x[6]=5.0;
	y[1]=127.0;  y[2]=151.0;  y[3]=379.0;  y[4]=421.0;
	y[5]=460.0;  y[6]=426.0;
	par[1]=580.0;  par[2] = -180.0;  par[3] = -0.160;
	marquardt(6,3,par,rv,jjinv,expfunct,jacobian,in,out);
	printf("Parameters:\n   %9.4e   %9.4e   %9.4e\n\nOUT:\n"
		" %14.6e\n %14.6e\n %14.6e\n %14.6e\n %14.6e\n %14.6e\n"
		" %14.6e\n\nLast residual vector:\n"
		" %6.1f  %6.1f  %6.1f  %6.1f  %6.1f  %6.1f\n",
		par[1],par[2],par[3],out[7],out[2],out[6],out[3],out[4],
		out[5],out[1],rv[1],rv[2],rv[3],rv[4],rv[5],rv[6]);
	free_real_matrix(jjinv,1,3,1);
}

