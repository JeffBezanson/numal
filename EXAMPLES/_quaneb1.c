#include <stdio.h>

int fun(int n, int l, int u, real_t x[], real_t f[])
{
	int i;
	real_t x1,x2,x3;

	x1 = (l == 1) ? 0.0 : x[l-1];
	x2=x[l];
	x3 = (l == n) ? 0.0 : x[l+1];
	for (i=l; i<=u; i++) {
		f[i]=(3.0-2.0*x2)*x2+1.0-x1-x3*2.0;
		x1=x2;
		x2=x3;
		x3 = (i <= n-2) ? x[i+2] : 0.0;
	}
	return (1);
}

void main ()
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void quanewbnd1(int, int, int, real_t [], real_t [],
						int (*)(int, int, int, real_t[], real_t[]),
						real_t [], real_t []);
	int i;
	real_t *x,*f,in[6],out[6];

	x=allocate_real_vector(1,600);
	f=allocate_real_vector(1,600);
	for (i=1; i<=600; i++) x[i] = -1.0;
	in[0]=1.0e-6;  in[1]=in[2]=in[3]=1.0e-5;  in[4]=20000.0;
	in[5]=0.001;
	quanewbnd1(600,1,1,x,f,fun,in,out);
	printf("Norm Residual vector: %e\n"
			"Length of last step:  %e\n"
			"Number of function component evaluations: %6.0f\n"
			"Number of iterations: %3.0f\nReport: %3.0f\n",
			out[2],out[1],out[3],out[4],out[5]);
	free_real_vector(x,1);
	free_real_vector(f,1);
}

