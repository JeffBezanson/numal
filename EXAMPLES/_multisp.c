#include <stdio.h>

void der(real_t f[], int n, real_t x, real_t y[])
{
	real_t r;

	f[2]=r=3.0e7*y[1]*y[1];
	f[1]=0.04*(1-y[1]-y[2])-1.0e4*y[1]*y[2]-r;
}

int avail(int n, real_t x, real_t y[], real_t **jac)
{
	real_t r;

	jac[2][1]=r=6.0e7*y[1];
	jac[1][1] = -0.04-1.0e4*y[2]-r;
	jac[1][2] = -0.04-1.0e4*y[1];
	jac[2][2]=0.0;
	return 1;
}

void out(real_t h, int k, int n, real_t x, real_t y[])
{
	return;
}

void main ()
{
	real_t *allocate_real_vector(int, int);
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_vector(real_t *, int);
	void free_real_matrix(real_t **, int, int, int);
	int multistep(real_t *, real_t, real_t [], real_t, real_t,
				real_t [], real_t, int *, real_t [],
				void (*)(real_t [], int, real_t, real_t[]),
				int (*)(int, real_t, real_t [], real_t **), real_t **, int,
				int, void (*)(real_t, int, int, real_t, real_t []));
	int i,first;
	real_t x,xend,hmin,eps,y[13],ymax[3],*d,**jac;

	d=allocate_real_vector(-40,12);
	jac=allocate_real_matrix(1,2,1,2);
	hmin=1.0e-6;
	eps=1.0e-6;
	first=1;
	x=0.0;
	y[1]=y[2]=0.0;
	ymax[1]=0.0001;
	ymax[2]=1.0;
	printf("Delivers with hmin = %5.1e  and eps = %5.1e\n",hmin,eps);
	for (i=1; i<=10; i+=9) {
		xend=i;
		multistep(&x,xend,y,hmin,5,ymax,eps,&first,d,
					der,avail,jac,1,2,out);
		printf("   %e   %e\n",y[1],y[2]);
	}
	free_real_vector(d,-40);
	free_real_matrix(jac,1,2,1);
}

