#include <stdio.h>

real_t f(real_t x)
{
	return 1.0/(x-10.0);
}

void compute(int n, real_t a, real_t b, real_t (*f)(real_t))
{
	int k,l,m;
	real_t r,idm,*coef,em[4],*y,*fy;
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void minmaxpol(int, int, real_t [], real_t [], real_t [], real_t []);

	coef=allocate_real_vector(0,n);
	em[2]=10*n+5;
	m=100*n+10;
	y=allocate_real_vector(0,m);
	fy=allocate_real_vector(0,m);
	idm=(b-a)/m;
	r=y[0]=a;
	fy[0]=f(r);
	r=y[m]=b;
	fy[m]=f(r);
	l=m-1;
	for (k=1; k<=l; k++) {
		r=y[k]=a+k*idm;
		fy[k]=f(r);
	}
	minmaxpol(n,m,y,fy,coef,em);
	printf(" COEF:   ");
	for (k=0; k<=n; k++) printf("  %e",coef[k]);
	printf("\n EM[0:3]:  %e   %e   %4.0f   %3.0f\n",
			em[0],em[1],em[2],em[3]);
	free_real_vector(coef,0);
	free_real_vector(y,0);
	free_real_vector(fy,0);
}

void main ()
{
	int n;

	printf("MINMAXPOL delivers:\n");
	n=1;
	printf(" Degree = %d\n",n);
	compute(n,-1.0,1.0,f);
}

