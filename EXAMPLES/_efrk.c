#define real_t double

#include <math.h>
#include <stdio.h>

unsigned int passes;

void der(int m0, int m, real_t x, real_t y[])
{
	real_t y1,y2;

	y1=y[1];
	y2=y[2];
	y[1]=(y1+0.99)*(y2-1.0)+0.99;
	y[2]=1000.0*((1.0+y1)*(1.0-y2)-1.0);
	passes++;
}

void out(int m0, int m, real_t x, real_t xe, real_t y[],
			real_t *sigma, real_t *phi, real_t *diameter, int k,
			real_t *step, int r, int l)
{
	real_t s;

	s=(-1000.0*y[1]-1001.0+y[2])/2.0;
	*sigma=fabs(s-sqrt(s*s+10.0*(y[2]-1.0)));
	*diameter=2.0*(*step)*fabs(1000.0*
					(1.99*y[2]-2.0*y[1]*(1.0-y[2])));
	if (x == 50.0)
		printf(" %2d  %2d  %4d  %5u   %e   %e\n",
				r,l,k,passes,y[1],y[2]);
}

void main ()
{
	void efrk(real_t *, real_t, int, int, real_t [], real_t *,
			real_t *, real_t *, void (*)(int, int, real_t, real_t[]),
			int *, real_t *, real_t, real_t, real_t [], int, real_t,
			void (*)(int, int, real_t, real_t, real_t [], real_t *,
						real_t *, real_t *, int, real_t *, int, int));
	int k,r,l,i;
	real_t x,xe,sigma,phi,step,diameter,y[3],beta[7];

	printf("The results with EFRK are:\n\n"
		"  R   L    K   DER.EV.     Y[1]         Y[2]\n");
	phi=4.0*atan(1.0);
	beta[0]=beta[1]=1.0;
	for (r=1; r<=3; r++)
		for (l=1; l<=3; l++) {
			for (k=2; k<=r; k++) beta[k]=beta[k-1]/k;
			for (i=1; i<=2; i++) {
				step = ((i == 1) ? 1.0 : 0.1);
				passes=k=0;
				x=y[2]=0.0;
				y[1]=1.0;
				out(1,2,x,xe,y,&sigma,&phi,&diameter,k,&step,r,l);
				efrk(&x,50.0,1,2,y,&sigma,&phi,&diameter,der,&k,
						&step,r,l,beta,r>=3,1.0e-3,out);
			}
			printf("\n");
		}
}

