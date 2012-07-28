#include <math.h>
#include <stdio.h>

real_t fxyj(int n, int j, real_t x, real_t y[])
{
	return ((j == 1) ? y[2] : -y[1]);
}

void main ()
{
	void rk3n(real_t *, real_t, real_t, real_t [], real_t [], real_t [],
				real_t [], real_t (*)(int, int, real_t, real_t[]),
				real_t [], real_t [], int, int);
	int k,fi,i,n,j;
	real_t b,x,y[3],ya[3],z[3],e[9],d[8],x2,term;

	printf("Results from RK3N :\n");
	for (k=1; k<=8; k++) e[k]=1.0e-5;
	fi=1;
	y[1]=y[2]=1.0;
	z[1]=z[2]=0.0;
	b=0.0;
	do {
		b += 1.0;
		rk3n(&x,0.0,b,y,y,z,z,fxyj,e,d,fi,2);
		ya[1]=ya[2]=0.0;
		term=1.0;
		x2=x*x*0.5;
		n=1;
		do {
			for (i=1; i<=2; i++) {
				j=(i+n-2)/2;
				ya[i] += term*((j%2 == 0) ? 1 : -1);
			}
			term=term*x2/n/(n*2-1);
			n++;
		} while (fabs(term) > 1.0e-14);
		printf(" ABS(YEXACT[1]-Y[1])+ABS(YEXACT[2]-Y[2]) = %e\n",
				fabs(y[1]-ya[1])+fabs(ya[2]-y[2]));
		fi=0;
	} while (b < 5.0);
}

