#include <stdio.h>

real_t fxyzj(int n, int j, real_t x, real_t y[], real_t z[])
{
	return -5.0*(y[j]+z[j])+((j == 1) ? y[2] : y[1]);
}

void main ()
{
	void rk2n(real_t *, real_t, real_t, real_t [], real_t [], real_t [],
			real_t [], real_t (*)(int, int, real_t, real_t [], real_t []),
			real_t [], real_t [], int, int);
	int k,fi;
	real_t b,x,expx,y[3],ya[3],z[3],za[3],e[9],d[8];

	printf("Results from RK2N :\n");
	for (k=1; k<=8; k++) e[k]=1.0e-5;
	ya[1]=za[2]=1.0;
	ya[2]=za[1]=0.0;
	b=1.0;
	do {
		fi=(b == 1.0);
		rk2n(&x,0.0,b,y,ya,z,za,fxyzj,e,d,fi,2);
		expx=exp(-x);
		ya[1] = -expx*(expx*(expx*(expx/3.0+0.5)-1.0)-5.0/6.0);
		ya[2] = -expx*(expx*(expx*(expx/3.0-0.5)+1.0)-5.0/6.0);
		za[1]=expx*(expx*(expx*(expx/0.75+1.5)-2.0)-5.0/6.0);
		za[2]=expx*(expx*(expx*(expx/0.75-1.5)+2.0)-5.0/6.0);
		printf("\n  X = %6.4f\nY[1]-YEXACT[1] = %+6.2e    "
			"Y[2]-YEXACT[2] = %+6.2e\nZ[1]-ZEXACT[1] = %+6.2e    "
			"Z[2]-ZEXACT[2] = %+6.2e\n",x,y[1]-ya[1],y[2]-ya[2],
			z[1]-za[1],z[2]-za[2]);
		b += 1.0;
	} while (b < 5.0);
}

