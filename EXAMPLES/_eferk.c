#include <math.h>
#include <stdio.h>

int passes,pasjac;

void der(int m, real_t y[])
{
	real_t y1,y2;

	y1=y[1];
	y2=y[2];
	y[1]=(y1+0.99)*(y2-1.0)+0.99;
	y[2]=1000.0*((1.0+y1)*(1.0-y2)-1.0);
	passes++;
}

void jacobian(int m, real_t **j, real_t y[], real_t *sigma)
{
	j[1][1]=y[2]-1.0;
	j[1][2]=0.99+y[1];
	j[2][1]=1000.0*(1.0-y[2]);
	j[2][2] = -1000.0*(1.0+y[1]);
	*sigma=fabs(j[2][2]+j[1][1]-sqrt((j[2][2]-j[1][1])*
					(j[2][2]-j[1][1])+4.0*j[2][1]*j[1][2]))/2.0;
	pasjac++;
}

void out(real_t x, real_t xe, int m, real_t y[], real_t **j, int k)
{
	if (x == 50.0)
		printf("%3d  %4d    %4d     %e   %e\n",
				k,passes,pasjac,y[1],y[2]);
}

void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void eferk(real_t *, real_t, int, real_t [], real_t *, real_t,
			void (*)(int, real_t[]), real_t **,
			void (*)(int, real_t **, real_t [], real_t *),
			int *, int, int, real_t, real_t, real_t, real_t, int,
			void (*)(real_t, real_t, int, real_t [], real_t **, int));
	int i,k;
	real_t x,xe,sigma,phi,tol,y[3],**j;

	j=allocate_real_matrix(1,2,1,2);
	printf("The results with EFERK are:\n\n"
		" K   DER.EV. JAC.EV.     Y[1]         Y[2]\n");
	phi=4.0*atan(1.0);
	tol=1.0;
	for (i=1; i<=4; i++) {
		passes=pasjac=0;
		x=y[2]=0.0;
		y[1]=1.0;
		xe=50.0;
		eferk(&x,xe,2,y,&sigma,phi,der,j,jacobian,&k,1,1,tol,tol,
				1.0e-6,50.0,0,out);
		tol *= 1.0e-1;
	}
	free_real_matrix(j,1,2,1);
}

