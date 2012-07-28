#include <math.h>
#include <stdio.h>

int passes,pasjac;

real_t f(int m, real_t y[], int i, real_t *sigma1, real_t *sigma2)
{
	if (i == 1)
		return (y[1]+0.99)*(y[2]-1.0)+0.99;
	else {
		passes++;
		return 1000.0*((1.0+y[1])*(1.0-y[2])-1.0);
	}
}

void jacobian(int m, real_t **j, real_t y[],
					real_t *sigma1, real_t *sigma2)
{
	j[1][1]=y[2]-1.0;
	j[1][2]=0.99+y[1];
	j[2][1]=1000.0*(1.0-y[2]);
	j[2][2] = -1000.0*(1.0+y[1]);
	*sigma1=fabs(j[2][2]+j[1][1]-sqrt((j[2][2]-j[1][1])*
					(j[2][2]-j[1][1])+4.0*j[2][1]*j[1][2]))/2.0;
	pasjac++;
}

int evaluate1(int i)
{
	return (i == 1);
}

int evaluate2(int i)
{
	return 1;
}

void out(real_t x, real_t xe, int m, real_t y[],
			real_t sigma1, real_t sigma2, real_t **j, int k)
{
	if (x == 50.0)
		printf("%3d  %4d    %4d     %e   %e\n",
				k,passes,pasjac,y[1],y[2]);
}

void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void liniger2(real_t *, real_t, int, real_t [], real_t *, real_t *,
			real_t (*)(int, real_t[], int, real_t *, real_t *),
			int (*)(int), real_t **,
			void (*)(int, real_t **, real_t [], real_t *, real_t *),
			int *, int, real_t, real_t, real_t,
			void (*)(real_t, real_t, int, real_t [], real_t, real_t,
						real_t **, int));
	int i,k,itmax;
	real_t x,sigma1,sigma2,step,y[3],**j;

	j=allocate_real_matrix(1,2,1,2);
	printf("The results with LINIGER2 (second order) are:\n"
		" K   DER.EV. JAC.EV.     Y[1]         Y[2]\n");
	for (i=1; i<=2; i++) {
		step = (i == 1) ? 10.0 : 1.0;
		for (itmax=1; itmax<=3; itmax += 2) {
			passes=pasjac=0;
			x=y[2]=0.0;
			y[1]=1.0;
			sigma2=0.0;
			liniger2(&x,50.0,2,y,&sigma1,&sigma2,f,evaluate1,j,
						jacobian,&k,itmax,step,1.0e-4,1.0e-4,out);
		}
	}
	printf("\nThe results with LINIGER2 (third order) are:\n"
		" K   DER.EV. JAC.EV.     Y[1]         Y[2]\n");
	for (i=1; i<=2; i++) {
		step = (i == 1) ? 10.0 : 1.0;
		for (itmax=1; itmax<=3; itmax += 2) {
			passes=pasjac=0;
			x=y[2]=0.0;
			y[1]=1.0;
			sigma2=0.0;
			liniger2(&x,50.0,2,y,&sigma1,&sigma2,f,evaluate2,j,
						jacobian,&k,itmax,step,1.0e-4,1.0e-4,out);
		}
	}
	free_real_matrix(j,1,2,1);
}

