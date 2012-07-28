#include <math.h>
#include <stdio.h>

int passes,k;

void der(int n, real_t x, real_t y[], real_t dy[])
{
	real_t mu,mu1,y1,y2,y3,y4,s1,s2;

	mu=1.0/82.45;
	mu1=1.0-mu;
	passes++;
	y1=y[1];
	y2=dy[1]=y[2];
	y3=y[3];
	y4=dy[3]=y[4];
	s1=(y1+mu)*(y1+mu)+y3*y3;
	s2=(y1-mu1)*(y1-mu1)+y3*y3;
	s1 *= sqrt(s1);
	s2 *= sqrt(s2);
	dy[2]=y1+2.0*y4-mu1*(y1+mu)/s1-mu*(y1-mu1)/s2;
	dy[4]=y3-2.0*y2-mu1*y3/s1-mu*y3/s2;
}

void out(int n, real_t x, real_t xe, real_t y[], real_t s[])
{
	k++;
	if (x >= xe)
		printf(" %3d   %4d    %e   %e\n",k,passes,y[1],y[3]);
}

void main ()
{
	int i;
	real_t x,xe,tol,h0,y[5],s[5];

	printf("Results with DIFFSYS are :\n"
			"  K   DER.EV.      Y[1]          Y[3]\n");
	tol=1.0e-2;
	for (i=1; i<=2; i++) {
		tol *= 1.0e-2;
		passes=k=0;
		x=0.0;
		xe=6.192169331396;
		y[1]=1.2;
		y[2]=y[3]=0.0;
		y[4] = -1.04935750983;
		s[1]=s[2]=s[3]=s[4]=0.0;
		h0=0.2;
		diffsys(&x,xe,4,y,der,tol,tol,s,h0,out);
	}
}

