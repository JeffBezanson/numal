#include <math.h>
#include <stdio.h>

int nfe,nje,point;
real_t print[6];

void f(real_t t, real_t y[], real_t f1[], int n)
{
	nfe++;
	f1[1]=0.2*(y[2]-y[1]);
	f1[2]=10.0*y[1]-(60.0-0.125*y[3])*y[2]+0.125*y[3];
	f1[3]=1.0;
}

int available(real_t t, real_t y[], real_t **a, int n)
{
	nje++;
	a[1][1] = -0.2;
	a[1][2] = 0.2;
	a[1][3]=a[3][1]=a[3][2]=a[3][3]=0.0;
	a[2][1]=10.0;
	a[2][2]=0.125*y[3]-60.0;
	a[2][3]=0.125*(1.0+y[2]);
	return 1;
}

void update(real_t sw[], real_t r1[], int n)
{
	int i;
	real_t s1,s2;

	for (i=1; i<=n; i++) {
		s1=1.0/sw[i];
		s2=fabs(r1[i]);
		if (s1 < s2) sw[i]=1.0/s2;
	}
}

void control(real_t *tp, real_t t, real_t h, real_t hnew, real_t **y,
				real_t err[], int n, real_t tend)
{
	int i;
	real_t c[6],*x,s,s2,s3,s4;

	x=allocate_real_vector(1,n);
	while (1) {
		s=(t-(*tp))/h;
		s2=s*s;
		s3=s2*s;
		s4=s3*s;
		c[3]=(s2-s)/2.0;
		c[4] = -s3/6.0+s2/2.0-s/3.0;
		c[5]=s4/24.0-s3/4.0+11.0*s2/24.0-s/4.0;
		for (i=1; i<=n; i++)
			x[i]=y[1][i]-s*y[2][i]+c[3]*y[3][i]+
					c[4]*y[4][i]+c[5]*y[5][i];
		printf(" %6.2f  %7.2e  %e   %e   %4d  %3d\n",
			*tp,err[3],x[1],x[2],nfe,nje);
		if (*tp >= tend) break;
		point++;
		*tp = print[point];
		if (*tp > t) break;
	}
	free_real_vector(x,1);
}

void main ()
{
	void impex(int, real_t, real_t, real_t [],
				void (*)(real_t, real_t [], real_t [], int),
				int (*)(real_t, real_t [], real_t **, int),
				real_t, real_t, int, real_t, real_t [],
				void (*)(real_t [], real_t [], int), int *fail,
				void (*)(real_t *, real_t, real_t, real_t, real_t **,
							real_t [], int, real_t));
	void dupvec(int, int, int, real_t [], real_t []);
	real_t vecvec(int, int, int, real_t [], real_t []);
	void elmvec(int, int, int, real_t [], real_t [], real_t);
	int n,fail,i,it;
	real_t t,tend,eps,hmax,l,h2,y[4],sw[4],f1[4],f2[4],z[4],x[4],
			n1,n2;

	printf("The results with IMPEX are:\n\n");
	n=3;
	nje=nfe=0;
	t=0.0;
	tend=400.0;
	eps=1.0e-5;
	hmax=400.0;
	y[1]=y[2]=y[3]=0.0;
	sw[1]=sw[2]=sw[3]=1.0;
	print[1]=0.1;  print[2]=1.0;  print[3]=10.0;  print[4]=100.0;
	print[5]=400.0;
	dupvec(1,n,0,z,y);
	for (i=1; i<=n; i++)
		x[i] = (y[i] == 0.0) ? eps : (1.0+eps)*y[i];
	n1=sqrt(vecvec(1,n,0,x,x))*eps;
	f(t,x,f1,n);
	for (it=1; it<=5; it++) {
		f(t,z,f2,n);
		elmvec(1,n,0,f2,f1,-1.0);
		n2=n1/sqrt(vecvec(1,n,0,f2,f2));
		dupvec(1,n,0,z,x);
		elmvec(1,n,0,z,f2,n2);
	}
	f(t,z,f2,n);
	elmvec(1,n,0,f2,f1,-1.0);
	l=sqrt(vecvec(1,n,0,f2,f2))/n1;
	h2=pow(eps*320.0,1.0/5.0)/(4.0*l);
	printf("EPS = %e\nInterval of integration = (%1.0f,%3.0f)\n"
		"Maximally allowed stepsize = %e\n\nLipschconst = %e\n"
		"Starting stepsize = %e\nFunctional eval = %2d\n\n"
		"    X    ERROR       Y[1]          Y[2]        NFE  NJE\n",
		eps,t,tend,hmax,l,h2,nfe);
	impex(n,t,tend,y,f,available,h2,hmax,0,eps,sw,update,&fail,
			control);
	printf("\nNumber of functional evaluations =%4d\n"
		"Number of Jacobian evaluations   = %3d\n",nfe,nje);
}

