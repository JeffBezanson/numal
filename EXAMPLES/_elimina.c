#include <stdio.h>

int nn;
real_t h,h2,d1,d2;

void residual(int lj, int uj, int ll, int ul, real_t **u)
{
	int ujmin1,ulmin1,ljplus1,j,l;
	real_t u2,u1[12];

	ujmin1=uj-1;
	ulmin1=ul-1;
	ljplus1=lj+1;
	for (j=lj; j<=uj; j++) {
		u1[j]=u[j][ll];
		u[j][ll]=0.0;
	}
	for (l=ll+1; l<=ulmin1; l++) {
		u1[lj]=u[lj][l];
		u[lj][l]=0.0;
		for (j=ljplus1; j<=ujmin1; j++) {
			u2=u[j][l];
			u[j][l]=(4.0*u2-u1[j-1]-u1[j]-u[j+1][l]-u[j][l+1])+
						2.0*((j*h)*(j*h)+(l*h)*(l*h))*h2;
			u1[j]=u2;
		}
		u[uj][l]=0.0;
	}
	for (j=lj; j<=uj; j++) u[j][ul]=0.0;
}

void out2(real_t **u, int lj, int uj, int ll, int ul, int *p,
			real_t discr[], int k, real_t rateconv, real_t domeigval)
{
	if (k == *p)
		printf(" %2d  %e   %e   %e\n",k,discr[1],discr[2],rateconv);
}

void out1(real_t **u, int lj, int uj, int ll, int ul, int *n,
			real_t discr[], int k, real_t rateconv, real_t domeigval)
{
	if (k == 0)
		d1=d2=1.0;
	else {
		d2=d1;
		d1=domeigval;
		*n = (fabs((d1-d2)/d2) < 1.0e-4) ? k : nn;
		if (k == *n)
			printf("  K    DISCR[1]      DISCR[2]      RATECONV\n"
				" %2d  %e   %e   %e\n",k,discr[1],discr[2],rateconv);
	}
}

void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void richardson(real_t **, int, int, int, int,
			int, void (*)(int, int, int, int, real_t **),
			real_t, real_t, int *, real_t [], int *, real_t *, real_t *,
			void (*)(real_t **, int, int, int, int, int *, real_t [],
						int, real_t, real_t));
	void elimination(real_t **, int, int, int, int,
			void (*)(int, int, int, int, real_t **),
			real_t, real_t, int *, real_t [], int *, real_t *, real_t *,
			void (*)(real_t **, int, int, int, int, int *, real_t [],
							int, real_t, real_t));
	int j,l,lj,uj,ll,ul,n,p,k;
	real_t pi,domeigval,rateconvr,rateconve,rateconv,a,b,discr[3],**u;

	u=allocate_real_matrix(0,11,0,11);
	printf("RICHARDSON and ELIMINATION deliver:\n\n");
	pi=3.14159265358979;
	lj=0;  uj=11;  ll=0;  ul=11;  n=50;
	a=0.326;  b=7.83;
	h=pi/(uj-lj);
	h2=h*h;
	for (j=lj; j<=uj; j++)
		for (l=ll; l<=ul; l++)
			u[j][l] = (j==lj || j==uj || l==ll || l==ul) ?
							(j*h)*(j*h)*(l*h)*(l*h) : 1.0;
	nn=n;
	richardson(u,lj,uj,ll,ul,1,residual,a,b,&n,discr,&k,&rateconv,
					&domeigval,out1);
	rateconvr=rateconv;
	printf("\n dominant eigenvalue:  %e\n\n",domeigval);
	elimination(u,lj,uj,ll,ul,residual,a,b,&p,discr,&k,&rateconv,
					&domeigval,out2);
	rateconve=rateconv;
	nn=n+p;
	printf("\nTotal number of iterations: %2d\n"
		"Rate of convergence with respect to\n"
		"   the zeroth iterand of RICHARDSON:   %e\n",
		nn,(n*rateconvr+p*rateconve)/nn);
	free_real_matrix(u,0,11,0);
}

