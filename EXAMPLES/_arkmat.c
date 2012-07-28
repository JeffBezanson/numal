int tel;
real_t h1k,h2k,hpi,h1,h2;

#include <stdio.h>
#include <math.h>

void der(int m, int n, real_t t, real_t **u, real_t **du)
{
	int i,j;

	n=n/2;
	for (i=2; i<=n-1; i++)
		for (j=2; j<=m-1; j++) {
			du[i][j]=u[i+n][j];
			du[i+n][j]=(u[i][j+1]-2.0*u[i][j]+u[i][j-1])/h1k+
							(u[i+1][j]-2.0*u[i][j]+u[i-1][j])/h2k;
		}
	for (j=1; j<=m; j += m-1) {
		inimat(n+1,n+n,j,j,du,0.0);
		for (i=1; i<=n; i++) du[i][j]=u[n+1][j];
	}
	for (i=1; i<=n; i += n-1)
		for (j=2; j<=m-1; j++) {
			du[i][j]=u[i+n][j];
			if (i == 1)
				du[n+1][j]=(u[1][j+1]-2.0*u[1][j]+u[1][j-1])/h1k+
								(2.0*u[2][j]-2.0*u[1][j])/h2k;
			else
				du[2*n][j]=0.0;
		}
}

void out(real_t t, real_t te, int m, int n, real_t **u, int type,
			int order, real_t *spr)
{
	int i;

	tel++;
	if (t == te) {
		for (i=1; i<=10; i++)
			printf("%6.3f  %6.3f  %9.6f   %9.6f\n",
				(i-1)*h1,(i-1)*h2,u[i][i],
			 sin(h1*(i-1))*cos(hpi*h2*(i-1))*cos(t*sqrt(1.0+hpi*hpi)));
		printf("\nThe number of integration steps: %2d\n"
			" Type is %1d   Order is %1d\n",tel,type,order);
	}
}

void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void arkmat(real_t *, real_t, int, int, real_t **,
		void (*)(int, int, real_t, real_t **, real_t **),
		int, int *, real_t *,
		void (*)(real_t, real_t, int, int, real_t **, int, int, real_t *));
	int i,j,n,m,typ,orde;
	real_t **u,t,te,cos1,spr;

	u=allocate_real_matrix(1,20,1,10);
	hpi=2.0*atan(1.0);
	h2=1.0/9.0;
	h1=(2.0*hpi)/9.0;
	n=m=10;
	h1k=h1*h1;
	h2k=h2*h2;
	tel=0;
	t=0.0;
	te=1.0;
	for (j=1; j<=m; j++) u[n][j]=sin(h1*(j-1));
	for (i=1; i<=n; i++) {
		cos1=cos(h2*hpi*(i-1));
		for (j=1; j<=m; j++) u[i][j]=u[n][j]*cos1;
	}
	inimat(n+1,n+n,1,m,u,0.0);
	typ=3;
	orde=2;
	spr=80.0;
	arkmat(&t,te,m,n+n,u,der,typ,&orde,&spr,out);
	free_real_matrix(u,1,20,1);
}

