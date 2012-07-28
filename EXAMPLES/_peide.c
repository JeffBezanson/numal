#include <math.h>
#include <stdio.h>

void communication(int post, real_t fa, int n, int m, int nobs,
		int nbp, real_t par[], real_t res[], int bp[], real_t **jtjinv,
		real_t in[], real_t out[], int weight, int nis)
{
	real_t vecvec(int, int, int, real_t [], real_t []);
	int i,j;
	real_t c,*conf;

	conf=allocate_real_vector(1,m);
	if (post == 5) {
		printf("\nThe first residual vector\n  I      RES[I]\n");
		for (i=1; i<=nobs; i++) printf(" %2d   %+e\n",i,res[i]);
	} else if (post == 3) {
		printf("\nThe Euclidean norm of the residual vector : %e\n"
				"Calculated parameters:\n",
				sqrt(vecvec(1,nobs,0,res,res)));
		for (i=1; i<=m; i++) printf("  %+e\n",par[i]);
		printf("Number of integration steps performed: %d\n",nis);
	} else if (post == 4) {
		if (nbp == 0)
			printf("Minimization is started without break-points\n");
		else {
			printf("Minimization is started with weight =%3d\n"
				"The extra parameters are the observations:\n",
					weight);
			for (i=1; i<=nbp; i++) printf("  %5d\n",bp[i]);
		}
		printf("Starting values of the parameters:\n");
		for (i=1; i<=m; i++) printf("  %+e\n",par[i]);
		printf("Rel. tol. for the Eucl. norm of the res. vector:"
			" %e\nAbs. tol. for the Eucl. norm of the res. vector:"
			" %e\nRelative starting value of lambda: %e\n",
			in[3],in[4],in[6]);
	} else if (post == 1) {
		printf("\nStarting values of the parameters:\n");
		for (i=1; i<=m; i++) printf("  %+e\n",par[i]);
		printf("Numer of equations    : %2d\n"
			"Number of observations: %2d\n\n"
			"Machine precision                         : %6.2e\n"
			"Relative local error bound for integration: %6.2e\n"
			"Relative tolerance for residual           : %6.2e\n"
			"Absolute tolerance for residual           : %6.2e\n"
			"Maximum number of integrations to perform : %3.0f\n"
			"Relative starting value of lambda         : %6.2e\n"
			"Relative minimal steplength               : %6.2e\n",
			n,nobs,in[0],in[2],in[3],in[4],in[5],in[6],in[1]);
		if (nbp == 0)
			printf("There are no break-points\n");
		else {
			printf("Break-points are the observations:");
			for (i=1; i<=nbp; i++) printf("  %3d",bp[i]);
		}
		printf("\nThe alpha-point of the f-distribution: %5.2f\n",
				fa);
	} else if (post == 2) {
		if (out[1] == 0.0)
			printf("\nNormal termination of the process\n");
		else if (out[1] == 1.0)
			printf("Number of integrations allowed was exceeded\n");
		else if (out[1] == 2.0)
			printf("Minimal steplength was decreased four times\n");
		else if (out[1] == 3.0)
			printf("A call of deriv delivered false\n");
		else if (out[1] == 4.0)
			printf("A call of jacdfdy delivered false\n");
		else if (out[1] == 5.0)
			printf("A call of jacdfdp delivered false\n");
		else if (out[1] == 6.0)
			printf("Precision asked for may not be attained\n");
		if (nbp == 0)
			printf("Last integration was performed "
					"without break-points\n");
		else {
			printf("The process stopped with Break-points:");
			for (i=1; i<=nbp; i++) printf("  %3d",bp[i]);
		}
		printf("\nEucl. norm of the last residual vector : %e\n"
				"Eucl. norm of the first residual vector: %e\n"
				"Number of integrations performed       : %2.0f\n"
				"Last improvement of the euclidean norm : %e\n"
				"Condition number of J'*J               : %e\n"
				"Local error bound was exceeded (maxim.): %2.0f\n",
				out[2],out[3],out[4],out[6],out[7],out[5]);
		printf("\n   Parameters    Confidence Interval\n");
		for (i=1; i<=m; i++) {
			conf[i]=sqrt(m*fa*jtjinv[i][i]/(nobs-m))*out[2];
			printf("  %+e      %+e\n",par[i],conf[i]);
		}
		c = (nobs == m) ? 0.0 : out[2]*out[2]/(nobs-m);
		printf("\nCorrelation matrix         Covariance matrix\n");
		for (i=1; i<=m; i++) {
			for (j=1; j<=m; j++) {
				if (i == j)
					printf("                     ");
				if (i > j)
					printf("%9.3e  ",
						jtjinv[i][j]/sqrt(jtjinv[i][i]*jtjinv[j][j]));
				else
					printf("%9.3e  ",jtjinv[i][j]*c);
			}
			printf("\n");
		}
		printf("\nThe last residual vector\n  I      RES[I]\n");
		for (i=1; i<=nobs; i++) printf(" %2d   %+e\n",i,res[i]);
	}
	free_real_vector(conf,1);
}

int jacdfdp(int n, int m, real_t par[], real_t y[], real_t x,
				real_t **fp)
{
	real_t y2;

	y2=y[2];
	fp[1][1]=fp[1][3]=0.0;
	fp[1][2]=y2*exp(par[2]);
	fp[2][1]=exp(par[1])*(y[1]*(1.0-y2)-
				(exp(par[2])+exp(par[3]))*y2);
	fp[2][2] = -exp(par[1]+par[2])*y2;
	fp[2][3] = -exp(par[1]+par[3])*y2;
	return 1;
}

void data(int nobs, real_t tobs[], real_t obs[], int cobs[])
{
	int i;
	static real_t a[23]={0.0002, 0.0004, 0.0006, 0.0008, 0.001,
			0.0012, 0.0014, 0.0016, 0.0018, 0.002, 0.02, 0.04, 0.06,
			0.08, 0.1, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0};
	static real_t b[23]={0.1648, 0.2753, 0.3493, 0.3990, 0.4322,
			0.4545, 0.4695, 0.4795, 0.4862, 0.4907, 0.4999, 0.4998,
			0.4998, 0.4998, 0.4998, 0.4986, 0.4973, 0.4936, 0.4872,
			0.4808, 0.4743, 0.4677, 0.4610};

	tobs[0]=0.0;
	for (i=1; i<=nobs; i++) {
		tobs[i]=a[i-1];
		obs[i]=b[i-1];
		cobs[i]=2;
	}
	printf("\nThe observations were:\n"
		"  I   TOBS[I]   COBS[I]   OBS[I]\n");
	for (i=1; i<=nobs; i++)
		printf(" %2d   %7.4f     %1d       %6.4f\n",
				i,tobs[i],cobs[i],obs[i]);
}

void callystart(int n, int m, real_t par[], real_t y[], real_t ymax[])
{
	y[1]=ymax[1]=ymax[2]=1.0;
	y[2]=0.0;
}

int deriv(int n, int m, real_t par[], real_t y[], real_t x,
			real_t df[])
{
	real_t y2;

	y2=y[2];
	df[1] = -(1.0-y2)*y[1]+exp(par[2])*y2;
	df[2]=exp(par[1])*((1.0-y2)*y[1]-(exp(par[2])+exp(par[3]))*y2);
	return 1;
}

int jacdfdy(int n, int m, real_t par[], real_t y[], real_t x,
				real_t **fy)
{
	fy[1][1] = -1.0+y[2];
	fy[1][2]=exp(par[2])+y[1];
	fy[2][1]=exp(par[1])*(1.0-y[2]);
	fy[2][2] = -exp(par[1])*(exp(par[2])+exp(par[3])+y[1]);
	return 1;
}

void monitor(int post, int ncol, int nrow, real_t par[],
		real_t res[], int weight, int nis)
{
}

void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void peide(int, int, int, int *, real_t [], real_t [],
			int [], real_t **, real_t [], real_t [],
			int (*)(int,int,real_t [],real_t [],real_t,real_t []),
			int (*)(int,int,real_t [],real_t [],real_t,real_t **),
			int (*)(int,int,real_t [],real_t [],real_t,real_t **),
			void (*)(int,int,real_t [],real_t [],real_t[]),
			void (*)(int,real_t [],real_t [],int[]),
			void (*)(int,int,int,real_t [],real_t [],int,int));
	int m,n,nobs,nbp,bp[4];
	real_t fa,par[7],res[27],**jtjinv,in[7],out[8];

	jtjinv=allocate_real_matrix(1,3,1,3);
	printf("        E S C E P - problem\n");
	m=3;  n=2;  nobs=23;  nbp=3;
	par[1]=log(1600.0);  par[2]=log(0.8);  par[3]=log(1.2);
	in[0]=1.0e-14;  in[3]=in[4]=1.0e-4;  in[5]=50.0;  in[6]=1.0e-2;
	in[1]=1.0e-4;  in[2]=1.0e-5;
	bp[1]=17;  bp[2]=19;  bp[3]=21;
	fa=4.94;
	communication(1,fa,n,m,nobs,nbp,par,res,bp,jtjinv,in,out,0,0);
	peide(n,m,nobs,&nbp,par,res,bp,jtjinv,in,out,deriv,jacdfdy,
			jacdfdp,callystart,data,monitor);
	communication(2,fa,n,m,nobs,nbp,par,res,bp,jtjinv,in,out,0,0);
	free_real_matrix(jtjinv,1,3,1);
}

