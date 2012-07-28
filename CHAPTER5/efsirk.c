#include "../real.h"


void efsirk(real_t *x, real_t xe, int m, real_t y[],
			real_t *delta, void (*derivative)(int, real_t[], real_t *),
			void (*jacobian)(int, real_t **, real_t [], real_t *),
			real_t **j, int *n, real_t aeta, real_t reta, real_t hmin,
			real_t hmax, int linear,
			void (*output)(real_t, real_t, int, real_t [],
								real_t, real_t **, int))
{
	int *allocate_integer_vector(int, int);
	real_t *allocate_real_vector(int, int);
	real_t **allocate_real_matrix(int, int, int, int);
	void free_integer_vector(int *, int);
	void free_real_vector(real_t *, int);
	void free_real_matrix(real_t **, int, int, int);
	real_t vecvec(int, int, int, real_t [], real_t []);
	real_t matmat(int, int, int, int, real_t **, real_t **);
	real_t matvec(int, int, int, real_t **, real_t []);
	void gsselm(real_t **, int, real_t [], int [], int []);
	void solelm(real_t **, int, int [], int [], real_t []);
	int k,l,lin,*ri,*ci;
	real_t step,h,mu0,mu1,mu2,theta0,theta1,nu1,nu2,nu3,yk,fk,c1,c2,
			d,*f,*k0,*labda,**j1,aux[8],discr,eta,s,z1,z2,e,alpha1,a,b;

	ri=allocate_integer_vector(1,m);
	ci=allocate_integer_vector(1,m);
	f=allocate_real_vector(1,m);
	k0=allocate_real_vector(1,m);
	labda=allocate_real_vector(1,m);
	j1=allocate_real_matrix(1,m,1,m);

	aux[2]=FLT_EPSILON;
	aux[4]=8.0;
	for (k=1; k<=m; k++) f[k]=y[k];
	*n = 0;
	(*output)(*x,xe,m,y,*delta,j,*n);
	step=0.0;
	do {
		(*n)++;
		/* difference scheme */
		(*derivative)(m,f,delta);
		/* step size */
		if (linear)
			s=h=hmax;
		else
			if (*n == 1 || hmin == hmax)
				s=h=hmin;
			else {
				eta=aeta+reta*sqrt(vecvec(1,m,0,y,y));
				c1=nu3*step;
				for (k=1; k<=m; k++) labda[k] += c1*f[k]-y[k];
				discr=sqrt(vecvec(1,m,0,labda,labda));
				s=h=(eta/(0.75*(eta+discr))+0.33)*h;
				if (h < hmin)
					s=h=hmin;
				else
					if (h > hmax) s=h=hmax;
			}
		if ((*x)+s > xe) s=xe-(*x);
		lin=((step == s) && linear);
		step=s;
		if (!linear || *n == 1) (*jacobian)(m,j,y,delta);
		if (!lin) {
			/* coefficient */
			z1=step*(*delta);
			if (*n == 1) z2=z1+z1;
			if (fabs(z2-z1) > 1.0e-6*fabs(z1) || z2 > -1.0) {
				a=z1*z1+12.0;
				b=6.0*z1;
				if (fabs(z1) < 0.1)
					alpha1=(z1*z1/140.0-1.0)*z1/30.0;
				else if (z1 < 1.0e-14)
					alpha1=1.0/3.0;
				else if (z1 < -33.0)
					alpha1=(a+b)/(3.0*z1*(2.0+z1));
				else {
					e=((z1 < 230.0) ? exp(z1) : FLT_MAX);
					alpha1=((a-b)*e-a-b)/(((2.0-z1)*e-2.0-z1)*3.0*z1);
				}
				mu2=(1.0/3.0+alpha1)*0.25;
				mu1 = -(1.0+alpha1)*0.5;
				mu0=(6.0*mu1+2.0)/9.0;
				theta0=0.25;
				theta1=0.75;
				a=3.0*alpha1;
				nu3=(1.0+a)/(5.0-a)*0.5;
				a=nu3+nu3;
				nu1=0.5-a;
				nu2=(1.0+a)*0.75;
				z2=z1;
			}
			c1=step*mu1;
			d=step*step*mu2;
			for (k=1; k<=m; k++) {
				for (l=1; l<=m; l++)
					j1[k][l]=d*matmat(1,m,k,l,j,j)+c1*j[k][l];
				j1[k][k] += 1.0;
			}
			gsselm(j1,m,aux,ri,ci);
		}
		c1=step*step*mu0;
		d=step*2.0/3.0;
		for (k=1; k<=m; k++) {
			k0[k]=fk=f[k];
			labda[k]=d*fk+c1*matvec(1,m,k,j,f);
		}
		solelm(j1,m,ri,ci,labda);
		for (k=1; k<=m; k++) f[k]=y[k]+labda[k];
		(*derivative)(m,f,delta);
		c1=theta0*step;
		c2=theta1*step;
		d=nu1*step;
		for (k=1; k<=m; k++) {
			yk=y[k];
			fk=f[k];
			labda[k]=yk+d*fk+nu2*labda[k];
			y[k]=f[k]=yk+c1*k0[k]+c2*fk;
		}
		(*x) += step;
		(*output)(*x,xe,m,y,*delta,j,*n);
	} while (*x < xe);
	free_integer_vector(ri,1);
	free_integer_vector(ci,1);
	free_real_vector(f,1);
	free_real_vector(k0,1);
	free_real_vector(labda,1);
	free_real_matrix(j1,1,m,1);
}

