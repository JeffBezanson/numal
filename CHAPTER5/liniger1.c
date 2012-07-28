#include "../real.h"


int *allocate_integer_vector(int, int);
real_t *allocate_real_vector(int, int);
real_t **allocate_real_matrix(int, int, int, int);
void free_integer_vector(int *, int);
void free_real_vector(real_t *, int);
void free_real_matrix(real_t **, int, int, int);
void inivec(int, int, real_t [], real_t);
void mulvec(int, int, int, real_t [], real_t [], real_t);
void mulrow(int, int, int, int, real_t **, real_t **, real_t);
void dupvec(int, int, int, real_t [], real_t []);
real_t vecvec(int, int, int, real_t [], real_t []);
real_t matvec(int, int, int, real_t **, real_t []);
void elmvec(int, int, int, real_t [], real_t [], real_t);
void dec(real_t **, int, real_t [], int []);
void sol(real_t **, int, int [], real_t []);
void liniger1vscoef(int, real_t **, real_t **, real_t [], int [],
					real_t, real_t, real_t *, real_t *, real_t *, real_t *);

void liniger1vs(real_t *x, real_t xe, int m, real_t y[], real_t *sigma,
			void (*derivative)(int, real_t[], real_t *), real_t **j,
			void (*jacobian)(int, real_t **, real_t [], real_t *),
			int itmax, real_t hmin, real_t hmax, real_t aeta, real_t reta,
			real_t info[],
			void (*output)(real_t, real_t, int, real_t [], real_t,
								real_t **, real_t []))
{
	int i,st,lastjac,last,first,evaljac,evalcoef,*pi,k,q;
	real_t h,hnew,mu,mu1,beta,p,e,e1,eta,eta1,discr,*dy,*yl,*yr,*f,
			**a,aux[4],hl,b;

	pi=allocate_integer_vector(1,m);
	dy=allocate_real_vector(1,m);
	yl=allocate_real_vector(1,m);
	yr=allocate_real_vector(1,m);
	f=allocate_real_vector(1,m);
	a=allocate_real_matrix(1,m,1,m);

	first=evaljac=1;
	last=evalcoef=0;
	inivec(1,9,info,0.0);
	eta=reta*sqrt(vecvec(1,m,0,y,y))+aeta;
	eta1=eta/sqrt(fabs(reta));
	dupvec(1,m,0,f,y);
	(*derivative)(m,f,sigma);
	info[2]=1.0;
	st=1;
	do {
		/* step size */
		if (eta < 0.0) {
			hl=h;
			h=hnew=hmax;
			info[5] += 1.0;
			if (1.1*hnew > xe-(*x)) {
				last=1;
				h=hnew=xe-(*x);
			}
			evalcoef=(h != hl);
		} else if (first) {
			h=hnew=hmin;
			first=0;
			info[4] += 1.0;
		} else {
			b=discr/eta;
			hl=h;
			if (b < 0.01) b=0.01;
			hnew = (b > 0.0) ? h*pow(b,-1.0/p) : hmax;
			if (hnew < hmin) {
				hnew=hmin;
				info[4] += 1.0;
			} else
				if (hnew > hmax)  {
					hnew=hmax;
					info[5] += 1.0;
				}
			if (1.1*hnew >= xe-(*x)) {
				last=1;
				h=hnew=xe-(*x);
			} else
				if (fabs(h/hnew-1.0) > 0.1) h=hnew;
			evalcoef=(h !=hl);
		}
		info[1] += 1.0;
		if (evaljac) {
			(*jacobian)(m,j,y,sigma);
			info[3] += 1.0;
			h=hnew;
			liniger1vscoef(m,a,j,aux,pi,h,*sigma,&mu,&mu1,&beta,&p);
			evaljac=0;
			lastjac=st;
		} else
			if (evalcoef)
				liniger1vscoef(m,a,j,aux,pi,h,*sigma,&mu,&mu1,&beta,&p);
		i=1;
		do {
			/* iteration */
			if (reta < 0.0) {
				if (i == 1) {
					mulvec(1,m,0,dy,f,h);
					for (k=1; k<=m; k++) yl[k]=y[k]+mu*dy[k];
					sol(a,m,pi,dy);
					e=1.0;
				} else {
					for (k=1; k<=m; k++) dy[k]=yl[k]-y[k]+mu1*f[k];
					if (e*sqrt(vecvec(1,m,0,y,y)) > e1*e1) {
						evaljac=(i >= 3);
						if (i > 3) {
							info[3] += 1.0;
							(*jacobian)(m,j,y,sigma);
							for (q=1; q<=m; q++) {
								mulrow(1,m,q,q,a,j,-mu1);
								a[q][q] += 1.0;
							}
							aux[2]=0.0;
							dec(a,m,aux,pi);
						}
					}
					sol(a,m,pi,dy);
				}
				e1=e;
				e=sqrt(vecvec(1,m,0,dy,dy));
				elmvec(1,m,0,y,dy,1.0);
				eta=sqrt(vecvec(1,m,0,y,y))*reta+aeta;
				discr=0.0;
				dupvec(1,m,0,f,y);
				(*derivative)(m,f,sigma);
				info[2] += 1.0;
			} else {
				if (i == 1) {
					/* linearity */
					for (k=1; k<=m; k++) dy[k]=y[k]-mu1*f[k];
					sol(a,m,pi,dy);
					elmvec(1,m,0,dy,y,-1.0);
					e=sqrt(vecvec(1,m,0,dy,dy));
					if (e*(st-lastjac) > eta) {
						(*jacobian)(m,j,y,sigma);
						lastjac=st;
						info[3] += 1.0;
						h=hnew;
						liniger1vscoef(m,a,j,aux,pi,h,*sigma,&mu,&mu1,
											&beta,&p);
						/* linearity */
						for (k=1; k<=m; k++) dy[k]=y[k]-mu1*f[k];
						sol(a,m,pi,dy);
						elmvec(1,m,0,dy,y,-1.0);
						e=sqrt(vecvec(1,m,0,dy,dy));
					}
					evaljac=(e*(st+1-lastjac) > eta);
					mulvec(1,m,0,dy,f,h);
					for (k=1; k<=m; k++) yl[k]=y[k]+mu*dy[k];
					sol(a,m,pi,dy);
					for (k=1; k<=m; k++) yr[k]=h*beta*matvec(1,m,k,j,dy);
					sol(a,m,pi,yr);
					elmvec(1,m,0,yr,dy,1.0);
				} else {
					for (k=1; k<=m; k++) dy[k]=yl[k]-y[k]+mu1*f[k];
					if (e > eta1 && discr > eta1) {
						info[3] += 1.0;
						(*jacobian)(m,j,y,sigma);
						for (q=1; q<=m; q++) {
							mulrow(1,m,q,q,a,j,-mu1);
							a[q][q] += 1.0;
						}
						aux[2]=0.0;
						dec(a,m,aux,pi);
					}
					sol(a,m,pi,dy);
					e=sqrt(vecvec(1,m,0,dy,dy));
				}
				elmvec(1,m,0,y,dy,1.0);
				eta=sqrt(vecvec(1,m,0,y,y))*reta+aeta;
				eta1=eta/sqrt(reta);
				dupvec(1,m,0,f,y);
				(*derivative)(m,f,sigma);
				info[2] += 1.0;
				for (k=1; k<=m; k++) dy[k]=yr[k]-h*f[k];
				discr=sqrt(vecvec(1,m,0,dy,dy))/2.0;
			}
			if (i > info[6]) info[6]=i;
			i++;
		} while (e > fabs(eta) && discr > 1.3*eta && i <= itmax);
		info[7]=eta;
		info[8]=discr;
		(*x) += h;
		if (discr > info[9]) info[9]=discr;
		(*output)(*x,xe,m,y,*sigma,j,info);
		st++;
	} while (!last);
	free_integer_vector(pi,1);
	free_real_vector(dy,1);
	free_real_vector(yl,1);
	free_real_vector(yr,1);
	free_real_vector(f,1);
	free_real_matrix(a,1,m,1);
}

void liniger1vscoef(int m, real_t **a, real_t **j, real_t aux[],
				int pi[], real_t h, real_t sigma, real_t *mu,
				real_t *mu1, real_t *beta, real_t *p)
{
	/* this function is internally used by LINIGER1VS */

	int q;
	real_t b,e;

	b=fabs(h*sigma);
	if (b > 40.0) {
		*mu = 1.0/b;
		*beta = 1.0;
		*p = 2.0+2.0/(b-2.0);
	} else if (b < 0.04) {
		e=b*b/30.0;
		*p = 3.0-e;
		*mu = 0.5-b/12.0*(1.0-e/2.0);
		*beta = 0.5+b/6.0*(1.0-e);
	} else {
		e=exp(b)-1.0;
		*mu = 1.0/b-1.0/e;
		*beta = (1.0-b/e)*(1.0+1.0/e);
		*p = ((*beta)-(*mu))/(0.5-(*mu));
	}
	*mu1 = h*(1.0-(*mu));
	for (q=1; q<=m; q++) {
		mulrow(1,m,q,q,a,j,-(*mu1));
		a[q][q] += 1.0;
	}
	aux[2]=0.0;
	dec(a,m,aux,pi);
}
