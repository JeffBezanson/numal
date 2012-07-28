#include "../real.h"


int multistep(real_t *x, real_t xend, real_t y[], real_t hmin, real_t hmax,
				real_t ymax[], real_t eps, int *first, real_t save[],
				void (*deriv)(real_t [], int, real_t, real_t[]),
				int (*available)(int, real_t, real_t [], real_t **),
				real_t **jacobian, int stiff, int n,
				void (*out)(real_t, int, int, real_t, real_t []))
{
	int *allocate_integer_vector(int, int);
	real_t *allocate_real_vector(int, int);
	real_t **allocate_real_matrix(int, int, int, int);
	void free_integer_vector(int *, int);
	void free_real_vector(real_t *, int);
	void free_real_matrix(real_t **, int, int, int);
	real_t matvec(int, int, int, real_t **, real_t []);
	void dec(real_t **, int, real_t [], int []);
	void sol(real_t **, int, int [], real_t []);
	void multistepreset(real_t [], real_t [], real_t *, real_t *,
			real_t *, real_t *, int *, real_t, real_t,
			real_t, real_t, int, int, int);
	void multisteporder(real_t [], real_t [], real_t *, real_t *,
				real_t *, real_t *, real_t *, int *, real_t, int, int);
	void multistepstep(int *, real_t *, real_t, real_t, real_t,
					real_t [], real_t,	real_t [], real_t [], real_t [],
					int, int, int, int);
	void multistepjacobian(int, real_t, real_t [], real_t,
			real_t [], real_t [], real_t [], real_t **,
			void (*)(real_t [], int, real_t, real_t[]),
			int (*)(int, real_t, real_t [], real_t **),
			int *, int *, int *);
	static real_t adams1[35] = {1.0, 1.0, 144.0, 4.0, 0.0, 0.5, 1.0,
		0.5, 576.0, 144.0, 1.0, 5.0/12.0, 1.0, 0.75, 1.0/6.0, 1436.0,
		576.0, 4.0, 0.375, 1.0, 11.0/12.0, 1.0/3.0, 1.0/24.0, 2844.0,
		1436.0, 1.0, 251.0/720.0, 1.0, 25.0/24.0, 35.0/72.0, 5.0/48.0,
		1.0/120.0, 0.0, 2844.0, 0.1};
	static real_t adams2[35] = {1.0, 1.0, 9.0, 4.0, 0.0, 2.0/3.0, 1.0,
		1.0/3.0, 36.0, 20.25, 1.0, 6.0/11.0, 1.0, 6.0/11.0, 1.0/11.0,
		84.028, 53.778, 0.25, 0.48, 1.0, 0.7, 0.2, 0.02, 156.25, 108.51,
		0.027778, 120.0/274.0, 1.0, 225.0/274.0, 85.0/274.0, 15.0/274.0,
		1.0/274.0, 0.0, 187.69, 0.0047361};
	static int adams,withjacobian,m,same,kold;
	static real_t xold,hold,a0,tolup,tol,toldwn,tolconv;
	int evaluate,evaluated,decompose,decomposed,conv,i,j,l,k,knew,
			fails,*p;
	real_t h,ch,chnew,error,dfi,c,a[6],*delta,*lastdelta,*df,**jac,
			aux[4],*fixy,*fixdy,*dy,ss,aa;

	p=allocate_integer_vector(1,n);
	delta=allocate_real_vector(1,n);
	lastdelta=allocate_real_vector(1,n);
	df=allocate_real_vector(1,n);
	fixy=allocate_real_vector(1,n);
	fixdy=allocate_real_vector(1,n);
	dy=allocate_real_vector(1,n);
	jac=allocate_real_matrix(1,n,1,n);

	if (*first) {
		*first = 0;
		m=n;
		save[-1]=save[-2]=save[-3]=0.0;
		(*out)(0.0,0,n,*x,y);
		adams=(!stiff);
		withjacobian=(!adams);
		if (withjacobian)
			multistepjacobian(n,*x,y,eps,fixy,fixdy,dy,jacobian,deriv,
							available,&evaluate,&decompose,&evaluated);
		if (adams)
			for (j=0; j<=34; j++) save[j-38]=adams1[j];
		else
			for (j=0; j<=34; j++) save[j-38]=adams2[j];
		NewStart:
		k=1;
		same=2;
		multisteporder(a,save,&tolup,&tol,&toldwn,&tolconv,
							&a0,&decompose,eps,k,n);
		(*deriv)(df,n,*x,y);
		if (!withjacobian)
			h=hmin;
		else {
			ss=FLT_MIN;
			for (i=1; i<=n; i++) {
				aa=matvec(1,n,i,jacobian,df)/ymax[i];
				ss += aa*aa;
			}
			h=sqrt(2.0*eps/sqrt(ss));
		}
		if (h > hmax)
			h=hmax;
		else
			if (h < hmin) h=hmin;
		xold=(*x);
		hold=h;
		kold=k;
		ch=1.0;
		for (i=1; i<=n; i++) {
			save[i]=y[i];
			save[m+i]=y[m+i]=df[i]*h;
		}
		(*out)(0.0,0,n,*x,y);
	} else {
		withjacobian=(!adams);
		ch=1.0;
		k=kold;
		multistepreset(y,save,x,&ch,&c,&h,&decomposed,hmin,hmax,
							hold,xold,m,k,n);
		multisteporder(a,save,&tolup,&tol,&toldwn,&tolconv,
							&a0,&decompose,eps,k,n);
		decompose=withjacobian;
	}
	fails=0;
	while (*x < xend) {
		if ((*x)+h <= xend)
			*x += h;
		else {
			h=xend-(*x);
			*x = xend;
			ch=h/hold;
			c=1.0;
			for (j=m; j<=k*m; j+=m) {
				c *= ch;
				for (i=j+1; i<=j+n; i++) y[i] *= c;
			}
			same = ((same < 3) ? 3 : same+1);
		}
		/* prediction */
		for (l=1; l<=n; l++) {
			for (i=l; i<=(k-1)*m+l; i+=m)
				for (j=(k-1)*m+l; j>=i; j-=m) y[j] += y[j+m];
			delta[l]=0.0;
		}
		evaluated=0;
		/* correction and estimation local error */
		for (l=1; l<=3; l++) {
			(*deriv)(df,n,*x,y);
			for (i=1; i<=n; i++) df[i]=df[i]*h-y[m+i];
			if (withjacobian) {
				if (evaluate)
					multistepjacobian(n,*x,y,eps,fixy,fixdy,dy,jacobian,
						deriv,available,&evaluate,&decompose,&evaluated);
				if (decompose) {
					/* decompose jacobian */
					decompose=0;
					decomposed=1;
					c = -a0*h;
					for (j=1; j<=n; j++) {
						for (i=1; i<=n; i++) jac[i][j]=jacobian[i][j]*c;
						jac[j][j] += 1.0;
					}
					aux[2]=FLT_EPSILON;
					dec(jac,n,aux,p);
				}
				sol(jac,n,p,df);
			}
			conv=1;
			for (i=1; i<=n; i++) {
				dfi=df[i];
				y[i] += a0*dfi;
				y[m+i] += dfi;
				delta[i] += dfi;
				conv=(conv && (fabs(dfi) < tolconv*ymax[i]));
			}
			if (conv) {
				ss=FLT_MIN;
				for (i=1; i<=n; i++) {
					aa=delta[i]/ymax[i];
					ss += aa*aa;
				}
				error=ss;
				break;
			}
		}
		/* acceptance or rejection */
		if (!conv) {
			if (!withjacobian) {
				evaluate=withjacobian=((same >= k) || (h < 1.1*hmin));
				if (!withjacobian) ch /= 4.0;
			}
			else if (!decomposed) decompose=1;
			else if (!evaluated) evaluate=1;
			else if (h > 1.1*hmin) ch /= 4.0;
			else if (adams) goto TryCurtiss;
			else {
				save[-1]=1.0;
				break;
			}
			multistepreset(y,save,x,&ch,&c,&h,&decomposed,hmin,hmax,
								hold,xold,m,k,n);
		} else
			if (error > tol) {
				fails++;
				if (h > 1.1*hmin) {
					if (fails > 2) {
						if (adams) {
							adams=0;
							for (j=0; j<=34; j++) save[j-38]=adams2[j];
						}
						kold=0;
						multistepreset(y,save,x,&ch,&c,&h,&decomposed,
											hmin,hmax,hold,xold,m,k,n);
						goto NewStart;
					} else {
						multistepstep(&knew,&chnew,tolup,tol,toldwn,delta,
									error,lastdelta,y,ymax,fails,m,k,n);
						if (knew != k) {
							k=knew;
							multisteporder(a,save,&tolup,&tol,&toldwn,
										&tolconv,&a0,&decompose,eps,k,n);
						}
						ch *= chnew;
						multistepreset(y,save,x,&ch,&c,&h,&decomposed,
											hmin,hmax,hold,xold,m,k,n);
					}
				} else {
					if (adams) {
						TryCurtiss:
						adams=0;
						for (j=0; j<=34; j++) save[j-38]=adams2[j];
					} else
						if (k == 1) {
							/* violate eps criterion */
							c=eps*sqrt(error/tol);
							if (c > save[-3]) save[-3]=c;
							save[-2] += 1.0;
							same=4;
							goto ErrorTestOk;
						}
					k=kold=1;
					multistepreset(y,save,x,&ch,&c,&h,&decomposed,
										hmin,hmax,hold,xold,m,k,n);
					multisteporder(a,save,&tolup,&tol,&toldwn,
										&tolconv,&a0,&decompose,eps,k,n);
					same=2;
				}
			} else {
				ErrorTestOk:
				fails=0;
				for (i=1; i<=n; i++) {
					c=delta[i];
					for (l=2; l<=k; l++) y[l*m+i] += a[l]*c;
					if (fabs(y[i]) > ymax[i]) ymax[i]=fabs(y[i]);
				}
				same--;
				if (same == 1)
					for (i=1; i<=n; i++) lastdelta[i]=delta[i];
				else
					if (same == 0) {
						multistepstep(&knew,&chnew,tolup,tol,toldwn,delta,
									error,lastdelta,y,ymax,fails,m,k,n);
						if (chnew > 1.1) {
							decomposed=0;
							if (k != knew) {
								if (knew > k)
									for (i=1; i<=n; i++)
										y[knew*m+i]=delta[i]*a[k]/knew;
								k=knew;
								multisteporder(a,save,&tolup,&tol,&toldwn,
											&tolconv,&a0,&decompose,eps,k,n);
							}
							same=k+1;
							if (chnew*h > hmax) chnew=hmax/h;
							h *= chnew;
							c=1.0;
							for (j=m; j<=k*m; j+=m) {
								c *= chnew;
								for (i=j+1; i<=j+n; i++) y[i] *= c;
							}
						} else
							same=10;
					}
				if (*x != xend) {
					xold=(*x);
					hold=h;
					kold=k;
					ch=1.0;
					for (i=k*m+n; i>=1; i--) save[i]=y[i];
					(*out)(h,k,n,*x,y);
				}
			}
	}
	save[0]=(adams ? 0.0 : 1.0);
	free_integer_vector(p,1);
	free_real_vector(delta,1);
	free_real_vector(lastdelta,1);
	free_real_vector(df,1);
	free_real_vector(fixy,1);
	free_real_vector(fixdy,1);
	free_real_vector(dy,1);
	free_real_matrix(jac,1,n,1);
	return ((save[-1] == 0.0) && (save[-2] == 0.0));
}

void multistepreset(real_t y[], real_t save[], real_t *x, real_t *ch,
			real_t *c, real_t *h, int *decomposed, real_t hmin, real_t hmax,
			real_t hold, real_t xold, int m, int k, int n)
{
	/* this function is internally used by MULTISTEP */

	int i,j;

	if (*ch < hmin/hold)
		*ch = hmin/hold;
	else
		if (*ch > hmax/hold) *ch = hmax/hold;
	*x = xold;
	*h = hold*(*ch);
	*c = 1.0;
	for (j=0; j<=k*m; j+=m) {
		for (i=1; i<=n; i++) y[j+i]=save[j+i]*(*c);
		(*c) *= (*ch);
	}
	*decomposed = 0;
}

void multisteporder(real_t a[], real_t save[], real_t *tolup,
			real_t *tol, real_t *toldwn, real_t *tolconv, real_t *a0,
			int *decompose, real_t eps, int k, int n)
{
	/* this function is internally used by MULTISTEP */

	int i,j;
	real_t c;

	c=eps*eps;
	j=(k-1)*(k+8)/2-38;
	for (i=0; i<=k; i++) a[i]=save[i+j];
	*tolup = c*save[j+k+1];
	*tol = c*save[j+k+2];
	*toldwn = c*save[j+k+3];
	*tolconv = eps/(2*n*(k+2));
	*a0 = a[0];
	*decompose = 1;
}

void multistepstep(int *knew, real_t *chnew, real_t tolup, real_t tol,
					real_t toldwn, real_t delta[], real_t error,
					real_t lastdelta[], real_t y[], real_t ymax[],
					int fails, int m, int k, int n)
{
	/* this function is internally used by MULTISTEP */

	int i;
	real_t a1,a2,a3,aa,ss;

	if (k <= 1)
		a1=0.0;
	else {
		ss=FLT_MIN;
		for (i=1; i<=n; i++) {
			aa=y[k*m+i]/ymax[i];
			ss += aa*aa;
		}
		a1=0.75*pow(toldwn/ss,0.5/k);
	}
	a2=0.80*pow(tol/error,0.5/(k+1));
	if (k >= 5 || fails != 0)
		a3=0.0;
	else {
		ss=FLT_MIN;
		for (i=1; i<=n; i++) {
			aa=(delta[i]-lastdelta[i])/ymax[i];
			ss += aa*aa;
		}
		a3=0.70*pow(tolup/ss,0.5/(k+2));
	}
	if (a1 > a2 && a1 > a3) {
		*knew = k-1;
		*chnew = a1;
	} else
		if (a2 > a3) {
			*knew = k;
			*chnew = a2;
		} else {
			*knew = k+1;
			*chnew = a3;
		}
}

void multistepjacobian(int n, real_t x, real_t y[], real_t eps,
			real_t fixy[], real_t fixdy[], real_t dy[], real_t **jacobian,
			void (*deriv)(real_t [], int, real_t, real_t[]),
			int (*available)(int, real_t, real_t [], real_t **),
			int *evaluate, int *decompose, int *evaluated)
{
	/* this function is internally used by MULTISTEP */

	int i,j;
	real_t d;

	*evaluate = 0;
	*decompose = *evaluated = 1;
	if (!(*available)(n,x,y,jacobian)) {
		for (i=1; i<=n; i++) fixy[i]=y[i];
		(*deriv)(fixdy,n,x,y);
		for (j=1; j<=n; j++) {
			d=((eps > fabs(fixy[j])) ? eps*eps : eps*fabs(fixy[j]));
			y[j] += d;
			(*deriv)(dy,n,x,y);
			for (i=1; i<=n; i++) jacobian[i][j]=(dy[i]-fixdy[i])/d;
			y[j]=fixy[j];
		}
	}
}

