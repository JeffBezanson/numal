#include "../real.h"


void eferk(real_t *x, real_t xe, int m, real_t y[],
			real_t *sigma, real_t phi,
			void (*derivative)(int, real_t[]),
			real_t **j,
			void (*jacobian)(int, real_t **, real_t [], real_t *),
			int *k, int l, int aut, real_t aeta, real_t reta, real_t hmin,
			real_t hmax, int linear,
			void (*output)(real_t, real_t, int, real_t [], real_t **, int))
{
	int *allocate_integer_vector(int, int);
	real_t *allocate_real_vector(int, int);
	real_t **allocate_real_matrix(int, int, int, int);
	void free_integer_vector(int *, int);
	void free_real_vector(real_t *, int);
	void free_real_matrix(real_t **, int, int, int);
	real_t vecvec(int, int, int, real_t [], real_t []);
	real_t matvec(int, int, int, real_t **, real_t []);
	void dec(real_t **, int, real_t [], int []);
	void sol(real_t **, int, int [], real_t []);
	int m1,i,change,last,*p,c1,q;
	real_t h,b,b0,phi0,cosphi,sinphi,eta,discr,fac,pi,*beta,*betha,
			*betac,*k0,*d,*d1,*d2,**a,aux[4],s,cos2phi,sina,e,zi,
			c2,cosiphi,siniphi,cosphil,*dd,emin1,b1,b2,a0,a1,a2,a3,
			c,ddd,betai,bethai;

	p=allocate_integer_vector(1,l);
	beta=allocate_real_vector(0,l);
	betha=allocate_real_vector(0,l);
	betac=allocate_real_vector(0,l+3);
	k0=allocate_real_vector(1,m);
	d=allocate_real_vector(1,m);
	d1=allocate_real_vector(1,m);
	d2=allocate_real_vector(1,m);
	dd=allocate_real_vector(1,l);
	a=allocate_real_matrix(1,l,1,l);

	b0 = phi0 = -1.0;
	pi=4.0*atan(1.0);
	betac[l]=betac[l+1]=betac[l+2]=betac[l+3]=0.0;
	beta[0]=1.0/6.0;
	betha[0]=0.5;
	fac=1.0;
	for (i=2; i<=l-1; i++) fac *= i;
	m1=(aut ? m : m-1);
	*k = 0;
	last=0;
	do {
		for (i=1; i<=m; i++) d[i]=y[i];
		(*derivative)(m,d);
		if (!linear || *k == 0) (*jacobian)(m,j,y,sigma);
		/* step size */
		eta=aeta+reta*sqrt(vecvec(1,m1,0,y,y));
		if (*k == 0) {
			discr=sqrt(vecvec(1,m1,0,d,d));
			h=eta/discr;
		} else {
			s=0.0;
			for (i=1; i<=m1; i++) {
				ddd=d[i]-d2[i];
				s += ddd*ddd;
			}
			discr=h*sqrt(s)/eta;
			h *= (linear ? 4.0/(4.0+discr)+0.5 :
								4.0/(3.0+discr)+1.0/3.0);
		}
		if (h < hmin) h=hmin;
		if (h > hmax) h=hmax;
		b=fabs(h*(*sigma));
		change=(fabs(1.0-b/b0) > 0.05 || phi != phi0);
		if (1.1*h >= xe-(*x)) {
			change=last=1;
			h=xe-(*x);
		}
		if (!change) h=h*b0/b;
		if (change) {
			/* coefficient */
			b0=b=fabs(h*(*sigma));
			if (b >= 0.1) {
				if (phi != pi && l == 2 || fabs(phi-pi) > 0.01) {
					/* solution of complex equations */
					if (l == 2) {
						phi0=phi;
						cosphi=cos(phi0);
						sinphi=sin(phi0);
						e=exp(b*cosphi);
						zi=b*sinphi-3.0*phi0;
						sina=((fabs(sinphi) < 1.0e-6) ? -e*(b+3.0) :
								e*sin(zi)/sinphi);
						cos2phi=2.0*cosphi*cosphi-1.0;
						betha[2]=(0.5+(2.0*cosphi+(1.0+2.0*cos2phi+
										sina)/b)/b)/b/b;
						sina=((fabs(sinphi) < 1.0e-6) ? e*(b+4.0) :
								sina*cosphi-e*cos(zi));
						betha[1] = -(cosphi+(1.0+2.0*cos2phi+
										(4.0*cosphi*cos2phi+sina)/b)/b)/b;
						beta[1]=betha[2]+2.0*cosphi*(betha[1]-1.0/6.0)/b;
						beta[2]=(1.0/6.0-betha[1])/b/b;
					} else {
						if (phi0 != phi) {
							/* elements of matrix */
							phi0=phi;
							cosphi=cos(phi0);
							sinphi=sin(phi0);
							cosiphi=1.0;
							siniphi=0.0;
							for (i=0; i<=l-1; i++) {
								c1=4+i;
								c2=1.0;
								for (q=l-1; q>=1; q-=2) {
									a[q][l-i]=c2*cosiphi;
									a[q+1][l-i]=c2*siniphi;
									c2 *= c1;
									c1--;
								}
								cosphil=cosiphi*cosphi-siniphi*sinphi;
								siniphi=cosiphi*sinphi+siniphi*cosphi;
								cosiphi=cosphil;
							}
							aux[2]=0.0;
							dec(a,l,aux,p);
						}
						/* right hand side */
						e=exp(b*cosphi);
						zi=b*sinphi-4.0*phi0;
						cosiphi=e*cos(zi);
						siniphi=e*sin(zi);
						zi=1.0/b/b/b;
						for (q=l; q>=2; q-=2) {
							dd[q]=zi*siniphi;
							dd[q-1]=zi*cosiphi;
							cosphil=cosiphi*cosphi-siniphi*sinphi;
							siniphi=cosiphi*sinphi+siniphi*cosphi;
							cosiphi=cosphil;
							zi *= b;
						}
						siniphi=2.0*sinphi*cosphi;
						cosiphi=2.0*cosphi*cosphi-1.0;
						cosphil=cosphi*(2.0*cosiphi-1.0);
						dd[l] += sinphi*(1.0/6.0+(cosphi+(1.0+2.0*
									cosiphi*(1.0+2.0*cosphi/b))/b)/b);
						dd[l-1] -= cosphi/6.0+(0.5*cosiphi+(cosphil+
										(2.0*cosiphi*cosiphi-1.0)/b)/b)/b;
						dd[l-2] += sinphi*(0.5+(2.0*cosphi+
										(2.0*cosiphi+1.0)/b)/b);
						dd[l-3] -= 0.5*cosphi-(cosiphi+cosphil/b)/b;
						if (l >= 5) {
							dd[l-4] += sinphi+siniphi/b;
							dd[l-5] -= cosphi+cosiphi/b;
							if (l >= 7) {
								dd[l-6] += sinphi;
								dd[l-7] -= cosphi;
							}
						}
						sol(a,l,p,dd);
						zi=1.0/b;
						for (i=1; i<=l; i++) {
							beta[i]=dd[l+1-i]*zi;
							betha[i]=(i+3)*beta[i];
							zi /= b;
						}
					}
				} else {
					/* form beta */
					if (l == 1) {
						betha[1]=(0.5-(1.0-(1.0-exp(-b))/b)/b)/b;
						beta[1]=(1.0/6.0-betha[1])/b;
					} else if (l == 2) {
						e=exp(-b);
						emin1=e-1.0;
						betha[1]=(1.0-(3.0+e+4.0*emin1/b)/b)/b;
						betha[2]=(0.5-(2.0+e+3.0*emin1/b)/b)/b/b;
						beta[2]=(1.0/6.0-betha[1])/b/b;
						beta[1]=(1.0/3.0-(1.5-(4.0+e+5.0*emin1/b)/b)/b)/b;
					} else {
						betac[l-1]=c=ddd=exp(-b)/fac;
						for (i=l-1; i>=1; i--) {
							c=i*b*c/(l-i);
							betac[i-1]=ddd=ddd*i+c;
						}
						b2=0.5-betac[2];
						b1=(1.0-betac[1])*(l+1)/b;
						b0=(1.0-betac[0])*(l+2)*(l+1)*0.5/b/b;
						a3=1.0/6.0-betac[3];
						a2=b2*(l+1)/b;
						a1=b1*(l+2)*0.5/b;
						a0=b0*(l+3)/3.0/b;
						ddd=l/b;
						for (i=1; i<=l; i++) {
							beta[i]=(a3/i-a2/(i+1)+a1/(i+2)-a0/(i+3))*ddd+
										betac[i+3];
							betha[i]=(b2/i-b1/(i+1)+b0/(i+2))*ddd+
										betac[i+2];
							ddd=ddd*(l-i)/i/b;
						}
					}
				}
			} else
				for (i=1; i<=l; i++) {
					betha[i]=beta[i-1];
					beta[i]=beta[i-1]/(i+3);
				}
		}
		(*output)(*x,xe,m,y,j,*k);
		/* difference scheme */
		if (m1 < m) {
			d2[m]=1.0;
			k0[m]=y[m]+2.0*h/3.0;
			y[m] += 0.25*h;
		}
		for (q=1; q<=m1; q++) {
			k0[q]=y[q]+2.0*h/3.0*d[q];
			y[q] += 0.25*h*d[q];
			d1[q]=h*matvec(1,m,q,j,d);
			d2[q]=d1[q]+d[q];
		}
		for (i=0; i<=l; i++) {
			betai=4.0*beta[i]/3.0;
			bethai=betha[i];
			for (q=1; q<=m1; q++) d[q]=h*d1[q];
			for (q=1; q<=m1; q++) {
				k0[q] += betai*d[q];
				d1[q]=matvec(1,m1,q,j,d);
				d2[q] += bethai*d1[q];
			}
		}
		(*derivative)(m,k0);
		for (q=1; q<=m; q++) y[q] += 0.75*h*k0[q];
		(*k)++;
		(*x) += h;
	} while (!last);
	(*output)(*x,xe,m,y,j,*k);
	free_integer_vector(p,1);
	free_real_vector(beta,0);
	free_real_vector(betha,0);
	free_real_vector(betac,0);
	free_real_vector(k0,1);
	free_real_vector(d,1);
	free_real_vector(d1,1);
	free_real_vector(d2,1);
	free_real_vector(dd,1);
	free_real_matrix(a,1,l,1);
}

