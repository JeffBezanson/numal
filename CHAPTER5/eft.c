#include "../real.h"


void eft(real_t *t, real_t te, int m0, int m, real_t u[],
			real_t (*sigma)(real_t, int, int), real_t phi,
			real_t (*diameter)(real_t, int, int),
			void (*derivative)(real_t, int, int, int, real_t []),
			int *k, real_t alfa, int norm,
			real_t (*aeta)(real_t, int, int),
			real_t (*reta)(real_t, int, int),
			real_t *eta, real_t *rho, real_t hmin, real_t *hstart,
			void (*out)(real_t, real_t, int, int, real_t [],
							int, real_t, real_t))
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void inivec(int, int, real_t [], real_t);
	void dupvec(int, int, int, real_t [], real_t []);
	real_t vecvec(int, int, int, real_t [], real_t []);
	void elmvec(int, int, int, real_t [], real_t [], real_t);
	int kl,last,start,j,i,ext,extrapolate;
	real_t q,ec0,ec1,ec2,h,hi,h0,h1,h2,betan,t2,sigmal,phil,*c,
			*ro,beta[4],betha[4],s,x,hacc,hstab,hcr,hmax,a,b,cc,
			b1,b2,bb,e,beta2,beta3,
			c0,fc,b0,fb,a0,fa,d0,fd,fdb,fda,w,mb,tol,mm,p0,q0;

	c=allocate_real_vector(m0,m);
	ro=allocate_real_vector(m0,m);

	start=1;
	last=0;
	dupvec(m0,m,0,c,u);
	(*derivative)(*t,m0,m,1,c);
	if (*k == 0) {
		/* local error bound */
		s=0.0;
		if (norm == 1)
			for (j=m0; j<=m; j++) {
				x=fabs(u[j]);
				if (x > s) s=x;
			}
		else
			s=sqrt(vecvec(m0,m,0,u,u));
		*eta = (*aeta)(*t,m0,m)+(*reta)(*t,m0,m)*s;
		s=0.0;
		if (norm == 1)
			for (j=m0; j<=m; j++) {
				x=fabs(c[j]);
				if (x > s) s=x;
			}
		else
			s=sqrt(vecvec(m0,m,0,c,c));
		*hstart = (*eta)/s;
	}
	do {
		/* difference scheme */
		hi=1.0;
		sigmal=(*sigma)(*t,m0,m);
		phil=phi;
		/* step size */
		if (!start) {
			/* local error bound */
			s=0.0;
			if (norm == 1)
				for (j=m0; j<=m; j++) {
					x=fabs(u[j]);
					if (x > s) s=x;
				}
			else
				s=sqrt(vecvec(m0,m,0,u,u));
			*eta = (*aeta)(*t,m0,m)+(*reta)(*t,m0,m)*s;
		}
		if (start) {
			h1=h2=hacc=(*hstart);
			ec2=ec1=1.0;
			kl=1;
			start=0;
		} else if (kl < 3) {
			hacc=pow((*eta)/(*rho),1.0/q)*h2;
			if (hacc > 10.0*h2)
				hacc=10.0*h2;
			else
				kl++;
		} else {
			a=(h0*(ec2-ec1)-h1*(ec1-ec0))/(h2*h0-h1*h1);
			h=h2*((*eta < *rho) ? pow((*eta)/(*rho),1.0/q) : alfa);
			if (a > 0.0) {
				b=(ec2-ec1-a*(h2-h1))/h1;
				cc=ec2-a*h2-b*t2;
				hacc=0.0;
				hmax=h;
				/* find zero */
				b0=hacc;
				fb=pow(hacc,q)*(a*hacc+b*(*t)+cc)-(*eta);
				a0=hacc=h;
				fa=pow(hacc,q)*(a*hacc+b*(*t)+cc)-(*eta);
				c0=a0;
				fc=fa;
				ext=0;
				extrapolate=1;
				while (extrapolate) {
					if (fabs(fc) < fabs(fb)) {
						if (c0 != a0) {
							d0=a0;
							fd=fa;
						}
						a0=b0;
						fa=fb;
						b0=hacc=c0;
						fb=fc;
						c0=a0;
						fc=fa;
					}
					tol=1.0e-3*h2;
					mm=(c0+b0)*0.5;
					mb=mm-b0;
					if (fabs(mb) > tol) {
						if (ext > 2)
							w=mb;
						else {
							if (mb == 0.0)
								tol=0.0;
							else
								if (mb < 0.0) tol = -tol;
							p0=(b0-a0)*fb;
							if (ext <= 1)
								q0=fa-fb;
							else {
								fdb=(fd-fb)/(d0-b0);
								fda=(fd-fa)/(d0-a0);
								p0 *= fda;
								q0=fdb*fa-fda*fb;
							}
							if (p0 < 0.0) {
								p0 = -p0;
								q0 = -q0;
							}
							w=(p0<FLT_MIN || p0<=q0*tol) ? tol :
										((p0<mb*q0) ? p0/q0 : mb);
						}
						d0=a0;
						fd=fa;
						a0=b0;
						fa=fb;
						hacc = b0 += w;
						fb=pow(hacc,q)*(a*hacc+b*(*t)+cc)-(*eta);
						if ((fc >= 0.0) ? (fb >= 0.0) : (fb <= 0.0)) {
							c0=a0;
							fc=fa;
							ext=0;
						} else
							ext = (w == mb) ? 0 : ext+1;
					} else
						break;
				}
				h=c0;
				if (!((fc >= 0.0) ? (fb <= 0.0) : (fb >= 0.0)))
					hacc=hmax;
			} else
				hacc=h;
			if (hacc < 0.5*h2) hacc=0.5*h2;
		}
		if (hacc < hmin) hacc=hmin;
		h=hacc;
		if (h*sigmal > 1.0) {
			a=fabs((*diameter)(*t,m0,m)/sigmal+FLT_EPSILON)/2.0;
			b=2.0*fabs(sin(phil));
			betan=((a > b) ? 1.0/a : 1.0/b)/a;
			hstab=fabs(betan/sigmal);
			if (hstab < 1.0e-14*(*t)) break;
			if (h > hstab) h=hstab;
		}
		hcr=h2*h2/h1;
		if (kl > 2 && fabs(h-hcr) < FLT_EPSILON*hcr)
			h = (h < hcr) ? hcr*(1.0-FLT_EPSILON) :
								hcr*(1.0+FLT_EPSILON);
		if ((*t)+h > te) {
			last=1;
			*hstart = h;
			h=te-(*t);
		}
		h0=h1;
		h1=h2;
		h2=h;
		/* coefficient */
		b=h*sigmal;
		b1=b*cos(phil);
		bb=b*b;
		if (fabs(b) < 1.0e-3) {
			beta2=0.5-bb/24.0;
			beta3=1.0/6.0+b1/12.0;
			betha[3]=0.5+b1/3.0;
		} else if (b1 < -40.0) {
			beta2=(-2.0*b1-4.0*b1*b1/bb+1.0)/bb;
			beta3=(1.0+2.0*b1/bb)/bb;
			betha[3]=1.0/bb;
		} else {
			e=exp(b1)/bb;
			b2=b*sin(phil);
			beta2=(-2.0*b1-4.0*b1*b1/bb+1.0)/bb;
			beta3=(1.0+2.0*b1/bb)/bb;
			if (fabs(b2/b) < 1.0e-5) {
				beta2 -= e*(b1-3.0);
				beta3 += e*(b1-2.0)/b1;
				betha[3]=1.0/bb+e*(b1-1.0);
			} else {
				beta2 -= e*sin(b2-3.0*phil)/b2*b;
				beta3 += e*sin(b2-2.0*phil)/b2;
				betha[3]=1.0/bb+e*sin(b2-phil)/b2*b;
			}
		}
		beta[1]=betha[1]=1.0;
		beta[2]=beta2;
		beta[3]=beta3;
		betha[2]=1.0-bb*beta3;
		b=fabs(b);
		q = (b < 1.5) ? 4.0-2.0*b/3.0 :
					((b < 6.0) ? (30.0-2.0*b)/9.0 : 2.0);
		for (i=1; i<=3; i++) {
			hi *= h;
			if (i > 1) (*derivative)(*t,m0,m,i,c);
			/* local error construction */
			if (i == 1) inivec(m0,m,ro,0.0);
			if (i < 4) elmvec(m0,m,0,ro,c,betha[i]*hi);
			if (i == 4) {
				elmvec(m0,m,0,ro,c,-h);
				s=0.0;
				if (norm == 1)
					for (j=m0; j<=m; j++) {
						x=fabs(ro[j]);
						if (x > s) s=x;
					}
				else
					s=sqrt(vecvec(m0,m,0,ro,ro));
				*rho=s;
				ec0=ec1;
				ec1=ec2;
				ec2=(*rho)/pow(h,q);
			}
			elmvec(m0,m,0,u,c,beta[i]*hi);
		}
		t2=(*t);
		(*k)++;
		if (last) {
			last=0;
			(*t)=te;
			start=1;
		} else
			(*t) += h;
		dupvec(m0,m,0,c,u);
		(*derivative)(*t,m0,m,1,c);
		/* local error construction */
		elmvec(m0,m,0,ro,c,-h);
		s=0.0;
		if (norm == 1)
			for (j=m0; j<=m; j++) {
				x=fabs(ro[j]);
				if (x > s) s=x;
			}
		else
			s=sqrt(vecvec(m0,m,0,ro,ro));
		*rho=s;
		ec0=ec1;
		ec1=ec2;
		ec2=(*rho)/pow(h,q);
		(*out)(*t,te,m0,m,u,*k,*eta,*rho);
	} while (*t != te);
	free_real_vector(ro,m0);
}

