#include "../real.h"


void modifiedtaylor(real_t *t, real_t te, int m0, int m, real_t u[],
			real_t (*sigma)(real_t, int, int), real_t taumin,
			void (*derivative)(real_t, int, int, int, real_t []),
			int *k, real_t data[], real_t alfa, int norm,
			real_t (*aeta)(real_t, int, int),
			real_t (*reta)(real_t, int, int),
			real_t *eta, real_t *rho,
			void (*out)(real_t, real_t, int, int, real_t [],
							int, real_t, real_t))
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	real_t vecvec(int, int, int, real_t [], real_t []);
	int i,n,p,q,start,step1,last,j;
	real_t ec0,ec1,ec2,tau0,tau1,tau2,taus,t2,t0,tau,taui,tauec,
			ecl,betan,gamma,*c,*beta,*betha,ifac,tauacc,taustab,
			aa,bb,cc,ec,s,x,b;

	n=data[-2];
	beta=allocate_real_vector(1,n);
	betha=allocate_real_vector(1,n);
	c=allocate_real_vector(m0,m);

	i=0;
	start=(*k == 0);
	t0=(*t);
	/* coefficient */
	ifac=1.0;
	gamma=0.5;
	p=data[-1];
	betan=data[0];
	q = (p < n) ? p+1 : n;
	for (j=1; j<=n; j++) {
		beta[j]=data[j];
		ifac /= j;
		betha[j]=ifac-beta[j];
	}
	if (p == n) betha[n]=ifac;
	last=0;
	do {
		/* step size */
		s=0.0;
		if (norm == 1)
			for (j=m0; j<=m; j++) {
				x=fabs(u[j]);
				if (x > s) s=x;
			}
		else
			s=sqrt(vecvec(m0,m,0,u,u));
		/* local error bound */
		*eta = (*aeta)(*t,m0,m)+(*reta)(*t,m0,m)*s;
		if (*eta > 0.0) {
			if (start) {
				if (*k == 0) {
					for (j=m0; j<=m; j++) c[j]=u[j];
					i=1;
					(*derivative)(*t,m0,m,i,c);
					s=0.0;
					if (norm == 1)
						for (j=m0; j<=m; j++) {
							x=fabs(c[j]);
							if (x > s) s=x;
						}
					else
						s=sqrt(vecvec(m0,m,0,c,c));
					tauacc=(*eta)/s;
					step1=1;
				} else if (step1) {
					tauacc=pow((*eta)/(*rho),1.0/q)*tau2;
					if (tauacc > 10.0*tau2)
						tauacc=10.0*tau2;
					else
						step1=0;
				} else {
					bb=(ec2-ec1)/tau1;
					cc=ec2-bb*t2;
					ec=bb*(*t)+cc;
					tauacc = (ec < 0.0) ? tau2 : pow((*eta)/ec,1.0/q);
					start=0;
				}
			} else {
				aa=((ec0-ec1)/tau0+(ec2-ec1)/tau1)/(tau1+tau0);
				bb=(ec2-ec1)/tau1-aa*(2.0*t2-tau1);
				cc=ec2-t2*(bb+aa*t2);
				ec=cc+(*t)*(bb+(*t)*aa);
				tauacc = (ec < 0.0) ? taus : pow((*eta)/ec,1.0/q);
				if (tauacc > alfa*taus) tauacc=alfa*taus;
				if (tauacc < gamma*taus) tauacc=gamma*taus;
			}
		} else
			tauacc=te-(*t);
		if (tauacc < taumin) tauacc=taumin;
		taustab=betan/(*sigma)(*t,m0,m);
		if (taustab < 1.0e-12*((*t)-t0)) {
			(*out)(*t,te,m0,m,u,*k,*eta,*rho);
			break;
		}
		tau = (tauacc > taustab) ? taustab : tauacc;
		taus=tau;
		if (tau >= te-(*t)) {
			tau=te-(*t);
			last=1;
		}
		tau0=tau1;
		tau1=tau2;
		tau2=tau;
		(*k)++;
		i=0;
		/* difference scheme */
		for (j=m0; j<=m; j++) c[j]=u[j];
		taui=1.0;
		do {
			i++;
			(*derivative)(*t,m0,m,i,c);
			taui *= tau;
			b=beta[i]*taui;
			if (*eta > 0.0 && i >= p) {
				/* local error construction */
				if (i == p) {
					ecl=0.0;
					tauec=1.0;
				}
				if (i > p+1) tauec *= tau;
				s=0.0;
				if (norm == 1)
					for (j=m0; j<=m; j++) {
						x=fabs(c[j]);
						if (x > s) s=x;
					}
				else
					s=sqrt(vecvec(m0,m,0,c,c));
				ecl += fabs(betha[i])*tauec*s;
				if (i == n) {
					ec0=ec1;
					ec1=ec2;
					ec2=ecl;
					*rho = ecl*pow(tau,q);
				}
			}
			for (j=m0; j<=m; j++) u[j] += b*c[j];
		} while (i < n);
		t2=(*t);
		if (last) {
			last=0;
			(*t)=te;
		} else
			(*t) += tau;
		(*out)(*t,te,m0,m,u,*k,*eta,*rho);
	} while (*t != te);
	free_real_vector(beta,1);
	free_real_vector(betha,1);
	free_real_vector(c,m0);
}

