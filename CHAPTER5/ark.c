#include "../real.h"


void ark(real_t *t, real_t *te, int *m0, int *m, real_t u[],
			void (*derivative)(int *, int *, real_t *, real_t[]),
			real_t data[],
			void (*out)(int *, int *, real_t *, real_t *, real_t [],
							real_t []))
{
	real_t *allocate_real_vector(int, int);
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_vector(real_t *, int);
	void free_real_matrix(real_t **, int, int, int);
	void inivec(int, int, real_t [], real_t);
	void mulvec(int, int, int, real_t [], real_t [], real_t);
	void dupvec(int, int, int, real_t [], real_t []);
	real_t vecvec(int, int, int, real_t [], real_t []);
	void elmvec(int, int, int, real_t [], real_t [], real_t);
	void decsol(real_t **, int, real_t [], real_t []);
	real_t arkmui(int, int, int, real_t []);
	real_t arklabda(int, int, int, int, real_t []);
	static real_t th1[8] = {1.0, 0.5, 1.0/6.0, 1.0/3.0, 1.0/24.0,
		1.0/12.0, 0.125, 0.25};
	static real_t ec0,ec1,ec2,tau0,tau1,tau2,taus,t2;
	int p,n,q,start,step1,last,i,j,k,l,n1,m00;
	real_t thetanm1,tau,betan,qinv,eta,*mu,*lambda,*thetha,*ro,*r,
			**alfa,th[9],aux[4],s,ss,theta0,tauacc,taustab,
			aa,bb,cc,ec,mt,lt;

	n=data[1];
	m00=(*m0);
	mu=allocate_real_vector(1,n);
	lambda=allocate_real_vector(1,n);
	thetha=allocate_real_vector(0,n);
	ro=allocate_real_vector(m00,*m);
	r=allocate_real_vector(m00,*m);
	alfa=allocate_real_matrix(1,8,1,n+1);

	p=data[2];
	ec1=ec2=0.0;
	betan=data[3];
	thetanm1 = (p == 3) ? 0.75 : 1.0;
	theta0=1.0-thetanm1;
	s=1.0;
	for (j=n-1; j>=1; j--) {
		s = -s*theta0+data[n+10-j];
		mu[j]=data[n+11-j]/s;
		lambda[j]=mu[j]-theta0;
	}
	for (i=1; i<=8; i++)
		for (j=0; j<=n; j++)
			if (i == 1) alfa[i][j+1]=1.0;
			else if (j == 0) alfa[i][j+1]=0.0;
			else if (i == 2 || i == 4 || i == 8)
					alfa[i][j+1]=pow(arkmui(j,n,p,lambda),(i+2)/3);
			else if ((i == 3 || i == 6) && j > 1) {
				s=0.0;
				for (l=1; l<=j-1; l++)
					s += arklabda(j,l,n,p,lambda)*
								pow(arkmui(l,n,p,lambda),i/3);
				alfa[i][j+1]=s;
			}
			else if (i == 5 && j > 2) {
				s=0.0;
				for (l=2; l<=j-1; l++) {
					ss=0.0;
					for (k=1; k<=l-1; k++)
						ss += arklabda(l,k,n,p,lambda)*
									arkmui(k,n,p,lambda);
					s += arklabda(j,l,n,p,lambda)*ss;
				}
				alfa[i][j+1]=s;
			}
			else if (i == 7 && j > 1) {
				s=0.0;
				for (l=1; l<=j-1; l++)
					s += arklabda(j,l,n,p,lambda)*arkmui(l,n,p,lambda);
				alfa[i][j+1]=s*arkmui(j,n,p,lambda);
			}
			else alfa[i][j+1]=0.0;
	n1 = ((n < 4) ? n+1 : ((n < 7) ? 4 : 8));
	for (i=1; i<=8; i++) th[i]=th1[i-1];
	if (p == 3 && n < 7) th[1]=th[2]=0.0;
	aux[2]=FLT_EPSILON;
	decsol(alfa,n1,aux,th);
	inivec(0,n,thetha,0.0);
	dupvec(0,n1-1,1,thetha,th);
	if (!(p == 3 && n < 7)) {
		thetha[0] -= theta0;
		thetha[n-1] -= thetanm1;
		q=p+1;
	} else
		q=3;
	qinv=1.0/q;
	start=(data[8] == 0.0);
	data[10]=0.0;
	last=0;
	dupvec(*m0,*m,0,r,u);
	(*derivative)(m0,m,t,r);
	do {
		/* stepsize */
		eta=sqrt(vecvec(*m0,*m,0,u,u))*data[7]+data[6];
		if (eta > 0.0) {
			if (start) {
				if (data[8] == 0) {
					tauacc=data[5];
					step1=1;
				} else
					if (step1) {
						tauacc=pow(eta/ec2,qinv);
						if (tauacc > 10.0*tau2)
							tauacc=10.0*tau2;
						else
							step1=0;
					} else {
						bb=(ec2-ec1)/tau1;
						cc = -bb*t2+ec2;
						ec=bb*(*t)+cc;
						tauacc = (ec < 0.0) ? tau2 : pow(eta/ec,qinv);
						start=0;
					}
			} else {
				aa=((ec0-ec1)/tau0+(ec2-ec1)/tau1)/(tau1+tau0);
				bb=(ec2-ec1)/tau1-(2.0*t2-tau1)*aa;
				cc = -(aa*t2+bb)*t2+ec2;
				ec=(aa*(*t)+bb)*(*t)+cc;
				tauacc = ((ec < 0.0) ? taus : pow(eta/ec,qinv));
				if (tauacc > 2.0*taus) tauacc=2.0*taus;
				if (tauacc < taus/2.0) tauacc=taus/2.0;
			}
		} else
			tauacc=data[5];
		if (tauacc < data[5]) tauacc=data[5];
		taustab=betan/data[4];
		if (taustab < data[5]) {
			data[10]=1.0;
			break;
		}
		tau = ((tauacc > taustab) ? taustab : tauacc);
		taus=tau;
		if (tau >= (*te)-(*t)) {
			tau=(*te)-(*t);
			last=1;
		}
		tau0=tau1;
		tau1=tau2;
		tau2=tau;
		/* difference scheme */
		mulvec(*m0,*m,0,ro,r,thetha[0]);
		if (p == 3) elmvec(*m0,*m,0,u,r,0.25*tau);
		for (i=1; i<=n-1; i++) {
			mt=mu[i]*tau;
			lt=lambda[i]*tau;
			for (j=(*m0); j<=(*m); j++) r[j]=lt*r[j]+u[j];
			s=(*t)+mt;
			(*derivative)(m0,m,&s,r);
			if (thetha[i] != 0.0) elmvec(*m0,*m,0,ro,r,thetha[i]);
			if (i == n) {
				data[9]=sqrt(vecvec(*m0,*m,0,ro,ro))*tau;
				ec0=ec1;
				ec1=ec2;
				ec2=data[9]/pow(tau,q);
			}
		}
		elmvec(*m0,*m,0,u,r,thetanm1*tau);
		dupvec(*m0,*m,0,r,u);
		s=(*t)+tau;
		(*derivative)(m0,m,&s,r);
		if (thetha[n] != 0.0) elmvec(*m0,*m,0,ro,r,thetha[n]);
		data[9]=sqrt(vecvec(*m0,*m,0,ro,ro))*tau;
		ec0=ec1;
		ec1=ec2;
		ec2=data[9]/pow(tau,q);
		t2=(*t);
		if (last) {
			last=0;
			(*t)=(*te);
		} else
			(*t) += tau;
		data[8] += 1.0;
		(*out)(m0,m,t,te,u,data);
	} while ((*t) != (*te));
	free_real_vector(mu,1);
	free_real_vector(lambda,1);
	free_real_vector(thetha,0);
	free_real_vector(ro,m00);
	free_real_vector(r,m00);
	free_real_matrix(alfa,1,8,1);
}

real_t arkmui(int i, int n, int p, real_t lambda[])
{
	/* this function is internally used by ARK */

	return ((i==n) ? 1.0 : ((i<1 || i>n) ? 0.0 :
			((p<3) ? lambda[i] : ((p==3) ? lambda[i]+0.25 : 0.0))));
}

real_t arklabda(int i, int j, int n, int p, real_t lambda[])
{
	/* this function is internally used by ARK */

	real_t arkmui(int, int, int, real_t []);

	return ((p<3) ? ((j==i-1) ? arkmui(i,n,p,lambda) : 0.0) :
				((p==3) ? ((i==n) ? ((j==0) ? 0.25 :
				((j==n-1) ? 0.75 : 0.0)) :
				((j==0) ? ((i==1) ? arkmui(1,n,p,lambda) : 0.25) :
				((j==i-1) ? lambda[i] : 0.0))) : 0.0));
}
