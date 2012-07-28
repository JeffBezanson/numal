#include "../real.h"


real_t rnk1min(int n, real_t x[], real_t g[], real_t h[],
					real_t (*funct)(int, real_t[], real_t[]),
					real_t in[], real_t out[])
{
	real_t *allocate_real_vector(int, int);
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_vector(real_t *, int);
	void free_real_matrix(real_t **, int, int, int);
	real_t vecvec(int, int, int, real_t [], real_t []);
	real_t matvec(int, int, int, real_t **, real_t []);
	real_t tamvec(int, int, int, real_t **, real_t []);
	void elmvec(int, int, int, real_t [], real_t [], real_t);
	real_t symmatvec(int, int, int, real_t [], real_t []);
	void inivec(int, int, real_t [], real_t);
	void inisymd(int, int, int, real_t [], real_t);
	void mulvec(int, int, int, real_t [], real_t [], real_t);
	void dupvec(int, int, int, real_t [], real_t []);
	void eigsym1(real_t [], int, int, real_t [], real_t **, real_t []);
	void linemin(int, real_t [], real_t [], real_t, real_t *, real_t [],
					real_t (*)(int, real_t[], real_t[]), real_t, real_t *,
					real_t, real_t *, int *, int, real_t []);
	void rnk1upd(real_t [], int, real_t [], real_t);
	void davupd(real_t [], int, real_t [], real_t [], real_t, real_t);
	void fleupd(real_t [], int, real_t [], real_t [], real_t, real_t);
	int i,it,n2,cntl,cnte,evl,evlmax,ok;
	real_t f,f0,fmin,mu,dg,dg0,ghg,gs,nrmdelta,alfa,macheps,reltol,
			abstol,eps,tolg,orth,aid,*v,*delta,*gamma,*s,*p,**vec,*th,
			em[10],temp1,temp2;

	v=allocate_real_vector(1,n);
	delta=allocate_real_vector(1,n);
	gamma=allocate_real_vector(1,n);
	s=allocate_real_vector(1,n);
	p=allocate_real_vector(1,n);
	vec=allocate_real_matrix(1,n,1,n);

	macheps=in[0];
	reltol=in[1];
	abstol=in[2];
	mu=in[3];
	tolg=in[4];
	fmin=in[5];
	it=0;
	alfa=in[6];
	evlmax=in[7];
	orth=in[8];
	n2=(n*(n+1))/2;
	cntl=cnte=0;
	if (alfa > 0.0) {
		inivec(1,n2,h,0.0);
		inisymd(1,n,0,h,alfa);
	}
	f=(*funct)(n,x,g);
	evl=1;
	dg=sqrt(vecvec(1,n,0,g,g));
	for (i=1; i<=n; i++) delta[i] = -symmatvec(1,n,i,h,g);
	nrmdelta=sqrt(vecvec(1,n,0,delta,delta));
	dg0=vecvec(1,n,0,delta,g);
	ok = dg0 < 0.0;
	eps=sqrt(vecvec(1,n,0,x,x))*reltol+abstol;
	it++;
	while ((nrmdelta > eps || dg > tolg || !ok) && (evl < evlmax)) {
		if (!ok) {
			/* calculating greenstadts direction */
			th=allocate_real_vector(1,n2);
			em[0]=macheps;
			em[2]=aid=sqrt(macheps*reltol);
			em[4]=orth;
			em[6]=aid*n;
			em[8]=5.0;
			cnte++;
			dupvec(1,n2,0,th,h);
			eigsym1(th,n,n,v,vec,em);
			for (i=1; i<=n; i++) {
				aid = -tamvec(1,n,i,vec,g);
				s[i]=aid*fabs(v[i]);
				v[i]=((v[i] == 0.0) ? 0.0 : ((v[i] > 0.0) ? aid : -aid));
			}
			for (i=1; i<=n; i++) {
				delta[i]=matvec(1,n,i,vec,s);
				p[i]=matvec(1,n,i,vec,v);
			}
			dg0=vecvec(1,n,0,delta,g);
			nrmdelta=sqrt(vecvec(1,n,0,delta,delta));
			free_real_vector(th,1);
		}
		dupvec(1,n,0,s,x);
		dupvec(1,n,0,v,g);
		if (it > n)
			alfa=1.0;
		else {
			if (it != 1)
				alfa /= nrmdelta;
			else {
				alfa=2.0*(fmin-f)/dg0;
				if (alfa > 1.0) alfa=1.0;
			}
		}
		elmvec(1,n,0,x,delta,alfa);
		f0=f;
		f=(*funct)(n,x,g);
		evl++;
		dg=vecvec(1,n,0,delta,g);
		if (it == 1 || f0-f < -mu*dg0*alfa) {
			/* line minimization */
			i=evlmax-evl;
			cntl++;
			linemin(n,s,delta,nrmdelta,&alfa,g,funct,f0,&f,
						dg0,&dg,&i,0,in);
			evl += i;
			dupvec(1,n,0,x,s);
		}
		dupvec(1,n,0,gamma,g);
		elmvec(1,n,0,gamma,v,-1.0);
		if (!ok) mulvec(1,n,0,v,p,-1.0);
		dg -= dg0;
		if (alfa != 1.0) {
			mulvec(1,n,0,delta,delta,alfa);
			mulvec(1,n,0,v,v,alfa);
			nrmdelta *= alfa;
			dg *= alfa;
		}
		dupvec(1,n,0,p,gamma);
		elmvec(1,n,0,p,v,1.0);
		for (i=1; i<=n; i++) v[i]=symmatvec(1,n,i,h,gamma);
		dupvec(1,n,0,s,delta);
		elmvec(1,n,0,s,v,-1.0);
		gs=vecvec(1,n,0,gamma,s);
		ghg=vecvec(1,n,0,v,gamma);
		aid=dg/gs;
		temp1=vecvec(1,n,0,delta,p);
		temp2=orth*nrmdelta;
		if (temp1*temp1 > vecvec(1,n,0,p,p)*temp2*temp2)
			rnk1upd(h,n,s,1.0/gs);
		else
			if (aid >= 0.0)
				fleupd(h,n,delta,v,1.0/dg,(1.0+ghg/dg)/dg);
			else
				davupd(h,n,delta,v,1.0/dg,1.0/ghg);
		for (i=1; i<=n; i++) delta[i] = -symmatvec(1,n,i,h,g);
		alfa=nrmdelta;
		nrmdelta=sqrt(vecvec(1,n,0,delta,delta));
		eps=sqrt(vecvec(1,n,0,x,x))*reltol+abstol;
		dg=sqrt(vecvec(1,n,0,g,g));
		dg0=vecvec(1,n,0,delta,g);
		ok = dg0 <= 0.0;
		it++;
	}
	out[0]=nrmdelta;
	out[1]=dg;
	out[2]=evl;
	out[3]=cntl;
	out[4]=cnte;
	free_real_vector(v,1);
	free_real_vector(delta,1);
	free_real_vector(gamma,1);
	free_real_vector(s,1);
	free_real_vector(p,1);
	free_real_matrix(vec,1,n,1);
	return f;
}
