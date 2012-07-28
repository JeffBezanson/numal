#include "../real.h"


real_t flemin(int n, real_t x[], real_t g[], real_t h[],
					real_t (*funct)(int, real_t[], real_t[]),
					real_t in[], real_t out[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	real_t vecvec(int, int, int, real_t [], real_t []);
	void elmvec(int, int, int, real_t [], real_t [], real_t);
	real_t symmatvec(int, int, int, real_t [], real_t []);
	void inivec(int, int, real_t [], real_t);
	void inisymd(int, int, int, real_t [], real_t);
	void mulvec(int, int, int, real_t [], real_t [], real_t);
	void dupvec(int, int, int, real_t [], real_t []);
	void linemin(int, real_t [], real_t [], real_t, real_t *, real_t [],
					real_t (*)(int, real_t[], real_t[]), real_t, real_t *,
					real_t, real_t *, int *, int, real_t []);
	void davupd(real_t [], int, real_t [], real_t [], real_t, real_t);
	void fleupd(real_t [], int, real_t [], real_t [], real_t, real_t);
	int i,it,cntl,evl,evlmax;
	real_t f,f0,fmin,mu,dg,dg0,nrmdelta,alfa,reltol,abstol,eps,tolg,
			aid,*v,*delta,*s;

	v=allocate_real_vector(1,n);
	delta=allocate_real_vector(1,n);
	s=allocate_real_vector(1,n);

	reltol=in[1];
	abstol=in[2];
	mu=in[3];
	tolg=in[4];
	fmin=in[5];
	alfa=in[6];
	evlmax=in[7];
	out[4]=0.0;
	it=0;
	f=(*funct)(n,x,g);
	evl=1;
	cntl=0;
	if (alfa > 0.0) {
		inivec(1,(n*(n+1))/2,h,0.0);
		inisymd(1,n,0,h,alfa);
	}
	for (i=1; i<=n; i++) delta[i] = -symmatvec(1,n,i,h,g);
	dg=sqrt(vecvec(1,n,0,g,g));
	nrmdelta=sqrt(vecvec(1,n,0,delta,delta));
	eps=sqrt(vecvec(1,n,0,x,x))*reltol+abstol;
	dg0=vecvec(1,n,0,delta,g);
	it++;
	while ((nrmdelta > eps || dg > tolg) && (evl < evlmax)) {
		dupvec(1,n,0,s,x);
		dupvec(1,n,0,v,g);
		if (it >= n)
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
		if (alfa != 1.0) mulvec(1,n,0,delta,delta,alfa);
		mulvec(1,n,0,v,v,-1.0);
		elmvec(1,n,0,v,g,1.0);
		for (i=1; i<=n; i++) s[i]=symmatvec(1,n,i,h,v);
		aid=vecvec(1,n,0,v,s);
		dg=(dg-dg0)*alfa;
		if (dg > 0.0)
			if (dg >= aid)
				fleupd(h,n,delta,s,1.0/dg,(1.0+aid/dg)/dg);
			else
				davupd(h,n,delta,s,1.0/dg,1.0/aid);
		for (i=1; i<=n; i++) delta[i] = -symmatvec(1,n,i,h,g);
		alfa *= nrmdelta;
		nrmdelta=sqrt(vecvec(1,n,0,delta,delta));
		eps=sqrt(vecvec(1,n,0,x,x))*reltol+abstol;
		dg=sqrt(vecvec(1,n,0,g,g));
		dg0=vecvec(1,n,0,delta,g);
		if (dg0 > 0.0) {
			out[4] = -1.0;
			break;
		}
		it++;
	}
	out[0]=nrmdelta;
	out[1]=dg;
	out[2]=evl;
	out[3]=cntl;
	free_real_vector(v,1);
	free_real_vector(delta,1);
	free_real_vector(s,1);
	return f;
}
