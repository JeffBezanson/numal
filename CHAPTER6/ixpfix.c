#include "../real.h"


void ixpfix(real_t x, real_t p, real_t q, int nmax, real_t eps,
				real_t i[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	real_t incbeta(real_t, real_t, real_t, real_t);
	void forward(real_t, real_t, real_t, real_t, real_t, int, real_t []);
	void backward(real_t, real_t, real_t, real_t, int, real_t, real_t []);
	int m,mmax;
	real_t s,p0,i0,i1,iq0,iq1,*ip;

	m=floor(p);
	s=p-m;
	p0 = (s > 0.0) ? s : s+1.0;
	mmax = (s > 0.0) ? m : m-1;
	i0=incbeta(x,p0,q,eps);
	i1=incbeta(x,p0,q+1.0,eps);
	ip=allocate_real_vector(0,mmax);
	backward(x,p0,q,i0,mmax,eps,ip);
	iq0=ip[mmax];
	backward(x,p0,q+1.0,i1,mmax,eps,ip);
	iq1=ip[mmax];
	free_real_vector(ip,0);
	forward(x,p,q,iq0,iq1,nmax,i);
}
