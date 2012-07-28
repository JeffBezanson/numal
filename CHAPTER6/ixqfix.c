#include "../real.h"


void ixqfix(real_t x, real_t p, real_t q, int nmax, real_t eps,
				real_t i[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	real_t incbeta(real_t, real_t, real_t, real_t);
	void forward(real_t, real_t, real_t, real_t, real_t, int, real_t []);
	void backward(real_t, real_t, real_t, real_t, int, real_t, real_t []);
	int m,mmax;
	real_t s,iq0,iq1,q0,*iq;

	m=floor(q);
	s=q-m;
	q0 = (s > 0.0) ? s : s+1.0;
	mmax = (s > 0.0) ? m : m-1;
	iq0=incbeta(x,p,q0,eps);
	if (mmax > 0) iq1=incbeta(x,p,q0+1.0,eps);
	iq=allocate_real_vector(0,mmax);
	forward(x,p,q0,iq0,iq1,mmax,iq);
	backward(x,p,q,iq[mmax],nmax,eps,i);
	free_real_vector(iq,0);
}
