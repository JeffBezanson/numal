#include "../real.h"
real_t *allocate_real_vector(int, int);
void free_real_vector(real_t *, int);
real_t quanewbnd1t(real_t, int);
void quanewbnd(int, int, int, real_t [], real_t [], real_t [],
			   int (*)(int, int, int, real_t[], real_t[]),
			   real_t [], real_t []);
void jacobnbndf(int, int, int, real_t [], real_t [],
				real_t [], real_t (*)(int),
				int (*)(int, int, int, real_t[], real_t[]));
real_t quanewbnd1s(int i);

void quanewbnd1(int n, int lw, int rw, real_t x[], real_t f[],
					int (*funct)(int, int, int, real_t[], real_t[]),
					real_t in[], real_t out[])
{
	int k;
	real_t *jac;

	jac=allocate_real_vector(1,(lw+rw)*(n-1)+n);
	(*funct)(n,1,n,x,f);
	k=(lw+rw)*(n-1)+n*2-((lw-1)*lw+(rw-1)*rw)/2;
	in[4] -= k;
	quanewbnd1t(in[5], 1);
	jacobnbndf(n,lw,rw,x,f,jac,quanewbnd1s,funct);
	quanewbnd(n,lw,rw,x,f,jac,funct,in,out);
	in[4] += k;
	out[3] += k;
	free_real_vector(jac,1);
}

real_t quanewbnd1s(int i)
{
	/* this function is used internally by QUANEWBND1 */

	real_t quanewbnd1t(real_t, int);

	return (quanewbnd1t(0.0,0));
}

real_t quanewbnd1t(real_t x, int i)
{
	/* this function is used internally by QUANEWBND1 */

	static real_t y;

	y = (i ? x : y);
	return (y);
}

