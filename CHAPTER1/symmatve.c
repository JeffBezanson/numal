#include "../real.h"
real_t symmatvec(int l, int u, int i, real_t a[], real_t b[])
{
	int k, m;
	real_t vecvec(int, int, int, real_t [], real_t []);
	real_t seqvec(int, int, int, int, real_t [], real_t []);

	m=(l>i) ? l : i;
	k=(m*(m-1))/2;
	return (vecvec(l, (i<=u) ? i-1 : u, k,b,a) + seqvec(m,u,k+i,0,a,b));
}
