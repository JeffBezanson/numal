#include "../real.h"
void chlsol1(real_t a[], int n, real_t b[])
{
	real_t vecvec(int, int, int, real_t [], real_t []);
	real_t seqvec(int, int, int, int, real_t [], real_t []);
	int i,ii;

	ii=0;
	for (i=1; i<=n; i++) {
		ii += i;
		b[i]=(b[i]-vecvec(1,i-1,ii-i,b,a))/a[ii];
	}
	for (i=n; i>=1; i--) {
		b[i]=(b[i]-seqvec(i+1,n,ii+i,0,a,b))/a[ii];
		ii -= i;
	}
}
