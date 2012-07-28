#include "../real.h"
void conjgrad(void (*matvec)(real_t [], real_t []), real_t x[],
					real_t r[], int l, int n, int (*goon)(int, real_t),
					int *iterate, real_t *norm2)
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	real_t vecvec(int, int, int, real_t [], real_t []);
	void elmvec(int, int, int, real_t [], real_t [], real_t);
	int i;
	real_t a,b,prr,rrp,*p,*ap;

	p=allocate_real_vector(l,n);
	ap=allocate_real_vector(l,n);
	*iterate=0;
	do {
		if (*iterate == 0) {
			(*matvec)(x,p);
			for (i=l; i<=n; i++) p[i] = r[i] -= p[i];
			prr=vecvec(l,n,0,r,r);
		} else {
			b=rrp/prr;
			prr=rrp;
			for (i=l; i<=n; i++) p[i]=r[i]+b*p[i];
		}
		(*matvec)(p,ap);
		a=prr/vecvec(l,n,0,p,ap);
		elmvec(l,n,0,x,p,a);
		elmvec(l,n,0,r,ap,-a);
		*norm2=rrp=vecvec(l,n,0,r,r);
		(*iterate)++;
	} while ((*goon)(*iterate,*norm2));
	free_real_vector(p,l);
	free_real_vector(ap,l);
}
