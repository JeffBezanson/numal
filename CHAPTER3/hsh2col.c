#include "../real.h"


void hsh2col(int la, int lb, int u, int i, real_t a1, real_t a2,
					real_t **a, real_t **b)
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void hshvecmat(int, int, int, int, real_t, real_t [], real_t **);
	real_t *v,d1,d2,s1,s2,r,d,c;

	if (a2 != 0.0) {
		v=allocate_real_vector(i,i+1);
		d1=fabs(a1);
		d2=fabs(a2);
		s1 = (a1 >= 0.0) ? 1.0 : -1.0;
		s2 = (a2 >= 0.0) ? 1.0 : -1.0;
		if (d2 <= d1) {
			r=d2/d1;
			d=sqrt(1.0+r*r);
			c = -1.0-1.0/d;
			v[i+1]=s1*s2*r/(1.0+d);
		} else {
			r=d1/d2;
			d=sqrt(1.0+r*r);
			c = -1.0-r/d;
			v[i+1]=s1*s2/(r+d);
		}
		v[i]=1.0;
		hshvecmat(i,i+1,la,u,c,v,a);
		hshvecmat(i,i+1,lb,u,c,v,b);
		free_real_vector(v,i);
	}
}
