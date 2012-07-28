#include "../real.h"
void lngintpower(int u[], int exponent, int result[])
{
	int *allocate_integer_vector(int, int);
	void free_integer_vector(int *, int);
	void lngintmult(int [], int [], int []);
	int max,i,n,exp,*y,*z,*h;

	exp=exponent;
	max=u[0]*exp;
	y=allocate_integer_vector(0,max);
	z=allocate_integer_vector(0,max);
	h=allocate_integer_vector(0,max);

	y[0]=y[1]=1;
	for (i=u[0]; i>=0; i--) z[i]=u[i];
	while (1) {
		n=exp/2;
		if (n+n != exp) {
			lngintmult(y,z,h);
			for (i=h[0]; i>=0; i--) y[i]=h[i];
			if (n == 0) {
				for (i=y[0]; i>=0; i--) result[i]=y[i];
				free_integer_vector(y,0);
				free_integer_vector(z,0);
				free_integer_vector(h,0);
				return;
			}
		}
		lngintmult(z,z,h);
		for (i=h[0]; i>=0; i--) z[i]=h[i];
		exp=n;
	}
}
