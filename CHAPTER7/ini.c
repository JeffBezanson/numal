#include "../real.h"


void ini(int n, int m, int s[])
{
	int i,j,k,l;
	real_t pin2,temp;

	pin2=atan(1.0)*2.0/n;
	k=1;
	l=n-1;
	j=s[0]=0;
	s[n]=m;
	while (k < l) {
		temp=sin(k*pin2);
		i=temp*temp*m;
		j = s[k] = ((i <= j) ? j+1 : i);
		s[l]=m-j;
		l--;
		k++;
	}
	if (l*2 == n) s[l]=m/2;
}
