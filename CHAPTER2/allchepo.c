#include "../real.h"
void allchepol(int n, real_t x, real_t t[])
{
	int i;
	real_t t1,t2,h,x2;

	if (n == 0) {
		t[0]=1.0;
		return;
	}
	if (n == 1) {
		t[0]=1.0;
		t[1]=x;
		return;
	}
	t[0]=t1=1.0;
	t[1]=t2=x;
	x2=x+x;
	for (i=2; i<=n; i++) {
		t[i]=h=x2*t2-t1;
		t1=t2;
		t2=h;;
	}
}
