#include "../real.h"
void nonexpspherbessk(real_t x, int n, real_t k[])
{
	int i;
	real_t ki,ki1,ki2;
	x=1.0/x;
	k[0]=ki2=x*1.5707963267949;
	if (n != 0) {
		k[1]=ki1=ki2*(1.0+x);
		for (i=2; i<=n; i++) {
			k[i]=ki=ki2+(i+i-1)*x*ki1;
			ki2=ki1;
			ki1=ki;
		}
	}
}
