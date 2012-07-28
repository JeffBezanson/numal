#include "../real.h"


void spherbessy(real_t x, int n, real_t y[])
{
	if (n == 0)
		y[0] = -cos(x)/x;
	else {
		int i;
		real_t yi,yi1,yi2;
		yi2 = y[0] = -cos(x)/x;
		yi1=y[1]=(yi2-sin(x))/x;
		for (i=2; i<=n; i++) {
			y[i] = yi = -yi2+(i+i-1)*yi1/x;
			yi2=yi1;
			yi1=yi;
		}
	}
}
