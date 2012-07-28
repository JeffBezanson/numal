#include "../real.h"


void bessy(real_t x, int n, real_t y[])
{
	void bessy01(real_t, real_t *, real_t *);
	int i;
	real_t y0,y1,y2;

	bessy01(x,&y0,&y1);
	y[0]=y0;
	if (n > 0) y[1]=y1;
	x=2.0/x;
	for (i=2; i<=n; i++) {
		y[i]=y2=(i-1)*x*y1-y0;
		y0=y1;
		y1=y2;
	}
}
