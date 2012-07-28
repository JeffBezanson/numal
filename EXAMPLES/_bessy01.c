#include <stdio.h>
void main ()
{
	void bessy01(real_t, real_t *, real_t *);
	real_t x,y0,y1;

	x=1.0;
	bessy01(x,&y0,&y1);
	printf("BESSY01 delivers:\n %3.1f   %e   %e\n",x,y0,y1);
}

