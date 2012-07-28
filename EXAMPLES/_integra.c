#include <stdio.h>

real_t fx(real_t x)
{
	return 10.0/(x*x);
}

void main ()
{
	real_t integral(real_t, real_t, real_t (*)(real_t), real_t [],
					int, int);
	int ua,ub,i;
	real_t e[7],a;
	static real_t b[4]={2.0, 4.0, 20.0, 100.0};

	printf("INTEGRAL delivers:\n");
	ua=1;
	e[1]=e[2]=1.0e-6;
	for (i=0; i<=3; i++) {
		ub=(b[i] < 50.0);
		a=integral(-1.0,-b[i],fx,e,ua,ub);
		printf(" %e  %3.0f   %e   %3.0f  %3.0f\n",
				a,e[3],e[4],e[5],e[6]);
		ua=0;
	}
}

