#include <stdio.h>
void main ()
{
	void besspq0(real_t, real_t *, real_t *);
	void besspq1(real_t, real_t *, real_t *);
	int i;
	real_t p,q,r,s;
	real_t x[4]={1.0, 3.0, 5.0, 10.0};

	printf("BESSPQ0 and BESSPQ1 deliver:\n");
	for (i=0; i<=3; i++) {
		besspq0(x[i],&p,&q);
		besspq1(x[i],&r,&s);
		printf("   %7.2e",fabs(p*r+q*s-1.0));
	}
}

