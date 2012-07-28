#include <stdio.h>
void main ()
{
	void comkwd(real_t, real_t, real_t, real_t,
					real_t *, real_t *, real_t *, real_t *);
	real_t gr,gi,kr,ki;

	comkwd(-0.1,0.3,0.11,0.02,&gr,&gi,&kr,&ki);
	printf("x**2-2(-0.1+0.3*i)*x-(0.11+0.02*i) has roots\n"
			" %6.2f+%4.2f*i\n %6.3f+%4.2f*i\n",gr,gi,kr,ki);
}

