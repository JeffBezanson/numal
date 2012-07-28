#include <stdio.h>
void main ()
{
	void fresnel(real_t, real_t *, real_t *);
	void fg(real_t, real_t *, real_t *);
	real_t c,s,f,g;

	fresnel(1.0,&c,&s);
	fg(1.0,&f,&g);
	printf("FRESNEL and FG deliver:\n\n C(1) = %e   S(1) = %e\n"
			" F(1) = %e   G(1) = %e\n",c,s,f,g);
}

