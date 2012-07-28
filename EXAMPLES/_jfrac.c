#include <stdio.h>
void main ()
{
	real_t fouser(int, real_t, real_t []);
	int i;
	real_t p[11],q[11];

	for (i=1; i<=10; i++) {
		p[i]=1.0;
		q[i]=2.0;
	}
	q[0]=1.0;
	printf("JFRAC delivers:  %11.7f%11.7f%11.7f%11.7f",
			jfrac(7,p,q),jfrac(8,p,q),jfrac(9,p,q),jfrac(10,p,q));
}

