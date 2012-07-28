#include <stdio.h>
void main ()
{
	void selzerortpol(int, int, int, real_t [], real_t [],
							real_t [], real_t []);
	int i;
	real_t b[5],c[5],zer[4],em[6];

	em[0]=em[2]=1.0e-6;
	for (i=0; i<=3; i++) {
		b[i]=0.0;
		c[i]=i*i/(4.0*i*i-1.0);
	}
	selzerortpol(4,3,3,b,c,zer,em);
	printf("The third zero:\n %13.6e\n\nEM[1]: %5.2f\nEM[5]: %3.0f\n",
			zer[3],em[1],em[5]);
}

