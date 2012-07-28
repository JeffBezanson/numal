#include <stdio.h>
void main ()
{
	void allzerortpol(int, real_t [], real_t [], real_t [], real_t []);
	real_t b[4],c[4],zer[4],em[6];

	em[0]=em[2]=1.0e-6;  em[4]=15.0;
	b[2]=b[1]=b[0]=0.0;
	c[0]=0.0;  c[1]=0.5;  c[2]=0.25;
	allzerortpol(3,b,c,zer,em);
	printf("The three zeros:\n %13.6e\n %13.6e\n %13.6e\n\n"
			"EM[1]: %5.2f\nEM[3]: %9.3e\nEM[5]: %2.0f\n",
			zer[1],zer[2],zer[3],em[1],em[3],em[5]);
}

