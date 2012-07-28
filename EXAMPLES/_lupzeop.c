#include <stdio.h>
void main ()
{
	void lupzerortpol(int, int, real_t [], real_t [], real_t [], real_t []);
	int i;
	real_t b[4],c[4],zer[3],em[7];

	em[0]=em[2]=1.0e-6;  em[4]=45.0;  em[6]=1.0;
	for (i=0; i<=2; i++) {
		b[i]=2*i+1;
		c[i]=i*i;
	}
	lupzerortpol(3,2,b,c,zer,em);
	printf("The two lower zeros:\n %13.6e\n %13.6e\n\n"
			"EM[1]: %5.2f\nEM[3]: %9.3e\nEM[5]: %3.0f\n",
			zer[1],zer[2],em[1],em[3],em[5]);
	em[6]=0.0;
	for (i=0; i<=2; i++) {
		b[i] = -2*i-1;
		c[i]=i*i;
	}
	lupzerortpol(3,2,b,c,zer,em);
	printf("\nThe two upper zeros:\n %13.6e\n %13.6e\n\n"
			"EM[1]: %5.2f\nEM[3]: %9.3e\nEM[5]: %3.0f\n",
			-zer[1],-zer[2],em[1],em[3],em[5]);
}

