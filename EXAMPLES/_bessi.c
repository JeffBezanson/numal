#include <stdio.h>
void main ()
{
	void bessi(real_t, int, real_t []);
	void bessk(real_t, int, real_t []);
	int j,n;
	real_t x,i[6],k[6];

	printf("BESSI and BESSK deliver:");
	for (j=1; j<=20; j++) {
		x=j;
		printf("\n %3.0f",x);
		bessi(x,5,i);
		bessk(x,5,k);
		for (n=1; n<=5; n++)
			printf(" %+9.2e",x*(i[n]*k[n-1]+i[n-1]*k[n])-1.0);
	}
}

