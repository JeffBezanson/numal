#include <math.h>
#include <stdio.h>
void main ()
{
	void decsoltri(real_t [], real_t [], real_t [], int,
						real_t [], real_t []);
	real_t vecvec(int, int, int, real_t [], real_t []);
	int i;
	real_t d[31],sub[31],super[31],b[31],aux[6];

	for (i=1; i<=30; i++) {
		sub[i]=i*2;
		super[i]=i;
		d[i]=i+10;
		b[i]=0.0;
	}
	b[1]=1.0;  b[2]=12.0;  b[3]=4.0;
	aux[2]=1.0e-6;
	decsoltri(sub,d,super,30,aux,b);
	b[2]--;
	printf("AUX[3] and AUX[5]:  %e   %e\n"
			"Error in the solution:  %e\n",
			aux[3],aux[5],sqrt(vecvec(1,30,0,b,b)));
}

