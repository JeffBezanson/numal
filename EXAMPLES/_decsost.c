#include <math.h>
#include <stdio.h>
void main ()
{
	void decsolsymtri(real_t [], real_t [], int, real_t [], real_t []);
	real_t vecvec(int, int, int, real_t [], real_t []);
	int i;
	real_t d[101],co[101],b[101],aux[6];

	for (i=1; i<=100; i++) {
		d[i]=i;
		co[i]=i*2;
		b[i]=0.0;
	}
	b[1]=b[2]=2.0;  b[3]=4.0;
	aux[2]=1.0e-6;
	decsolsymtri(d,co,100,aux,b);
	b[2]--;
	printf("AUX[3] and AUX[5]:  %e   %e\n"
			"Error in the solution:  %e\n",
			aux[3],aux[5],sqrt(vecvec(1,100,0,b,b)));
}

