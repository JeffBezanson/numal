#include <stdio.h>
void main ()
{
	void decsolbnd(real_t [], int, int, int, real_t [], real_t []);
	real_t determbnd(real_t [], int, int, int, int);
	int i;
	real_t band[14],right[6],aux[6];

	for (i=1; i<=13; i++)
		band[i] = (((i+1)/3)*3 < i) ? 2.0 : -1.0;
	right[1]=right[5]=1.0;
	right[2]=right[3]=right[4]=0.0;
	aux[2]=1.0e-12;
	decsolbnd(band,5,1,1,aux,right);
	if (aux[3] == 5)
		printf("Delivers: %8.4f %8.4f %8.4f %8.4f %8.4f\n"
				"Determinant is  %e\n",right[1],right[2],right[3],
				right[4],right[5],determbnd(band,5,1,1,aux[1]));
}

