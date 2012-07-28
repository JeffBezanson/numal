#include <stdio.h>
void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void eigvalsym2(real_t **, int, int, real_t [], real_t []);
	int i,j;
	real_t **a,val[3],em[4];

	a=allocate_real_matrix(1,4,1,4);
	em[0]=1.0e-6;
	em[2]=1.0e-5;
	for (i=1; i<=4; i++)
		for (j=i; j<=4; j++) a[i][j]=1.0/(i+j-1);
	eigvalsym2(a,4,2,val,em);
	printf("The eigenvalues:  %12.5e   %12.5e\n"
			"EM[1] = %e\nEM[3] =%3.0f\n",val[1],val[2],em[1],em[3]);
	free_real_matrix(a,1,4,1);
}

