#include <stdio.h>
void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void gsselm(real_t **, int, real_t [], int [], int []);
	real_t determ(real_t **, int, int);
	int i,j,ri[5],ci[5];
	real_t d,**a,aux[8];

	a=allocate_real_matrix(1,4,1,4);
	for (i=1; i<=4; i++)
		for (j=1; j<=4; j++) a[i][j]=1.0/(i+j-1);
	aux[2]=1.0e-5;
	aux[4]=8;
	gsselm(a,4,aux,ri,ci);
	d = (aux[3] == 4) ? determ(a,4,aux[1]) : 0.0;
	printf("Determinant = %e",d);
	free_real_matrix(a,1,4,1);
}

