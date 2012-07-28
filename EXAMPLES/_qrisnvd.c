#include <stdio.h>
void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	int qrisngvaldec(real_t **, int, int, real_t [], real_t **, real_t []);
	int i,j;
	real_t **a,**v,val[5],em[8];

	a=allocate_real_matrix(1,6,1,5);
	v=allocate_real_matrix(1,5,1,5);
	for (i=1; i<=6; i++)
		for (j=1; j<=5; j++) a[i][j]=1.0/(i+j-1);
	em[0]=1.0e-6;  em[2]=1.0e-5;  em[4]=25.0;  em[6]=1.0e-5;
	i=qrisngvaldec(a,6,5,val,v,em);
	printf("Number of singular values not found : %2d\n"
		"Infinity norm : %e\nMax neglected subdiagonal element : %e\n"
		"Number of iterations : %3.0f\nNumerical rank : %3.0f\n"
		"\nSingular values :\n",
		i,em[1],em[3],em[5],em[7]);
	for (i=1; i<=5; i++)
		printf("  %12.6e\n",val[i]);
	printf("\nMatrix U, first 3 columns :\n");
	for (i=1; i<=6; i++)
		printf("  %12.6e  %12.6e  %12.6e\n",a[i][1],a[i][2],a[i][3]);
	printf("\n          Last 2 columns :\n");
	for (i=1; i<=6; i++)
		printf("          %12.6e  %12.6e\n",a[i][4],a[i][5]);
	free_real_matrix(a,1,6,1);
	free_real_matrix(v,1,5,1);
}

