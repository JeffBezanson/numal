#include <stdio.h>
void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void gssinv(real_t **, int, real_t []);
	int i;
	real_t **a,aux[10];

	a=allocate_real_matrix(1,4,1,4);
	a[1][1]=4.0;  a[1][2]=2.0;  a[1][3]=4.0;  a[1][4]=1.0;
	a[2][1]=30.0; a[2][2]=20.0; a[2][3]=45.0; a[2][4]=12.0;
	a[3][1]=20.0; a[3][2]=15.0; a[3][3]=36.0; a[3][4]=10.0;
	a[4][1]=35.0; a[4][2]=28.0; a[4][3]=70.0; a[4][4]=20.0;
	aux[2]=1.0e-5;
	aux[4]=8;
	gssinv(a,4,aux);
	printf("Calculated inverse:\n");
	for (i=1; i<=4; i++)
		printf(" %4.0f%4.0f%4.0f%4.0f\n",a[i][1],a[i][2],a[i][3],a[i][4]);
	printf("\nAUX elements:\n%e  %e  %e  %e  %e\n",
			aux[1],aux[3],aux[5],aux[7],aux[9]);
	free_real_matrix(a,1,4,1);
}

