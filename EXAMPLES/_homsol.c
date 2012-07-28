#define real_t double

#include <stdio.h>
void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	int homsol(real_t **, int, int, real_t **, real_t []);
	int i,j;
	real_t **a,**v,em[8];

	a=allocate_real_matrix(1,8,1,5);
	v=allocate_real_matrix(1,5,1,5);
	a[1][1]=22;	a[1][2]=a[2][3]=10.0; a[1][3]=a[7][1]=a[8][5]=2.0;
	a[1][4]=a[3][5]=3.0; a[1][5]=a[2][2]=7.0; a[2][1]=14.0; a[2][5]=8.0;
	a[2][4]=a[8][3]=0.0; a[3][1]=a[3][3]=a[6][5] = -1.0; a[3][2]=13.0;
	a[3][4] = -11.0; a[4][1] = -3.0;
	a[4][2]=a[4][4]=a[5][4]=a[8][4] = -2.0;
	a[4][3]=13.0; a[4][5]=a[5][5]=a[8][1]=4.0; a[5][1]=a[6][1]=9.0;
	a[5][2]=8.0; a[5][3]=a[6][2]=a[7][5]=1.0; a[6][3] = -7.0;
	a[6][4]=a[7][4]=a[8][2]=5.0; a[7][2] = -6.0; a[7][3]=6.0;
	em[0]=1.0e-14; em[2]=1.0e-12; em[4]=80.0; em[6]=1.0e-10;
	i=homsol(a,8,5,v,em);
	printf("Number of singular values not found : %2d\n"
		"Norm : %e\nMaximal neglected subdiagonal element : %e\n"
		"Number of iterations : %3.0f\nRank : %3.0f\n",
		i,em[1],em[3],em[5],em[7]);
	for (j=em[7]+1; j<=5; j++) {
		printf("\nColumn number : %2d\n",j);
		for (i=1; i<=5; i++)
			printf(" %13.6e %13.6e\n",a[i][j],v[i][j]);
			printf(" %13.6e\n %13.6e\n %13.6e\n",
					a[6][j],a[7][j],a[8][j]);
	}
	free_real_matrix(a,1,8,1);
	free_real_matrix(v,1,5,1);
}

