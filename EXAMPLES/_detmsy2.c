#include <stdio.h>
void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void decsym2(real_t **, int, real_t, int [], int [], real_t []);
	real_t detmsym2(real_t [], int, int []);
	int i,j,aux[6],p[6];
	real_t tol,determinant,**a,detaux[6];

	a=allocate_real_matrix(1,5,1,5);
	a[1][1]=a[1][2] = -3.0; a[1][3] = -18.0;
	a[1][4] = -30.0;  a[1][5]=18.0;
	a[2][2] = -1.0;	a[2][3] = -4.0;    a[2][4] = -48.0;  a[2][5]=8.0;
	a[3][3] = -6.0;	a[3][4] = -274.0;  a[3][5]=6.0;
	a[4][4]=119.0;	a[4][5]=19.0;    a[5][5]=216.0;
	for (i=1; i<=5; i++)
		for (j=i+1; j<=5; j++) a[j][i]=a[i][j];
	tol=1.0e-6;
	decsym2(a,5,tol,aux,p,detaux);
	if (aux[2] == 1)
		printf("\nThe matrix is symmetric.");
	else
		printf("The matrix is asymmetric, results are meaningless.");
	determinant=determsym2(detaux,5,aux);
	printf("\nThe determinant of the matrix :%8.2f\n",determinant);
	free_real_matrix(a,1,5,1);
}

