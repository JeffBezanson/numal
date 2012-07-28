#include <stdio.h>
void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void decsolsym2(real_t **, int, real_t [], real_t, int []);
	int i,j,aux[6];
	real_t tol,**a,b[6];

	a=allocate_real_matrix(1,5,1,5);
	a[1][1]=a[1][2] = -3.0; 	a[1][3] = -18.0;
   a[1][4] = -30.0;	a[1][5]=18.0;
	a[2][2] = -1.0;	a[2][3] = -4.0;	a[2][4] = -48.0;	a[2][5]=8.0;
	a[3][3] = -6.0;	a[3][4] = -274.0;	a[3][5]=6.0;
	a[4][4]=119.0; 	a[4][5]=19.0;	a[5][5]=216.0;
	b[1]=327.0;	b[2]=291.0;	b[3]=1290.0; b[4]=275.0;	b[5]=1720.0;
	for (i=1; i<=5; i++)
		for (j=i+1; j<=5; j++) a[j][i]=a[i][j];
	tol=1.0e-6;
	decsolsym2(a,5,b,tol,aux);
	if (aux[2] == 1)
		printf("\nThe matrix is symmetric.");
	else
		printf("The matrix is asymmetric, results are meaningless.");
	printf("\nInertia : %2d,%2d,%2d\n",aux[3],aux[4],aux[5]);
	printf("\nThe computed solution:\n%10.5f%10.5f%10.5f%10.5f%10.5f\n",
			b[1],b[2],b[3],b[4],b[5]);
	free_real_matrix(a,1,5,1);
}

