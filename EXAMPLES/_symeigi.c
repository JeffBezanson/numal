#include <stdio.h>
void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void symeigimp(int, real_t **, real_t **, real_t [],
						real_t [], real_t [], real_t []);
	int qrisym(real_t **, int, real_t [], real_t []);
	int i,j;
	real_t **a,**x,val[5],lbound[5],ubound[5],em[6],aux[6];

	a=allocate_real_matrix(1,4,1,4);
	x=allocate_real_matrix(1,4,1,4);
	a[1][1]=a[2][2]=a[3][3]=a[4][4]=6.0;
	a[1][2]=a[2][1]=a[3][1]=a[1][3]=4.0;
	a[4][2]=a[2][4]=a[3][4]=a[4][3]=4.0;
	a[1][4]=a[4][1]=a[3][2]=a[2][3]=1.0;
	for (i=1; i<=4; i++)
		for (j=i; j<=4; j++) x[i][j]=x[j][i]=a[i][j];
	em[0]=1.0e-6;	em[4]=100.0;	em[2]=1.0e-5;
	qrisym(x,4,val,em);
	aux[0]=0.0; 	aux[4]=10.0;	aux[2]=1.0e-6;
	symeigimp(4,a,x,val,lbound,ubound,aux);
	printf("\nThe exact eigenvalues are:  -1,  5,  5,  15\n\n"
		"The computed eigenvalues:\n %12.6e\n %12.6e\n %12.6e\n"
		" %12.6e\n\n Lowerbounds   Upperbounds\n %e   %e\n"
		" %e   %e\n %e   %e\n %e   %e\n"
		"\nNumber of iterations =%3.0f\n"
		"Infinity norm of A   =%3.0f\n"
		"Maximum absolute element of residu = %e",
		val[1],val[2],val[3],val[4],
		lbound[1],ubound[1],lbound[2],ubound[2],
		lbound[3],ubound[3],lbound[4],ubound[4],
		aux[5],aux[1],aux[3]);
	free_real_matrix(a,1,4,1);
	free_real_matrix(x,1,4,1);
}

