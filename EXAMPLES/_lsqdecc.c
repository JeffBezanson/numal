#include <stdio.h>
void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void lsqdecomp(real_t **, int, int, int, real_t [], real_t [], int []);
	void lsqrefsol(real_t **, real_t **, int, int, int, real_t [], real_t [],
						int [], real_t [], real_t *, real_t [], real_t []);
	int n,m,n1,i,j,ci[4];
	real_t ldx,**a,**qr,aux[8],b[6],res[6],aid[4],x[4];

	a=allocate_real_matrix(1,5,1,3);
	qr=allocate_real_matrix(1,5,1,3);
	n=5;  m=3;  n1=1;
	a[1][1]=1.0;  a[1][2]=1000.0;  a[1][3]=5.0;
	a[2][1]=1.0;  a[2][2]=0.0;     a[2][3]=8.0;
	a[3][1]=0.0;  a[3][2]=3.0;     a[3][3]=2.0;
	a[4][1]=1.0;  a[4][2]=2.0;     a[4][3]=1.0e-5;
	a[5][1]=a[5][2]=a[5][3]=0.0;
	b[1]=2016.0;  b[2]=25.0;  b[3]=12.0;  b[4]=5.00003;  b[5]=1.0;
	aux[2]=1.0e-6;
	aux[6]=5.0;
	for (i=1; i<=5; i++)
		for (j=1; j<=3; j++) qr[i][j]=a[i][j];
	lsqdecomp(qr,n,m,n1,aux,aid,ci);
	lsqrefsol(a,qr,n,m,n1,aux,aid,ci,b,&ldx,x,res);
	printf("The solution vector:\n   %e   %e   %e\n"
			"\nThe residual vector:\n %e\n  %e\n  %e\n  %e\n"
			"Number of iterations: %3.0f\n"
			"Norm last correction of x:  %e\n",
			x[1],x[2],x[3],res[2],res[3],res[4],res[5],aux[7],ldx);
	free_real_matrix(a,1,5,1);
	free_real_matrix(qr,1,5,1);
}

