#include <stdio.h>
void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void qzi(int, real_t **, real_t **, real_t **, real_t [],
				real_t [], real_t [], int [], real_t []);
	int k,l,iter[5];
	real_t **a,**b,**x,alfr[5],alfi[5],beta[5],em[2];

	a=allocate_real_matrix(1,4,1,4);
	b=allocate_real_matrix(1,4,1,4);
	x=allocate_real_matrix(1,4,1,4);
	a[1][1]=2.0;  a[1][2]=3.0;     a[1][3] = -3.0;   a[1][4]=4.0;
	a[2][1]=1.0;  a[2][2] = -1.0;  a[2][3]=5.0;      a[2][4]=1.0;
	a[3][1]=0.0;  a[3][2]=2.0;     a[3][3]=6.0;      a[3][4]=8.0;
	a[4][1]=1.0;  a[4][2]=1.0;     a[4][3]=0.0;      a[4][4]=4.0;
	b[1][1]=1.0;  b[1][2]=5.0;     b[1][3]=9.0;      b[1][4]=0.0;
	b[2][1]=2.0;  b[2][2]=6.0;     b[2][3]=10.0;     b[2][4]=2.0;
	b[3][1]=3.0;  b[3][2]=7.0;     b[3][3]=11.0;     b[3][4] = -1.0;
	b[4][1]=4.0;  b[4][2]=8.0;     b[4][3]=12.0;     b[4][4]=3.0;
	for (k=1; k<=4; k++)
		for (l=1; l<=4; l++)	x[k][l] = (k == l) ? 1.0 : 0.0;
	em[0]=1.0e-35;
	em[1]=1.0e-6;
	qzi(4,a,b,x,alfr,alfi,beta,iter,em);
	for (k=1; k<=4; k++)
		printf("ITER[%1d]=%3d\n",k,iter[k]);
	printf("\nEigenvectors:\n");
	for (k=1; k<=4; k++)
		printf(" %12.6e  %12.6e  %12.6e  %12.6e\n",
				x[k][1],x[k][2],x[k][3],x[k][4]);
	printf("\nALFA(real part)    ALFA(imaginary part)      BETA\n");
	for (k=1; k<=4; k++)
		printf(" %12.6e  %16.6e  %21.6e\n",alfr[k],alfi[k],beta[k]);
	printf("\nLAMBDA(real part)  LAMBDA(imaginary part)\n");
	for (k=1; k<=4; k++)
		if (beta[k] == 0.0)
			printf("  INFINITE          INDEFINITE\n");
		else
			printf(" %12.6e  %16.6e\n",
					alfr[k]/beta[k],alfi[k]/beta[k]);
	free_real_matrix(a,1,4,1);
	free_real_matrix(b,1,4,1);
	free_real_matrix(x,1,4,1);
}

