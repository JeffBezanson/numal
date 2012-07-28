#include <stdio.h>
void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	int comvalqri(real_t **, int, real_t [], real_t [], real_t []);
	void comveches(real_t **, int, real_t, real_t,
						real_t [], real_t [], real_t []);
	int i,j,m,k;
	real_t **a,re[5],im[5],u[5],v[5],em[10];

	a=allocate_real_matrix(1,4,1,4);
	em[0]=1.0e-6;  em[2]=1.0e-6;  em[1]=4.0;  em[4]=40.0;
	em[6]=1.0e-5;  em[8]=5.0;
	for (i=1; i<=4; i++)
		for (j=1; j<=4; j++)
			a[i][j] = (i == 1) ? -1.0 : ((i-j == 1) ? 1.0 : 0.0);
	m=comvalqri(a,4,em,re,im);
	printf("The number of not calculated eigenvalues: %2d\n"
			"\nThe eigenvalues and eigenvectors:\n",m);
	for (j=m+1; j<=4; j++) {
		for (i=1; i<=4; i++)
			for (k=1; k<=4; k++)
				a[i][k] = (i == 1) ? -1.0 : ((i-k == 1) ? 1.0 : 0.0);
		comveches(a,4,re[j],im[j],em,u,v);
		printf("\n %e   %e\n\n"
			"         %12.6e   %12.6e\n         %12.6e   %12.6e\n"
			"         %12.6e   %12.6e\n         %12.6e   %12.6e\n",
			re[j],im[j],u[1],v[1],u[2],v[2],u[3],v[3],u[4],v[4]);
	}
	printf("\nEM[3] = %e\nEM[7] = %e\nEM[5] = %3.0f\nEM[9] = %3.0f\n",
			em[3],em[7],em[5],em[9]);
	free_real_matrix(a,1,4,1);
}

