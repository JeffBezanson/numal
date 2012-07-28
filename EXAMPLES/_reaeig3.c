#include <stdio.h>
void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	int reaeig3(real_t **, int, real_t [], real_t [], real_t **);
	int i,j,m;
	real_t **a,**vec,val[5],em[6];

	a=allocate_real_matrix(1,4,1,4);
	vec=allocate_real_matrix(1,4,1,4);
	for (i=1; i<=4; i++)
		for (j=1; j<=4; j++)
			a[i][j] = (i == 1) ? 1.0 : 1.0/(i+j-1);
	em[0]=1.0e-6;
	em[2]=1.0e-5;
	em[4]=40.0;
	m=reaeig3(a,4,em,val,vec);
	printf("The number of not calculated eigenvalues: %3.0f\n\n"
		"The eigenvalues and corresponding eigenvectors:\n",m);
	for (i=m+1; i<=4; i++)
		printf("\n %12.6e   %12.6e\n                %12.6e"
				"\n                %12.6e\n                %12.6e\n",
				val[i],vec[1][i],vec[2][i],vec[3][i],vec[4][i]);
	printf("\nEM[1] = %e\nEM[3] = %e\nEM[5] = %3.0f\n",
			em[1],em[3],em[5]);
	free_real_matrix(a,1,4,1);
	free_real_matrix(vec,1,4,1);
}

