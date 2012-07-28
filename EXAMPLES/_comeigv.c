#include <stdio.h>
void main ()
{
	real_t *allocate_real_vector(int, int);
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_vector(real_t *, int);
	void free_real_matrix(real_t **, int, int, int);
	int comeigval(real_t **, int, real_t [], real_t [], real_t []);
	int i,m;
	real_t **a,*re,*im,em[6];

	re=allocate_real_vector(1,3);
	im=allocate_real_vector(1,3);
	a=allocate_real_matrix(1,3,1,3);

	em[0]=1.0e-6;
	em[2]=1.0e-5;
	em[4]=30.0;
	a[1][1] = 8.0;   a[1][2] = -1.0;  a[1][3] = -5.0;
	a[2][1] = -4.0;  a[2][2] = 4.0;   a[2][3] = -2.0;
	a[3][1] = 18.0;  a[3][2] = -5.0;  a[3][3] = -7.0;
	m=comeigval(a,3,em,re,im);
	printf("The number of not calculated eigenvalues: %3.0f\n\n"
		"The eigenvalues:\n",m);
	for (i=m+1; i<=3; i++)
		printf("   %12.6e   %12.6e\n",re[i],im[i]);
	printf("\nEM[1] = %e\nEM[3] = %e\nEM[5] = %3.0f\n",
			em[1],em[3],em[5]);
	free_real_vector(re,1);
	free_real_vector(im,1);
	free_real_matrix(a,1,3,1);
}

