#include <stdio.h>
void main ()
{
	real_t *allocate_real_vector(int, int);
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_vector(real_t *, int);
	void free_real_matrix(real_t **, int, int, int);
	int eigcom(real_t **, real_t **, int, real_t [], real_t [],
					real_t [], real_t **, real_t **);
	int i;
	real_t **ar,**ai,**vr,**vi,*valr,*vali,em[8];

	valr=allocate_real_vector(1,4);
	vali=allocate_real_vector(1,4);
	ar=allocate_real_matrix(1,4,1,4);
	ai=allocate_real_matrix(1,4,1,4);
	vr=allocate_real_matrix(1,4,1,4);
	vi=allocate_real_matrix(1,4,1,4);

	ar[1][1]=ar[1][4]=ar[2][2]=ar[3][2]=ar[4][1]=ar[4][3]=
	ai[1][2]=ai[1][4]=ai[2][3]=ai[3][3]=ai[4][2]=1.0;
	ar[1][2]=ar[2][3]=ar[3][1]=ai[1][3]=ai[2][2]=ai[3][4]=ai[4][1]=2.0;
	ar[1][3]=ar[2][1]=ar[3][3]=ar[4][2]=
	ai[1][1]=ai[2][4]=ai[3][1]=ai[4][4]=3.0;
	ar[2][4]=ai[2][1]=ai[4][3]=4.0;
	ar[3][4]=ar[4][4]=ai[3][2]=5.0;
	em[0]=5.0e-6;  em[2]=1.0e-5;  em[4]=10.0;  em[6]=10.0;
	printf("EIGCOM: %2d\n",eigcom(ar,ai,4,em,valr,vali,vr,vi));
	printf("\n Eigenvalues:\n");
	for (i=1; i<=4; i++)
		printf(" %12.5e%12.5e * I\n",valr[i],vali[i]);
	printf("\nFirst eigenvector:\n");
	for (i=1; i<=4; i++)
		printf(" %12.5e%12.5e\n",vr[i][1],vi[i][1]);
	printf("\nEM[1]: %6.2f\nEM[3]:   %e\nEM[5]: %3.0f\nEM[7]: %3.0f\n",
			em[1],em[3],em[5],em[7]);
	free_real_vector(valr,1);
	free_real_vector(vali,1);
	free_real_matrix(ar,1,4,1);
	free_real_matrix(ai,1,4,1);
	free_real_matrix(vr,1,4,1);
	free_real_matrix(vi,1,4,1);
}

