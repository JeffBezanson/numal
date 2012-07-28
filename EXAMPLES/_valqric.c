#include <stdio.h>
void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void inimat(int, int, int, int, real_t **, real_t);
	int valqricom(real_t **, real_t **, real_t [], int, real_t [],
						real_t [], real_t []);
	int i;
	real_t **a1,**a2,b[4],val1[5],val2[5],em[6];

	a1=allocate_real_matrix(1,4,1,4);
	a2=allocate_real_matrix(1,4,1,4);
	inimat(1,4,1,4,a1,0.0);
	inimat(1,4,1,4,a2,0.0);
	a1[1][1] = -4.0;  a1[1][2] = -5.0;
	a1[1][3] = a2[1][1] = a2[1][4] = -2.0;
	a2[1][2] = a2[1][3] = -6.0;
	b[1]=b[2]=b[3]=1.0;
	em[0]=5.0e-6;  em[1]=27.0;  em[2]=1.0e-6;  em[4]=15.0;
	printf("VALQRICOM: %2d\n",valqricom(a1,a2,b,4,em,val1,val2));
	printf("\n Eigenvalues:\n    Real part     Imaginary part\n");
	for (i=1; i<=4; i++)
		printf(" %12.4e    %12.4e\n",val1[i],val2[i]);
	printf("\nEM[3]: %e\nEM[5]: %3.0f\n",em[3],em[5]);
	free_real_matrix(a1,1,4,1);
	free_real_matrix(a2,1,4,1);
}

