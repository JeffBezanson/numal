#include <stdio.h>
void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void sclcom(real_t **, real_t **, int, int, int);
	void inimat(int, int, int, int, real_t **, real_t);
	void eighrm(real_t **, int, int, real_t [], real_t **, real_t **,
					real_t []);
	int i;
	real_t **a,**vecr,**veci,val[2],em[10];

	a=allocate_real_matrix(1,4,1,4);
	vecr=allocate_real_matrix(1,4,1,1);
	veci=allocate_real_matrix(1,4,1,1);

	inimat(1,4,1,4,a,0.0);
	a[1][1]=a[2][2]=3.0;
	a[1][2]=a[3][3]=a[3][4]=a[4][4]=1.0;
	a[3][2]=2.0;  a[4][1] = -2.0;
	em[0]=5.0e-6;
	em[2]=1.0e-5;
	em[4]=0.01;
	em[6]=1.0e-5;
	em[8]=5.0;
	eighrm(a,4,1,val,vecr,veci,em);
	sclcom(vecr,veci,4,1,1);
	printf("Largest eigenvalue: %e\n"
			"\nCorresponding eigenvector:\n",val[1]);
	for (i=1; i<=4; i++)
		printf("  %5.2f %6.3f*I\n",vecr[i][1],veci[i][1]);
	printf("\nEM[1]:  %e\nEM[3]:  %e\nEM[5]:  %e\nEM[7]:  %e"
			"\nEM[9]:  %e\n",em[1],em[3],em[5],em[7],em[9]);
	free_real_matrix(a,1,4,1);
	free_real_matrix(vecr,1,4,1);
	free_real_matrix(veci,1,4,1);
}

