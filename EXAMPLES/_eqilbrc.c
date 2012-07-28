#include <math.h>
#include <stdio.h>
void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void inimat(int, int, int, int, real_t **, real_t);
	void eqilbrcom(real_t **, real_t **, int, real_t [], real_t [], int []);
	int inter[4];
	real_t **a1,**a2,em[8],d[4];

	a1=allocate_real_matrix(1,3,1,3);
	a2=allocate_real_matrix(1,3,1,3);
	em[0]=5.0e-7;
	em[6]=10.0;
	inimat(1,3,1,3,a1,0.0);
	inimat(1,3,1,3,a2,0.0);
	a1[1][1]=a1[2][2]=1.0;
	a1[3][3]=2.0;
	a2[1][3]=pow(2.0,10.0);
	a2[3][1]=pow(2.0,-10.0);
	eqilbrcom(a1,a2,3,em,d,inter);
	printf("Equilibrated matrix:\n %2.0f %2.0f %2.0f\n"
			" %2.0f %2.0f  I\n %2.0f  I %2.0f\n"
			"\nEM[7]:      %2.0f\nD[1:3]:     %2.0f %5.0f %5.0f\n"
			"INTER[1:3]: %2d %5d %5d\n",a1[1][1],a1[1][2],
			a1[1][3],a1[2][1],a1[2][2],a1[3][1],a1[3][3],em[7],
			d[1],d[2],d[3],inter[1],inter[2],inter[3]);
	free_real_matrix(a1,1,3,1);
	free_real_matrix(a2,1,3,1);
}

