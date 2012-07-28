#include <stdio.h>
void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void hshhrmtri(real_t **, int, real_t [], real_t [], real_t [],
						real_t [], real_t [], real_t []);
	void inimat(int, int, int, int, real_t **, real_t);
	real_t **a,d[5],b[5],bb[5],tr[4],ti[4],em[2];

	a=allocate_real_matrix(1,4,1,4);
	inimat(1,4,1,4,a,0.0);
	a[1][1]=a[2][2]=3.0;
	a[1][2]=a[3][3]=a[3][4]=a[4][4]=1.0;
	a[3][2]=2.0;
	a[4][1] = -2.0;
	em[0]=1.0e-6;
	hshhrmtri(a,4,d,b,bb,em,tr,ti);
	printf("HSHHRMTRI delivers\n\nD[1:4]:  %7.3f %7.3f %7.3f %7.3f\n"
			"B[1:3]:  %7.3f %7.3f %7.3f\n"
			"BB[1:3]: %7.3f %7.3f %7.3f\n"
			"EM[1]:   %7.3f\n",
			d[1],d[2],d[3],d[4],b[1],b[2],b[3],bb[1],bb[2],bb[3],em[1]);
	free_real_matrix(a,1,4,1);
}

