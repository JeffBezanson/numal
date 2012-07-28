#include <stdio.h>
void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void sclcom(real_t **, real_t **, int, int, int);
	void inimat(int, int, int, int, real_t **, real_t);
	int qrihrm(real_t **, int, real_t [], real_t **, real_t **, real_t []);
	int i;
	real_t **a,**vr,**vi,val[5],em[6];

	a=allocate_real_matrix(1,4,1,4);
	vr=allocate_real_matrix(1,4,1,4);
	vi=allocate_real_matrix(1,4,1,4);

	inimat(1,4,1,4,a,0.0);
	a[1][1]=a[2][2]=3.0;
	a[3][2]=2.0;  a[4][1] = -2.0;
	a[1][2]=a[3][3]=a[3][4]=a[4][4]=1.0;
	em[0]=em[2]=5.0e-5;
	em[4]=20.0;
	printf("QRIHRM: %2d\n",qrihrm(a,4,val,vr,vi,em));
	sclcom(vr,vi,4,2,3);
	printf("\nEigenvalues:\n VAL[1]: %7.3f\n VAL[2]: %7.3f\n"
		" VAL[3]: %7.3f\n VAL[4]: %7.3f\n"
		"\nEigenvectors corresponding to\n   VAL[2]  ,    VAL[3]\n",
		val[1],val[2],val[3],val[4]);
	for (i=1; i<=4; i++)
		printf(" %3.0f%3.0f*I  ,  %3.0f%3.0f*I\n",
		vr[i][2],vi[i][2],vr[i][3],vi[i][3]);
	printf("\nEM[1]:  %e\nEM[3]:  %e\nEM[5]:  %e\n",
			em[1],em[3],em[5]);
	free_real_matrix(a,1,4,1);
	free_real_matrix(vr,1,4,1);
	free_real_matrix(vi,1,4,1);
}

