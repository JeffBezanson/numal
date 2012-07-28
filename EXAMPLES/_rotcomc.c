#include <stdio.h>

void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void rotcomcol(int, int, int, int, real_t **, real_t **,
						real_t, real_t, real_t);
	int i,j;
	real_t **ar,**ai;

	ar=allocate_real_matrix(1,2,1,2);
	ai=allocate_real_matrix(1,2,1,2);
	ar[1][1]=4.0;	ar[1][2]=5.0;	ar[2][1] = -5.0;	ar[2][2]=4.0;
	ai[1][1]=3.0;	ai[1][2]=ai[2][1]=0.0;	ai[2][2] = -3.0;
	rotcomcol(1,2,1,2,ar,ai,0.08,0.06,-0.1);
	printf("After postmultiplication:\n"
		" %+3.1f%+3.1f*I  %+3.1f%+3.1f*I\n"
		" %+3.1f%+3.1f*I  %+3.1f%+3.1f*I\n",
		ar[1][1],ai[1][1],ar[1][2],ai[1][2],
		ar[2][1],ai[2][1],ar[2][2],ai[2][2]);
	free_real_matrix(ar,1,2,1);
	free_real_matrix(ai,1,2,1);
}

