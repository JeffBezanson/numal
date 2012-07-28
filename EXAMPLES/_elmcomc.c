#include <stdio.h>

void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void elmcomcol(int, int, int, int, real_t **, real_t **,
						real_t **, real_t **, real_t, real_t);
	real_t **ar,**ai;

	ar=allocate_real_matrix(1,2,1,2);
	ai=allocate_real_matrix(1,2,1,2);
	ar[1][1]=1.0;	ar[1][2] = -9.0;	ar[2][1]=ar[2][2] = -1.0;
	ai[1][1]=ai[1][2]=ai[2][1]=2.0;	ai[2][2] = -2.0;
	elmcomcol(1,2,2,1,ar,ai,ar,ai,1,-4);
	printf("Matrix after elimination:\n"
		"  %-3.1f+%-3.1f*I  %-3.1f+%-3.1f*I\n"
		" %-3.1f+%-3.1f*I  %-3.1f+%-3.1f*I\n",
		ar[1][1],ai[1][1],ar[1][2],ai[1][2],
		ar[2][1],ai[2][1],ar[2][2],ai[2][2]);
	free_real_matrix(ar,1,2,1);
	free_real_matrix(ai,1,2,1);
}

