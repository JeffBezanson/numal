#include <stdio.h>

void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	int hshcomcol(int, int, int, real_t **, real_t **, real_t,
						real_t *, real_t *, real_t *, real_t *);
	void hshcomprd(int, int, int, int, int, real_t **,
						real_t **, real_t **, real_t **, real_t);
	real_t k,c,s,t;
	real_t **ar,**ai;

	ar=allocate_real_matrix(1,2,1,2);
	ai=allocate_real_matrix(1,2,1,2);
	ar[1][1]=3.0;
	ar[1][2]=ar[2][1]=0.0;
	ar[2][2]=5.0;
	ai[1][1]=0.0;
	ai[1][2]=ai[2][1]=4.0;
	ai[2][2]=0.0;
	if (hshcomcol(1,2,1,ar,ai,25.0e-28,&k,&c,&s,&t))
		hshcomprd(1,2,2,2,1,ar,ai,ar,ai,t);
	printf("After using hshcomcol and hshcomprd:\n"
		" %-3.1f+%-3.1f*I  %-3.1f+%-3.1f*I\n"
		" %-3.1f+%-3.1f*I  %-3.1f+%-3.1f*I\n"
		"k, c, s, t\n %6.1f  %6.1f  %6.1f  %6.1f",
		ar[1][1],ai[1][1],ar[1][2],ai[1][2],
		ar[2][1],ai[2][1],ar[2][2],ai[2][2],k,c,s,t);
	free_real_matrix(ar,1,2,1);
	free_real_matrix(ai,1,2,1);
}

