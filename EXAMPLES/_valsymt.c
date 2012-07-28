#include <stdio.h>
void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void valsymtri(real_t [], real_t [], int, int, int,
						real_t [], real_t []);
	void vecsymtri(real_t [], real_t [], int, int, int,
						real_t [], real_t **, real_t []);
	real_t b[5],d[5],bb[4],val[3],em[10],**vec;

	vec=allocate_real_matrix(1,4,1,2);
	em[0]=1.0e-6;	em[1]=4.0;	em[2]=1.0e-5;
   em[4]=1.0e-3;	em[6]=1.0e-6;	em[8]=5.0;
	d[1]=d[2]=d[3]=d[4]=2.0;
	b[4]=0.0;	b[1]=b[2]=b[3] = -1.0;
	bb[1]=bb[2]=bb[3]=1.0;
	valsymtri(d,bb,4,1,2,val,em);
	vecsymtri(d,b,4,1,2,val,vec,em);
	printf("The eigenvalues:\n  %12.5e   %12.5e\n\nThe eigenvectors:\n"
			"  %12.5e   %12.5e\n  %12.5e   %12.5e\n  %12.5e   %12.5e\n"
			"  %12.5e   %12.5e\n"
			"\nEM[7] = %e\nEM[3] =%3.0f\nEM[5] =%3.0f\nEM[9] =%3.0f\n",
			val[1],val[2],vec[1][1],vec[1][2],vec[2][1],vec[2][2],
			vec[3][1],vec[3][2],vec[4][1],vec[4][2],
			em[7],em[3],em[5],em[9]);
	free_real_matrix(vec,1,4,1);
}

