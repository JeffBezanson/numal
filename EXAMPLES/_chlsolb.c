#include <stdio.h>
void main ()
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void chldecbnd(real_t [], int, int, real_t []);
	void chlsolbnd(real_t [], int, int, real_t []);
	real_t chldetermbnd(real_t [], int, int);
	int i;
	real_t *symband,*right,*aux;

	symband=allocate_real_vector(1,9);
	right=allocate_real_vector(1,5);
	aux=allocate_real_vector(2,3);

	for (i=1; i<=9; i++)
		symband[i] = ((i/2)*2 < i) ? 2.0 : -1.0;
	right[1]=right[5]=1.0;
	right[2]=right[3]=right[4]=0.0;
	aux[2]=1.0e-12;
	chldecbnd(symband,5,1,aux);
	if (aux[3] == 5) {
		chlsolbnd(symband,5,1,right);
		printf("Delivers: %8.4f %8.4f %8.4f %8.4f %8.4f\n"
				"Determinant is  %e\n",right[1],right[2],right[3],
				right[4],right[5],chldetermbnd(symband,5,1));
	}
	free_real_vector(symband,1);
	free_real_vector(right,1);
	free_real_vector(aux,2);
}

