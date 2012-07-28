#include "../real.h"
void decsolsym2(real_t **a, int n, real_t b[], real_t tol, int aux[])
{
	int *allocate_integer_vector(int, int);
	real_t *allocate_real_vector(int, int);
	void free_integer_vector(int *, int);
	void free_real_vector(real_t *, int);
	void decsym2(real_t **, int, real_t, int [], int [], real_t []);
	void solsym2(real_t **, int, real_t [], int [], real_t []);
	int *p;
	real_t *detaux;

	p=allocate_integer_vector(1,n);
	detaux=allocate_real_vector(1,n);
	decsym2(a,n,tol,aux,p,detaux);
	if (aux[5] == 0) solsym2(a,n,b,p,detaux);
	free_integer_vector(p,1);
	free_real_vector(detaux,1);
}
