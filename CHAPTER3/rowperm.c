#include "../real.h"
void rowperm(int perm[], int low, int upp, int i, real_t **mat)
{
	int *allocate_integer_vector(int, int);
	void free_integer_vector(int *, int);
	int t,j,k,*todo;
	real_t a;

	todo=allocate_integer_vector(low,upp);
	for (t=low; t<=upp; t++) todo[t]=1;
	for (t=low; t<=upp; t++)
		if (todo[t]) {
			k=t;
			a=mat[i][k];
			j=perm[k];
			while (j != t) {
				mat[i][k]=mat[i][j];
				todo[k]=0;
				k=j;
				j=perm[k];
			}
			mat[i][k]=a;
			todo[k]=0;
		}
	free_integer_vector(todo,low);
}

