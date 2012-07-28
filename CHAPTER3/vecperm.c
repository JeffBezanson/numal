#include "../real.h"
void vecperm(int perm[], int low, int upp, real_t vector[])
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
			a=vector[k];
			j=perm[k];
			while (j != t) {
				vector[k]=vector[j];
				todo[k]=0;
				k=j;
				j=perm[k];
			}
			vector[k]=a;
			todo[k]=0;
		}
	free_integer_vector(todo,low);
}

