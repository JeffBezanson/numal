#include "../real.h"
#include <stdlib.h>

int **allocate_integer_matrix(int lr, int ur, int lc, int uc)
{
	/*  Allocates an integer matrix of range [lr..ur][lc..uc].  */

	void system_error(char *);
	int i, **p;

	p=(int **)malloc((unsigned) (ur-lr+1)*sizeof(int*));
	if (!p) system_error("Failure in allocate_integer_matrix().");
	p -= lr;

	for (i=lr; i<=ur; i++){
		p[i]=(int *)malloc((unsigned) (uc-lc+1)*sizeof(int));
		if (!p[i]) {
            for(--i; i >= lr; i--)
                free(p[i]);
            free(p);
            system_error("Failure in allocate_integer_matrix().");
        }
		p[i] -= lc;
	}
	return p;
}
