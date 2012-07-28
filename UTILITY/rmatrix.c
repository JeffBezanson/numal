#include "../real.h"
#include <stdlib.h>

real_t **allocate_real_matrix(int lr, int ur, int lc, int uc)
{
	/*  Allocates a real matrix of range [lr..ur][lc..uc].  */

	void system_error(char *);
	int i;
	real_t **p;

	p=(real_t **)malloc((unsigned) (ur-lr+1)*sizeof(real_t*));
	if (!p) system_error("Failure in allocate_real_matrix().");
	p -= lr;

	for (i=lr; i<=ur; i++){
		p[i]=(real_t *)malloc((unsigned) (uc-lc+1)*sizeof(real_t));
		if (!p[i]) {
            for(--i; i>=lr; i--)
                free(p[i]);
            free(p);
            system_error("Failure in allocate_real_matrix().");
        }
		p[i] -= lc;
	}
	return p;
}
