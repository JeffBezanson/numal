#include "../real.h"
#include <stdlib.h>

int *allocate_integer_vector(int l, int u)
{
	/*  Allocates an integer vector of range [l..u].  */

	void system_error(char *);
	int *p;

	p=(int *)malloc((unsigned) (u-l+1)*sizeof(int));
	if (!p) system_error("Failure in allocate_integer_vector().");
	return p-l;
}
