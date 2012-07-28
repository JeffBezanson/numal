#include "../real.h"


real_t determ(real_t **a, int n, int sign)
{
	int i;
	real_t det;

	det=1.0;
	for (i=1; i<=n; i++) det *= a[i][i];
	return (sign*fabs(det));
}
