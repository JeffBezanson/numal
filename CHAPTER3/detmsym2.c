#include "../real.h"
real_t determsym2(real_t detaux[], int n, int aux[])
{
	int i;
	real_t det;

	if (aux[5] > 0)
		det=0.0;
	else {
		det=1.0;
		for (i=1; i<=n; i++) det *= detaux[i];
	}
	return (det);
}
