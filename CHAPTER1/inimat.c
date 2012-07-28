#include "../real.h"
void inimat(int lr, int ur, int lc, int uc, real_t **a, real_t x)
{
	int j;

	for (; lr<=ur; lr++)
		for (j=lc; j<=uc; j++) a[lr][j]=x;
}
