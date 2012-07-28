#include "../real.h"
void decsoltri(real_t sub[], real_t diag[], real_t super[], int n,
					real_t aux[], real_t b[])
{
	void dectri(real_t [], real_t [], real_t [],	int, real_t []);
	void soltri(real_t [], real_t [], real_t [], int, real_t []);

	dectri(sub,diag,super,n,aux);
	if (aux[3] == n) soltri(sub,diag,super,n,b);
}
