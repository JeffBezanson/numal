#include "../real.h"
void decsolsymtri(real_t diag[], real_t co[], int n,
						real_t aux[], real_t b[])
{
	void decsymtri(real_t [], real_t [], int, real_t []);
	void solsymtri(real_t [], real_t [], int, real_t []);

	decsymtri(diag,co,n,aux);
	if (aux[3] == n) solsymtri(diag,co,n,b);
}
