#include "../real.h"
void polshtchs(int n, real_t a[])
{
	void lintfmpol(real_t, real_t, int, real_t []);
	void polchs(int, real_t []);

	lintfmpol(0.5,0.5,n,a);
	polchs(n,a);
}
