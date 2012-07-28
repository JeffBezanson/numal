#include "../real.h"
void shtchspol(int n, real_t a[])
{
	void chspol(int, real_t []);
	void lintfmpol(real_t, real_t, int, real_t []);

	chspol(n,a);
	lintfmpol(2.0,-1.0,n,a);
}
