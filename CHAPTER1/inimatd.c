#include "../real.h"
void inimatd(int lr, int ur, int shift, real_t **a, real_t x)
{
	for (; lr<=ur; lr++) a[lr][lr+shift]=x;
}
