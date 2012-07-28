#include "../real.h"
void commul(real_t ar, real_t ai, real_t br, real_t bi, real_t *rr, real_t *ri)
{
	*rr=ar*br-ai*bi;
	*ri=ar*bi+ai*br;
}
