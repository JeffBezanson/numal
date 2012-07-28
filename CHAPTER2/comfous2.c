#include "../real.h"


void comfouser2(int n, real_t theta, real_t ar[], real_t ai[],
					real_t *rr, real_t *ri)
{
	void comfouser(int, real_t, real_t [], real_t *, real_t *);
	real_t car,cai,sar,sai;

	comfouser(n,theta,ar,&car,&sar);
	comfouser(n,theta,ai,&cai,&sai);
	*rr=car-sai;
	*ri=cai+sar;
}
