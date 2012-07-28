#include "../real.h"
void baklbrcom(int n, int n1, int n2, real_t d[], int inter[],
					real_t **vr, real_t **vi)
{
	void baklbr(int, int, int, real_t [], int [], real_t **);

	baklbr(n,n1,n2,d,inter,vr);
	baklbr(n,n1,n2,d,inter,vi);
}
