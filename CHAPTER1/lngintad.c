#include "../real.h"
#define BASE 100

void lngintadd(int u[], int v[], int sum[])
{
	int lu,lv,diff,carry,i,t,max;

	lu=u[0];
	lv=v[0];
	if (lu >= lv) {
		max=lu;
		diff=lu-lv+1;
		carry=0;
		for (i=lu; i>=diff; i--) {
			t=u[i]+v[i-diff+1]+carry;
			carry = (t < BASE) ? 0 : 1;
			sum[i]=t-carry*BASE;
		}
		for (i=diff-1; i>=1; i--) {
			t=u[i]+carry;
			carry = (t < BASE) ? 0 : 1;
			sum[i]=t-carry*BASE;
		}
	} else {
		max=lv;
		diff=lv-lu+1;
		carry=0;
		for (i=lv; i>=diff; i--) {
			t=v[i]+u[i-diff+1]+carry;
			carry = (t < BASE) ? 0 : 1;
			sum[i]=t-carry*BASE;
		}
		for (i=diff-1; i>=1; i--) {
			t=v[i]+carry;
			carry = (t < BASE) ? 0 : 1;
			sum[i]=t-carry*BASE;
		}
	}
	if (carry == 1) {
		for (i=max; i>=1; i--) sum[i+1]=sum[i];
		sum[1]=1;
		max=max+1;
	}
	sum[0]=max;
}
