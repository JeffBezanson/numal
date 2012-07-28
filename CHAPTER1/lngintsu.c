#include "../real.h"
#define BASE 100

void lngintsubtract(int u[], int v[], int difference[])
{
	int lu,lv,diff,i,t,j,carry;

	lu=u[0];
	lv=v[0];
	if ((lu < lv) || ((lu == lv) && (u[1] < v[1]))) {
		difference[0]=0;
		return;
	}
	diff=lu-lv+1;
	carry=0;
	for (i=lu; i>=diff; i--) {
		t=u[i]-v[i-diff+1]+carry;
		carry = (t < 0) ? -1 : 0;
		difference[i]=t-carry*BASE;
	}
	for (i=diff-1; i>=1; i--) {
		t=u[i]+carry;
		carry = (t < 0) ? -1 : 0;
		difference[i]=t-carry*BASE;
	}
	if (carry == -1) {
		difference[0]=0;
		return;
	}
	i=1;
	j=lu;
	while ((difference[i] == 0) && (j > 1)) {
		j--;
		i++;
	}
	difference[0]=j;
	if (j < lu)
		for (i=1; i<=j; i++) difference[i]=difference[lu+i-j];
}
