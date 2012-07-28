#include "../real.h"
#define BASE 100

void lngintmult(int u[], int v[], int product[])
{
	int lu,lv,luv,i,j,carry,t;

	lu=u[0];
	lv=v[0];
	luv=lu+lv;
	for (i=lu+1; i<=luv; i++) product[i]=0;
	for (j=lu; j>=1; j--) {
		carry=0;
		for (i=lv; i>=1; i--) {
			t=u[j]*v[i]+product[j+i]+carry;
			carry=t/BASE;
			product[j+i]=t-carry*BASE;
		}
		product[j]=carry;
	}
	if (product[1] == 0) {
		for (i=2; i<=luv; i++) product[i-1]=product[i];
		luv--;
	}
	product[0]=luv;
}
