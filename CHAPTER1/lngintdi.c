#include "../real.h"
#define BASE 100

void lngintdivide(int u[], int v[], int quotient[], int remainder[])
{
	int *allocate_integer_vector(int, int);
	void free_integer_vector(int *, int);
	int lu,lv,v1,diff,i,t,scale,d,q1,j,carry,*uu,*a;

	lu=u[0];
	lv=v[0];
	v1=v[1];
	diff=lu-lv;

	if (lv == 1) {
		carry=0;
		for (i=1; i<=lu; i++) {
			t=carry*BASE+u[i];
			quotient[i]=t/v1;
			carry=t-quotient[i]*v1;
		}
		remainder[0]=1;
		remainder[1]=carry;
		if (quotient[1] == 0) {
			for (i=2; i<=lu; i++) quotient[i-1]=quotient[i];
			quotient[0]=lu - ((lu == 1) ? 0 : 1);
		} else
			quotient[0]=lu;
		return;
	}

	if (lu < lv) {
		quotient[0]=1;
		quotient[1]=0;
		for (i=0; i<=lu; i++) remainder[i]=u[i];
		return;
	}

	uu=allocate_integer_vector(0,lu);
	a=allocate_integer_vector(0,lv);
	for (i=0; i<=lu; i++) uu[i]=u[i];
	scale=BASE/(v1+1);
	if (scale > 1) {
		/* normalize u */
		carry=0;
		for (i=lu; i>=1; i--) {
			t=scale*uu[i]+carry;
			carry=t/BASE;
			uu[i]=t-carry*BASE;
		}
		uu[0]=carry;
		/* normalize v */
		carry=0;
		for (i=lv; i>=1; i--) {
			t=scale*v[i]+carry;
			carry=t/BASE;
			v[i]=t-carry*BASE;
		}
		v1=v[1];
	} else
		uu[0]=0;

	/* compute quotient and remainder */
	for (i=0; i<=diff; i++) {
		d=uu[i]*BASE+uu[i+1];
		q1 = (uu[i] == v1) ? BASE-1 : d/v1;
		if (v[2]*q1 > (d-q1*v1)*BASE+uu[i+2]) {
			q1--;
			if (v[2]*q1 > (d-q1*v1)*BASE+uu[i+2]) q1--;
		}
		/* uu[i:i+lv]=u[i:i+lv]-q1*v[1:lv] */
		carry=0;
		for (j=lv; j>=1; j--) {
			t=q1*v[j]+carry;
			carry=t/BASE;
			a[j]=t-carry*BASE;
		}
		a[0]=carry;
		carry=0;
		for (j=lv; j>=0; j--) {
			t=uu[i+j]-a[j]+carry;
			carry = (t < 0) ? -1 : 0;
			uu[i+j]=t-carry*BASE;
		}
		/* if carry=-1 then q1 is one too large,
			and v must be added back to uu[i:i+lv] */
		if (carry == -1) {
			q1--;
			carry=0;
			for (j=lv; j>=1; j--) {
				t=uu[i+j]+v[j]+carry;
				carry = (t < BASE) ? 0 :1;
				uu[i+j]=t-carry*BASE;
			}
		}
		quotient[i]=q1;
	}

	/* correct storage of quotient */
	if (quotient[0] != 0) {
		for (i=diff; i>=0; i--) quotient[i+1]=quotient[i];
		quotient[0]=diff+1;
	} else
		if (quotient[1] != 0)
			quotient[0]=diff;
		else {
			for (i=1; i<=diff-1; i++) quotient[i]=quotient[i+1];
			quotient[0]=diff-1;
		}

	/* remainder=uu[diff+1:lu]/scale */
	if (scale > 1) {
		carry=0;
		for (i=1; i<=lv; i++) {
			t=carry*BASE+uu[diff+i];
			remainder[i]=t/scale;
			carry=t-remainder[i]*scale;
		}
	} else
		for (i=1; i<=lv; i++) remainder[i]=uu[diff+i];

	/* correct storage of remainder */
	i=1;
	j=lv;
	while (remainder[i] == 0 && j > 1) {
		j--;
		i++;
	}
	remainder[0]=j;
	if (j < lv)
		for (i=1; i<=j; i++) remainder[i]=remainder[lv+i-j];

	/* unnormalize the divisor v */
	if (scale > 1) {
		carry=0;
		for (i=1; i<=lv; i++) {
			t=carry*BASE+v[i];
			v[i]=t/scale;
			carry=t-v[i]*scale;
		}
	}

	free_integer_vector(uu,0);
	free_integer_vector(a,0);
}
