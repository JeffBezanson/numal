#include "../real.h"


void nonexpenx(real_t x, int n1, int n2, real_t a[])
{
	int i,n;
	real_t w,an;

	n = (x <= 1.5) ? 1 : ceil(x);
	if (n <= 10) {
		real_t *allocate_real_vector(int, int);
		void free_real_vector(real_t *, int);
		void enx(real_t, int, int, real_t []);
		real_t *b;
		b=allocate_real_vector(n,n);
		enx(x,n,n,b);
		w=b[n]*exp(x);
		free_real_vector(b,n);
	} else {
		int k,k1;
		real_t ue,ve,we,we1,uo,vo,wo,wo1,r,s;
		ue=1.0;
		ve=we=1.0/(x+n);
		we1=0.0;
		uo=1.0;
		vo = -n/(x*(x+n+1.0));
		wo1=1.0/x;
		wo=vo+wo1;
		w=(we+wo)/2.0;
		k1=1;
		k=k1;
		while (wo-we > 1.0e-15*w && we > we1 && wo < wo1) {
			we1=we;
			wo1=wo;
			r=n+k;
			s=r+x+k;
			ue=1.0/(1.0-k*(r-1.0)*ue/((s-2.0)*s));
			uo=1.0/(1.0-k*r*uo/(s*s-1.0));
			ve *= (ue-1.0);
			vo *= (uo-1.0);
			we += ve;
			wo += vo;
			w=(we+wo)/2.0;
			k1++;
			k=k1;
		}
	}
	an=w;
	if (n <= n2 && n >= n1) a[n]=w;
	for (i=n-1; i>=n1; i--) {
		w=(1.0-i*w)/x;
		if (i <= n2) a[i]=w;
	}
	w=an;
	for (i=n+1; i<=n2; i++) {
		w=(1.0-x*w)/(i-1);
		if (i >= n1) a[i]=w;
	}
}
