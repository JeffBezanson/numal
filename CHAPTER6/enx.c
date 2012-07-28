#include "../real.h"


void enx(real_t x, int n1, int n2, real_t a[])
{
	if (x <= 1.5) {
		real_t ei(real_t);
      int i;
		real_t w,e;
		w = -ei(-x);
		if (n1 == 1) a[1]=w;
		if (n2 > 1) e=exp(-x);
		for (i=2; i<=n2; i++) {
			w=(e-x*w)/(i-1);
			if (i >= n1) a[i]=w;
		}
	} else {
		int i,n;
		real_t w,e,an;
		n=ceil(x);
		if (n <= 10) {
			real_t f,w1,t,h,p[20];
			p[2] =0.37534261820491e-1;  p[11]=0.135335283236613;
			p[3] =0.89306465560228e-2;  p[12]=0.497870683678639e-1;
			p[4] =0.24233983686581e-2;  p[13]=0.183156388887342e-1;
			p[5] =0.70576069342458e-3;  p[14]=0.673794699908547e-2;
			p[6] =0.21480277819013e-3;  p[15]=0.247875217666636e-2;
			p[7] =0.67375807781018e-4;  p[16]=0.911881965554516e-3;
			p[8] =0.21600730159975e-4;  p[17]=0.335462627902512e-3;
			p[9] =0.70411579854292e-5;  p[18]=0.123409804086680e-3;
			p[10]=0.23253026570282e-5;  p[19]=0.453999297624848e-4;
			f=w=p[n];
			e=p[n+9];
			w1=t=1.0;
			h=x-n;
			i=n-1;
			do {
				f=(e-i*f)/n;
				t = -h*t/(n-i);
				w1=t*f;
				w += w1;
				i--;
			} while (fabs(w1) > 1.0e-15*w);
		} else {
			real_t *allocate_real_vector(int, int);
			void free_real_vector(real_t *, int);
			void nonexpenx(real_t, int, int, real_t []);
			real_t *b;
			b=allocate_real_vector(n,n);
			nonexpenx(x,n,n,b);
			w=b[n]*exp(-x);
			free_real_vector(b,n);
		}
		if (n1 == n2 && n1 == n)
			a[n]=w;
		else {
			e=exp(-x);
			an=w;
			if (n <= n2 && n >= n1) a[n]=w;
			for (i=n-1; i>=n1; i--) {
				w=(e-i*w)/x;
				if (i <= n2) a[i]=w;
			}
			w=an;
			for (i=n+1; i<=n2; i++) {
				w=(e-x*w)/(i-1);
				if (i >= n1) a[i]=w;
			}
		}
	}
}
