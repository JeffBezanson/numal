#include "../real.h"


int qrisngvalbid(real_t d[], real_t b[], int n, real_t em[])
{
	int n1,k,k1,i,i1,count,max,rnk;
	real_t tol,bmax,z,x,y,g,h,f,c,s,min;

	tol=em[2]*em[1];
	count=0;
	bmax=0.0;
	max=em[4];
	min=em[6];
	rnk=n;
	do {
		k=n;
		n1=n-1;
		while (1) {
			k--;
			if (k <= 0) break;
			if (fabs(b[k]) >= tol) {
				if (fabs(d[k]) < tol) {
					c=0.0;
					s=1.0;
					for (i=k; i<=n1; i++) {
						f=s*b[i];
						b[i] *= c;
						i1=i+1;
						if (fabs(f) < tol) break;
						g=d[i1];
						d[i1]=h=sqrt(f*f+g*g);
						c=g/h;
						s = -f/h;
					}
					break;
				}
			} else {
				if (fabs(b[k]) > bmax) bmax=fabs(b[k]);
				break;
			}
		}
		if (k == n1) {
			if (d[n] < 0.0) d[n] = -d[n];
			if (d[n] <= min) rnk--;
			n=n1;
		} else {
			count++;
			if (count > max) break;
			k1=k+1;
			z=d[n];
			x=d[k1];
			y=d[n1];
			g = (n1 == 1) ? 0.0 : b[n1-1];
			h=b[n1];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=sqrt(f*f+1.0);
			f=((x-z)*(x+z)+h*(y/((f < 0.0) ? f-g : f+g)-h))/x;
			c=s=1.0;
			for (i=k1+1; i<=n; i++) {
				i1=i-1;
				g=b[i1];
				y=d[i];
				h=s*g;
				g *= c;
				z=sqrt(f*f+h*h);
				c=f/z;
				s=h/z;
				if (i1 != k1) b[i1-1]=z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y *= c;
				d[i1]=z=sqrt(f*f+h*h);
				c=f/z;
				s=h/z;
				f=c*g+s*y;
				x=c*y-s*g;
			}
			b[n1]=f;
			d[n]=x;
		}
	} while (n > 0);
	em[3]=bmax;
	em[5]=count;
	em[7]=rnk;
	return n;
}

