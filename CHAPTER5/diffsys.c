#include "../real.h"


void diffsys(real_t *x, real_t xe, int n, real_t y[],
				void (*derivative)(int, real_t, real_t [], real_t []),
				real_t aeta, real_t reta, real_t s[], real_t h0,
				void (*output)(int, real_t, real_t, real_t [], real_t []))
{
	real_t *allocate_real_vector(int, int);
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_vector(real_t *, int);
	void free_real_matrix(real_t **, int, int, int);
	int i,j,k,kk,jj,l,m,r,sr,konv,b0,bh,last,next;
	real_t a,b,b1,c,g,h,u,v,ta,fc,*ya,*yl,*ym,*dy,*dz,**dt,d[7],
			**yg,**yh;

	ya=allocate_real_vector(1,n);
	yl=allocate_real_vector(1,n);
	ym=allocate_real_vector(1,n);
	dy=allocate_real_vector(1,n);
	dz=allocate_real_vector(1,n);
	dt=allocate_real_matrix(1,n,0,6);
	yg=allocate_real_matrix(0,7,1,n);
	yh=allocate_real_matrix(0,7,1,n);

	last=0;
	h=h0;
	do {
		next=0;
		if (h*1.1 >= xe-(*x)) {
			last=1;
			h0=h;
			h=xe-(*x)+FLT_EPSILON;
		}
		(*derivative)(n,*x,y,dz);
		bh=0;
		for (i=1; i<=n; i++) ya[i]=y[i];
		while (1) {
			a=h+(*x);
			fc=1.5;
			b0=0;
			m=1;
			r=2;
			sr=3;
			jj = -1;
			for (j=0; j<=9; j++) {
				if (b0) {
					d[1]=16.0/9.0;
					d[3]=64.0/9.0;
					d[5]=256.0/9.0;
				} else {
					d[1]=9.0/4.0;
					d[3]=9.0;
					d[5]=36.0;
				}
				konv=1;
				if (j > 6) {
					l=6;
					d[6]=64.0;
					fc *= 0.6;
				} else {
					l=j;
					d[l]=m*m;
				}
				m *= 2;
				g=h/m;
				b=g*2.0;
				if (bh && j < 8)
					for (i=1; i<=n; i++) {
						ym[i]=yh[j][i];
						yl[i]=yg[j][i];
					}
				else {
					kk=(m-2)/2;
					m--;
					for (i=1; i<=n; i++) {
						yl[i]=ya[i];
						ym[i]=ya[i]+g*dz[i];
					}
					for (k=1; k<=m; k++) {
						(*derivative)(n,(*x)+k*g,ym,dy);
						for (i=1; i<=n; i++) {
							u=yl[i]+b*dy[i];
							yl[i]=ym[i];
							ym[i]=u;
							u=fabs(u);
							if (u > s[i]) s[i]=u;
						}
						if (k == kk && k != 2) {
							jj++;
							for (i=1; i<=n; i++) {
								yh[jj][i]=ym[i];
								yg[jj][i]=yl[i];
							}
						}
					}
				}
				(*derivative)(n,a,ym,dy);
				for (i=1; i<=n; i++) {
					v=dt[i][0];
					ta=c=dt[i][0]=(ym[i]+yl[i]+g*dy[i])/2.0;
					for (k=1; k<=l; k++) {
						b1=d[k]*v;
						b=b1-c;
						u=v;
						if (b != 0.0) {
							b=(c-v)/b;
							u=c*b;
							c=b1*b;
						}
						v=dt[i][k];
						dt[i][k]=u;
						ta += u;
					}
					if (fabs(y[i]-ta) > reta*s[i]+aeta) konv=0;
					y[i]=ta;
				}
				if (konv) {
					next=1;
					break;
				}
				d[2]=4.0;
				d[4]=16.0;
				b0 = !b0;
				m=r;
				r=sr;
				sr=m*2;
			}
			if (next) break;
			bh = !bh;
			last=0;
			h /= 2.0;
		}
		h *= fc;
		(*x)=a;
		(*output)(n,*x,xe,y,s);
	} while (!last);
	free_real_vector(ya,1);
	free_real_vector(yl,1);
	free_real_vector(ym,1);
	free_real_vector(dy,1);
	free_real_vector(dz,1);
	free_real_matrix(dt,1,n,0);
	free_real_matrix(yg,0,7,1);
	free_real_matrix(yh,0,7,1);
}

