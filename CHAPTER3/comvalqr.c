#include "../real.h"


int comvalqri(real_t **a, int n, real_t em[], real_t re[], real_t im[])
{
	int i,j,p,q,max,count,n1,p1,p2,imin1,i1,i2,i3,b;
	real_t disc,sigma,rho,g1,g2,g3,psi1,psi2,aa,e,k,s,norm,machtol2,
			tol,w;

	norm=em[1];
	w=em[0]*norm;
	machtol2=w*w;
	tol=em[2]*norm;
	max=em[4];
	count=0;
	w=0.0;
	do {
		i=n;
		do {
			q=i;
			i--;
		} while ((i >= 1) ? (fabs(a[i+1][i]) > tol) : 0);
		if (q > 1)
			if (fabs(a[q][q-1]) > w) w=fabs(a[q][q-1]);
		if (q >= n-1) {
			n1=n-1;
			if (q == n) {
				re[n]=a[n][n];
				im[n]=0.0;
				n=n1;
			} else {
				sigma=a[n][n]-a[n1][n1];
				rho = -a[n][n1]*a[n1][n];
				disc=sigma*sigma-4.0*rho;
				if (disc > 0.0) {
					disc=sqrt(disc);
					s = -2.0*rho/(sigma+((sigma >= 0.0) ? disc : -disc));
					re[n]=a[n][n]+s;
					re[n1]=a[n1][n1]-s;
					im[n]=im[n1]=0.0;
				} else {
					re[n]=re[n1]=(a[n1][n1]+a[n][n])/2.0;
					im[n1]=sqrt(-disc)/2.0;
					im[n] = -im[n1];
				}
				n -= 2;
			}
		} else {
			count++;
			if (count > max) break;
			n1=n-1;
			sigma=a[n][n]+a[n1][n1]+sqrt(fabs(a[n1][n-2]*a[n][n1])*em[0]);
			rho=a[n][n]*a[n1][n1]-a[n][n1]*a[n1][n];
			i=n-1;
			do {
				p1=i1=i;
				i--;
			} while ((i-1 >= q) ? (fabs(a[i][i-1]*a[i1][i]*
					(fabs(a[i][i]+a[i1][i1]-sigma)+fabs(a[i+2][i1]))) >
					fabs(a[i][i]*((a[i][i]-sigma)+
						a[i][i1]*a[i1][i]+rho))*tol) : 0);
			p=p1-1;
			p2=p+2;
			for (i=p; i<=n-1; i++) {
				imin1=i-1;
				i1=i+1;
				i2=i+2;
				if (i == p) {
					g1=a[p][p]*(a[p][p]-sigma)+a[p][p1]*a[p1][p]+rho;
					g2=a[p1][p]*(a[p][p]+a[p1][p1]-sigma);
					if (p1 <= n1) {
						g3=a[p1][p]*a[p2][p1];
						a[p2][p]=0.0;
					} else
						g3=0.0;
				} else {
					g1=a[i][imin1];
					g2=a[i1][imin1];
					g3 = (i2 <= n) ? a[i2][imin1] : 0.0;
				}
				k = (g1 >= 0.0) ? sqrt(g1*g1+g2*g2+g3*g3) :
										-sqrt(g1*g1+g2*g2+g3*g3);
				b = (fabs(k) > machtol2);
				aa = (b ? g1/k+1.0 : 2.0);
				psi1 = (b ? g2/(g1+k) : 0.0);
				psi2 = (b ? g3/(g1+k) : 0.0);
				if (i != q)
					a[i][imin1] = (i == p) ? -a[i][imin1] : -k;
				for (j=i; j<=n; j++) {
					e=aa*(a[i][j]+psi1*a[i1][j]+
							((i2 <= n) ? psi2*a[i2][j] : 0.0));
					a[i][j] -= e;
					a[i1][j] -= psi1*e;
					if (i2 <= n) a[i2][j] -= psi2*e;
				}
				for (j=q; j<=((i2 <= n) ? i2 : n); j++) {
					e=aa*(a[j][i]+psi1*a[j][i1]+
							((i2 <= n) ? psi2*a[j][i2] : 0.0));
					a[j][i] -= e;
					a[j][i1] -= psi1*e;
					if (i2 <= n) a[j][i2] -= psi2*e;
				}
				if (i2 <= n1) {
					i3=i+3;
					e=aa*psi2*a[i3][i2];
					a[i3][i] = -e;
					a[i3][i1] = -psi1*e;
					a[i3][i2] -= psi2*e;
				}
			}
		}
	} while (n > 0);
	em[3]=w;
	em[5]=count;
	return n;
}

