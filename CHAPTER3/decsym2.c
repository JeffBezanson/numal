#include "../real.h"


void decsym2(real_t **a, int n, real_t tol,
				int aux[], int p[], real_t detaux[])
{
	void elmrow(int, int, int, int, real_t **, real_t **, real_t);
	void ichrow(int, int, int, int, real_t **);
	void ichrowcol(int, int, int, int, real_t **);
	int i,j,m,ip1,ip2,onebyone,sym;
	real_t det,s,t,alpha,lambda,sigma,aii,aip1,aip1i,temp;

	aux[3]=aux[4]=0;
	sym=1;
	i=0;
	while (sym && (i < n)) {
		i++;
		j=i;
		while (sym && (j < n)) {
			j++;
			sym = sym && (a[i][j] == a[j][i]);
		}
	}
	if (sym)
		aux[2]=1;
	else {
		aux[2]=0;
		aux[5]=n;
		return;
	}
	alpha=(1.0+sqrt(17.0))/8.0;
	p[n]=n;
	i=1;
	while (i < n) {
		ip1=i+1;
		ip2=i+2;
		aii=fabs(a[i][i]);
		p[i]=i;
		lambda=fabs(a[i][ip1]);
		j=ip1;
		for (m=ip2; m<=n; m++)
			if (fabs(a[i][m]) > lambda) {
				j=m;
				lambda=fabs(a[i][m]);
			}
		t=alpha*lambda;
		onebyone=1;
		if (aii < t) {
			sigma=lambda;
			for (m=ip1; m<=j-1; m++)
				if (fabs(a[m][j]) > sigma) sigma=fabs(a[m][j]);
			for (m=j+1; m<=n; m++)
				if (fabs(a[j][m]) > sigma) sigma=fabs(a[j][m]);
			if (sigma*aii < lambda) {
				if (alpha*sigma < fabs(a[j][j])) {
					ichrow(j+1,n,i,j,a);
					ichrowcol(ip1,j-1,i,j,a);
					t=a[i][i];
					a[i][i]=a[j][j];
					a[j][j]=t;
					p[i]=j;
				} else {
					if (j > ip1) {
						ichrow(j+1,n,ip1,j,a);
						ichrowcol(ip2,j-1,ip1,j,a);
						t=a[i][i];
						a[i][i]=a[j][j];
						a[j][j]=t;
						t=a[i][j];
						a[i][j]=a[i][ip1];
						a[i][ip1]=t;
					}
					temp=a[i][ip1];
					det=a[i][i]*a[ip1][ip1]-temp*temp;
					aip1i=a[i][ip1]/det;
					aii=a[i][i]/det;
					aip1=a[ip1][ip1]/det;
					p[i]=j;
					p[ip1]=0;
					detaux[i]=1.0;
					detaux[ip1]=det;
					for (j=ip2; j<=n; j++) {
						s=aip1i*a[ip1][j]-aip1*a[i][j];
						t=aip1i*a[i][j]-aii*a[ip1][j];
						elmrow(j,n,j,i,a,a,s);
						elmrow(j,n,j,ip1,a,a,t);
						a[i][j]=s;
						a[ip1][j]=t;
					}
					aux[3]++;
					aux[4]++;
					i=ip2;
					onebyone=0;
				}
			}
		}
		if (onebyone) {
			if (tol < fabs(a[i][i])) {
				aii=a[i][i];
				detaux[i]=a[i][i];
				if (aii > 0.0)
					aux[3]++;
				else
					aux[4]++;
				for (j=ip1; j<=n; j++) {
					s = -a[i][j]/aii;
					elmrow(j,n,j,i,a,a,s);
					a[i][j]=s;
				}
			}
			i=ip1;
		}
	}
	if (i == n) {
		if (tol < fabs(a[n][n])) {
			if (a[n][n] > 0.0)
				aux[3]++;
			else
				aux[4]++;
		}
		detaux[n]=a[n][n];
	}
	aux[5]=n-aux[3]-aux[4];
}
